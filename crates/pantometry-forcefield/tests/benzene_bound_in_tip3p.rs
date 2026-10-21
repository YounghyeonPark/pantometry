//! **T4 lysozyme L99A and benzene (PDB 181L) in a box of TIP3P water**: the complex leg's system,
//! checked against closed forms and exact identities, and — ignored, in release — what each way
//! of holding the protein costs. Step W4 of the explicit-water track, phases 1 and 2: nothing here
//! runs the complex leg; `the_complex_leg_measured` is written and ignored, and
//! `the_cost_of_the_mobile_leg_measured` costs it.
//!
//! Default, unoptimised, on pockets cut out of 181L so that each takes seconds:
//!
//! - **the box is what it says**: neutral, the counter-ions one per charge, the lattice at its
//!   density, and every water the complex overlaps gone and no other, by brute force over 27
//!   images; each ion where the water farthest from everything was;
//! - **an ion is a point charge with Joung and Cheatham's Lennard-Jones term**, combined with a
//!   water oxygen by UFF's geometric rule;
//! - **a protein in a large box is the vacuum protein**: bonded terms and van der Waals the vacuum
//!   force field's, with the chain cut by every face and each atom wrapped on its own, and the
//!   electrostatics the vacuum's plus the closed-form terms of a charged cluster in a periodic box,
//!   `k_e [Q²ξ/2L + (2π/3V)(Q Σ q r² − μ²)]`, with what is left falling as `L⁻⁵`;
//! - **the end states are the systems they say**, with a counter-ion in the box and the restraint on;
//! - **`∂U/∂λ_r`, `∂U/∂λ_e`, `∂U/∂λ_v`** against central differences of the whole energy;
//! - **every state of `couplings` is its own `coupling`**, and **the decoupled reciprocal energy
//!   is the rest's own plus `λ_e` times the cross term** of three whole sums through the mesh;
//! - **each flexibility freezes what it says**, and its frozen atoms keep their bits;
//! - **constrained NVE falls as `h²`**, with the excess over four shrinking by four, on a pocket
//!   whose held bonds share atoms;
//! - **the Boresch release is its closed form**, at the anchors the rule chooses here;
//! - phase 2, **the leg's pieces**: its hydrogens relaxed in vacuum first, a stationary point in
//!   what was freed; phase 1's crevice water between Arg95 and Trp126 removed, by the rule for it;
//!   its box and its decoupling through a neighbour list, the same bits; its preparation the frozen
//!   melt and then the release, run as the stages say; and its schedule;
//! - after the leg, **benzene's coupling taken apart**: its van der Waals by protein UFF type, the
//!   waters and the ions, and its cross Coulomb by class, summing to the field's own coupling and
//!   to the whole evaluation's `U(coupled) − U(decoupled)` within their bounds; and **the refined
//!   window at λ_r 0.25** the one the insertion would run, to its header and its bits.
//! - **benzene's cross van der Waals reparametrised**: a variant at the current parameters
//!   changing nothing to the bit, Lorentz–Berthelot its formula for one pair, ff14SB's types
//!   the templates' for a handful of residues, each variant's `ΔU` two whole sums apart, frames
//!   round-tripping to the bit through their file, and the apportioning resuming to the same
//!   bits.
//! - **benzene's charges swapped**, QEq's for FreeSolv's AM1-BCC: a `.mol2` read with its charge
//!   column, FreeSolv's charges mapped by bonding and not by file order, a swap to benzene's own
//!   charges changing nothing to the bit, and a swap two whole evaluations apart through the mesh
//!   and through the classical sum.
//! - **variant 4's Hamiltonian staged**, `H(s) = H_UFF + s ΔU₄`: its forces the gradient of the
//!   frames' `ΔU₄`, the reciprocal sum's on an atom past the cutoff included; at `s` = 0 the field
//!   to the bit, its dynamics too; `H(s) − H(0)` the frames' `s ΔU₄` at every state; each state
//!   what its header says; and a staged run resuming to the same bits.
//!
//! Ignored, in release: the leg (`the_complex_leg_measured`, run once), its window at λ_r 0.25
//! (`the_restraint_interval_refined`), where its 3.85 kcal/mol comes from
//! (`where_the_binding_error_comes_from_measured`), and where the rest of it does, by Lorentz–
//! Berthelot, GAFF's benzene and ff14SB's protein on saved frames
//! (`the_rest_of_the_binding_error_apportioned_measured`), what AM1-BCC's charges change on
//! the same frames (`benzenes_charges_perturbed_measured`), and that one step checked by sampling
//! the states between (`variant_4_staged_measured`, `variant_4_staged_short_measured`, costed by
//! `the_staged_step_costed`).

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod amber;
mod freesolv;
mod gaff;
mod hydration;
mod protein;

use pantometry_forcefield::energy::Pair;
use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::free_energy::{bennett_chain, refine_schedule, window_seed};
use pantometry_forcefield::uff::{self, UffType, KCAL_PER_MOL};
use pantometry_forcefield::water::{self, WaterBox};
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Binding, Boresch, Complex, Component, Coupling, Element, Estimate,
    EwaldParameters, Flexibility, Ion, Lambda, MolecularDynamics, PeriodicBox, PeriodicDecoupling,
    PeriodicEnergy, PeriodicForceField, Protocol, SolvatedComplex, Solvation, Stage, System,
};
use pantometry_units::BOLTZMANN;
use std::sync::OnceLock;

const ANGSTROM: f64 = 1e-10;
const FS: f64 = 1e-15;
const EPS: f64 = f64::EPSILON;
/// W3's temperature, and the experiment's.
const KELVIN: f64 = 298.15;
/// The Wigner constant of the simple cubic lattice, as `tests/the_ewald_sum_against_closed_forms.rs`
/// checks it.
const XI: f64 = -2.837_297_479_480_62;

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn the_system() -> &'static System {
    static S: OnceLock<System> = OnceLock::new();
    S.get_or_init(protein::system)
}

/// 181L's pocket at 3 Å: four residues, 82 atoms with benzene, neutral.
fn pocket3() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| Binding::new(the_system(), 3.0 * ANGSTROM).expect("181L at 3 Å"))
}

/// 181L's pocket at 5 Å: thirteen residues, 245 atoms with benzene, +1 (Lys85).
fn pocket5() -> &'static Binding {
    static B: OnceLock<Binding> = OnceLock::new();
    B.get_or_init(|| Binding::new(the_system(), 5.0 * ANGSTROM).expect("181L at 5 Å"))
}

/// W4's solvation with a smaller margin and cutoff, for a pocket in a box of a few hundred waters.
fn small_solvation(margin: f64, cutoff: f64, seed: u64) -> Solvation {
    Solvation {
        margin: margin * ANGSTROM,
        cutoff: cutoff * ANGSTROM,
        seed,
        ..Solvation::w4()
    }
}

/// The 3 Å pocket in water: 3.5 Å margins, `r_c` = 6 Å.
fn small3(seed: u64) -> SolvatedComplex {
    SolvatedComplex::new(pocket3(), &small_solvation(3.5, 6.0, seed)).expect("the 3 Å box")
}

/// The 5 Å pocket in water with its chloride: 3.5 Å margins, `r_c` = 5.5 Å.
fn small5(seed: u64) -> SolvatedComplex {
    SolvatedComplex::new(pocket5(), &small_solvation(3.5, 5.5, seed)).expect("the 5 Å box")
}

/// The magnitudes the energy is a sum of, for its rounding.
fn magnitude(e: &PeriodicEnergy) -> f64 {
    let w = &e.ewald;
    e.bond.abs()
        + e.angle.abs()
        + e.torsion.abs()
        + e.inversion.abs()
        + e.van_der_waals.abs()
        + e.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs()
        + w.background.abs()
}

/// The largest force component, for a scale.
fn largest_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and the energy's rounding, at most `64 ε` of the sum of its parts'
/// magnitudes, divided by `h`. W1's and W3's.
fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    (d1, 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h)
}

/// The Boresch force constants of 3c-1 and 3c-3: 20 kcal mol⁻¹ Å⁻² and rad⁻² in `½ K (ξ − ξ₀)²`.
fn force_constants() -> [f64; 6] {
    let r = 20.0 * KCAL_PER_MOL / (ANGSTROM * ANGSTROM);
    let a = 20.0 * KCAL_PER_MOL;
    [r, a, a, a, a, a]
}

/// 3c-1's and 3c-3's anchor rule, the same code: receptor candidates `(N, CA, C)` and `(C, CA, N)`
/// of every residue of the binding, ligand candidates every path of three bonded heavy atoms, `r₀`
/// in [3, 8] Å, the four hinge angles furthest from straight.
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
    let heavy = |k: usize| b.elements()[k] != Element::H;
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

/// The box after `steps` of constrained BAOAB at 0.5 fs at `lambda` under `flexibility`, in a
/// 50 ps⁻¹ bath from 298 K velocities — W2's melt, because the lattice's waters start as close as
/// the 1.4 Å clearance to the complex and a 2 fs step from there turns a water past what SETTLE
/// can solve: a configuration constrained dynamics left.
fn moved(
    s: &SolvatedComplex,
    d: &PeriodicDecoupling,
    flexibility: Flexibility,
    lambda: Lambda,
    steps: usize,
    seed: u64,
) -> Vec<[f64; 3]> {
    let mut at = s.positions().to_vec();
    let mut md = s
        .dynamics(flexibility)
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 50e12,
            seed,
        })
        .thermalised(&at, KELVIN, seed);
    md.run(
        &AtLambda {
            hamiltonian: d,
            lambda,
        },
        &mut at,
        0.5 * FS,
        steps,
    );
    at
}

/// Benzene, the dictionary's: 181L's ligand atom for atom, with its own typing, checked equal to
/// the binding's.
fn benzene(b: &Binding) -> Component {
    let bnz = Component::from_ccd(protein::ccd("BNZ")).expect("BNZ");
    let lig = b.ligand_range();
    let elements: Vec<Element> = bnz.atoms().iter().map(|a| a.element).collect();
    assert_eq!(&b.elements()[lig.clone()], &elements[..]);
    assert_eq!(&b.types()[lig], &uff::assign(&bnz)[..]);
    bnz
}

// ---------------------------------------------------------------------------------------------
// The box

/// The waters and the ions of `s` as `solvation` built them, against brute force over the 27
/// nearest images of each pair: the lattice rebuilt and **its sites checked independently** — each
/// water's centre of mass at `(g + ½) L_a/k_a`, `g` its index read as `[x, y, z]` with `z` fastest,
/// to `16 ε L` — then every water with an atom within the clearance of a complex atom, or its
/// oxygen within the oxygen's clearance of two heavy atoms on opposite sides of it (each heavy
/// atom's nearest image, every pair's `(O − a)·(O − b)` below zero), found and removed; and, when
/// `every_branch`, **each rule must have removed a water the other would not**, and one water must
/// have an oxygen that close to a heavy atom and be kept, so that every branch is seen; **each ion
/// in turn on
/// the oxygen of the kept water farthest from every complex atom and every ion already placed**,
/// the first on a tie, that water gone; and the survivors in order, to the bit.
fn check_waters_and_ions(s: &SolvatedComplex, solvation: &Solvation, every_branch: bool) {
    let n = s.complex_atoms();
    let at = s.positions();
    let total = at.len();
    let ions = s.ions();
    let k = s.lattice();
    let lengths = s.cell().lengths();
    let lattice = WaterBox::lattice_box(k, solvation.density, solvation.seed);
    assert_eq!(lattice.cell().lengths(), lengths);
    let m = water::masses();
    let mass = m[0] + m[1] + m[2];
    for w in 0..lattice.count() {
        let g = [w / (k[1] * k[2]), (w / k[2]) % k[1], w % k[2]];
        let p = &lattice.positions()[3 * w..3 * w + 3];
        for a in 0..3 {
            let com = (m[0] * p[0][a] + m[1] * p[1][a] + m[2] * p[2][a]) / mass;
            let site = (g[a] as f64 + 0.5) * lengths[a] / k[a] as f64;
            assert!(
                (com - site).abs() <= 16.0 * EPS * lengths[a],
                "water {w} axis {a}: centre of mass {com:e} against the site {site:e}"
            );
        }
    }
    let nearest27 = |p: [f64; 3], others: &[[f64; 3]]| {
        let mut best = f64::INFINITY;
        for o in others {
            for i in -1..=1 {
                for j in -1..=1 {
                    for m in -1..=1 {
                        let image = [
                            o[0] + f64::from(i) * lengths[0],
                            o[1] + f64::from(j) * lengths[1],
                            o[2] + f64::from(m) * lengths[2],
                        ];
                        best = best.min(len(sub(p, image)));
                    }
                }
            }
        }
        best
    };
    let complex = &at[..n];
    let heavy: Vec<[f64; 3]> = (0..n)
        .filter(|&i| s.binding().elements()[i] != Element::H)
        .map(|i| at[i])
        .collect();
    // The nearest image of `o − h`, over the 27.
    let nearest_vector = |o: [f64; 3], h: [f64; 3]| {
        let mut best = [f64::INFINITY; 3];
        for i in -1..=1 {
            for j in -1..=1 {
                for m in -1..=1 {
                    let v = sub(
                        o,
                        [
                            h[0] + f64::from(i) * lengths[0],
                            h[1] + f64::from(j) * lengths[1],
                            h[2] + f64::from(m) * lengths[2],
                        ],
                    );
                    if dot(v, v) < dot(best, best) {
                        best = v;
                    }
                }
            }
        }
        best
    };
    let mut kept: Vec<usize> = Vec::new();
    let (mut by_clearance, mut between_only, mut close_kept) = (0, 0, 0);
    for w in 0..lattice.count() {
        let atoms = &lattice.positions()[3 * w..3 * w + 3];
        let clear = atoms
            .iter()
            .all(|&p| nearest27(p, complex) >= solvation.clearance);
        let close: Vec<[f64; 3]> = heavy
            .iter()
            .map(|&h| nearest_vector(atoms[0], h))
            .filter(|v| len(*v) < solvation.oxygen_clearance)
            .collect();
        let mut between = false;
        for a in 0..close.len() {
            for b in 0..a {
                between |= dot(close[a], close[b]) < 0.0;
            }
        }
        match (clear, between) {
            (true, false) => {
                kept.push(w);
                close_kept += usize::from(!close.is_empty());
            }
            (false, _) => by_clearance += 1,
            (true, true) => between_only += 1,
        }
    }
    println!(
        "removed: {by_clearance} by the clearance, {between_only} more between two heavy atoms; \
         kept {close_kept} with the oxygen that close to one side only"
    );
    if every_branch {
        assert!(
            by_clearance > 0 && between_only > 0 && close_kept > 0,
            "every branch of the rule is seen"
        );
    }
    assert_eq!(lattice.count() - kept.len(), s.overlapping());
    // Each ion in turn: the kept water whose oxygen is farthest from the complex and from every
    // ion placed before it.
    let mut placed: Vec<[f64; 3]> = Vec::new();
    let mut taken: Vec<usize> = Vec::new();
    for ion in 0..ions {
        let mut best: Option<(usize, f64)> = None;
        for &w in &kept {
            if taken.contains(&w) {
                continue;
            }
            let o = lattice.positions()[3 * w];
            let d = nearest27(o, complex).min(nearest27(o, &placed));
            let better = match best {
                None => true,
                Some((_, b)) => d > b,
            };
            if better {
                best = Some((w, d));
            }
        }
        let (w, d) = best.expect("a water for the ion");
        let site = lattice.positions()[3 * w];
        assert_eq!(at[total - ions + ion], site, "ion {ion}'s site");
        println!(
            "ion {ion}: {:.2} Å from the complex and the ions before it",
            d / ANGSTROM
        );
        placed.push(site);
        taken.push(w);
    }
    let expected: Vec<[f64; 3]> = kept
        .iter()
        .filter(|w| !taken.contains(w))
        .flat_map(|&w| lattice.positions()[3 * w..3 * w + 3].iter().copied())
        .collect();
    assert_eq!(
        &at[n..total - ions],
        &expected[..],
        "the waters kept, in order, to the bit"
    );
}

/// **Two ions are each where the water farthest from everything placed was**: the 5 Å pocket, +1,
/// neutralised by two half-charged chlorides — an ion of −0.5 e with the chloride's van der Waals —
/// so that the second ion's place depends on the first's. [`check_waters_and_ions`] holds both by
/// brute force, the distance to the ion already placed included, and the box is neutral. **The only
/// test with more than one ion**: with one, an update of the distances to the ions placed that did
/// nothing passed.
#[test]
fn two_ions_are_each_farthest_from_everything_placed() {
    let half = pantometry_forcefield::Ion {
        name: "Cl-/2",
        charge: -0.5,
        ..pantometry_forcefield::Ion::chloride()
    };
    let solvation = Solvation {
        counter_ion: half,
        ..small_solvation(3.5, 5.5, 0x2C1)
    };
    let s = SolvatedComplex::new(pocket5(), &solvation).unwrap();
    assert_eq!(s.ions(), 2);
    let q = s.field().charges();
    assert_eq!(&q[q.len() - 2..], &[-0.5, -0.5]);
    let sum: f64 = q.iter().sum();
    let scale: f64 = q.iter().map(|x| x.abs()).sum();
    assert!(sum.abs() <= q.len() as f64 * EPS * scale);
    check_waters_and_ions(&s, &solvation, false);
}

/// **The box is what it says**, on the 5 Å pocket, whose formal charge is +1:
///
/// - **neutral**: one chloride, the last atom, of charge −1, and every charge summing to zero
///   within `n ε Σ|q|`, the rounding of the sum, QEq's solve having hit its total;
/// - **at its density**: the lattice's waters over the box's volume is 33.00 nm⁻³ to `8 ε`, its edges
///   whole spacings, and each edge at least the complex's extent and both margins and less than a
///   spacing more, the complex's heavy atoms moved rigidly to the centre and each hydrogen on its
///   held bond's length;
/// - **no water overlaps the complex, and no other was removed**: the lattice rebuilt, every water
///   with an atom within the clearance of a complex atom, or its oxygen within the oxygen's
///   clearance of a heavy atom, found by brute force over the 27 nearest images of each pair, the
///   rest kept in order to the bit;
/// - **the lattice is where it should be, independently**: every lattice water's centre of mass at
///   `(g + ½) L_a/k_a`, written here from its index, to `16 ε L`;
/// - **the ion is where the water farthest from the complex was**, by the same brute force, and
///   that water is gone ([`check_waters_and_ions`]);
/// - the chloride's parameters are Joung and Cheatham's as two printed files carry them: σ
///   OpenMM's `amber14/tip3p.xml`'s 0.4477656957373345 nm as a literal, `x⁶ = 2σ⁶` to rounding,
///   `x/2` AMBER's `frcmod.ionsjc_tip3p`'s `R_min/2` of 2.513 Å and `D` its 0.035 591 0 kcal/mol,
///   each to its last printed digit.
#[test]
fn the_box_is_neutral_at_its_density_and_clear_of_the_complex() {
    let solvation = small_solvation(3.5, 5.5, 0xB0C);
    let s = SolvatedComplex::new(pocket5(), &solvation).unwrap();
    let n = s.complex_atoms();
    let (p, l) = s.binding().formal_charges();
    assert_eq!((p, l), (1, 0), "the 5 Å pocket is +1 and benzene neutral");
    assert_eq!(s.ions(), 1);
    assert_eq!(s.field().ion_atoms(), 1);
    let q = s.field().charges();
    let total = q.len();
    assert_eq!(total, n + 3 * s.waters() + 1);
    assert_eq!(q[total - 1], -1.0);
    let sum: f64 = q.iter().sum();
    let scale: f64 = q.iter().map(|x| x.abs()).sum();
    println!(
        "{} complex atoms, {} waters ({} overlapping removed, 1 for the ion), the net charge \
         {sum:+.3e} e against a bound {:.3e}",
        n,
        s.waters(),
        s.overlapping(),
        total as f64 * EPS * scale
    );
    assert!(sum.abs() <= total as f64 * EPS * scale);

    // The density and the box.
    let k = s.lattice();
    let cell = s.cell();
    let lengths = cell.lengths();
    let lattice_waters = (k[0] * k[1] * k[2]) as f64;
    let number = lattice_waters / cell.volume();
    println!(
        "lattice {k:?}, box {:.3} × {:.3} × {:.3} Å, {:.6e} waters per m³",
        lengths[0] / ANGSTROM,
        lengths[1] / ANGSTROM,
        lengths[2] / ANGSTROM,
        number
    );
    assert!((number / 33.0e27 - 1.0).abs() <= 8.0 * EPS);
    let spacing = lengths[0] / k[0] as f64;
    let at = s.positions();
    // The complex moved rigidly, its heavy atoms by one shift, and the box sized and centred on
    // the binding's extent; its hydrogens each on its held bond's length.
    let given = s.binding().positions();
    let heavy: Vec<usize> = (0..n)
        .filter(|&i| s.binding().elements()[i] != Element::H)
        .collect();
    let shift = sub(at[heavy[0]], given[heavy[0]]);
    for &i in &heavy {
        let d = sub(sub(at[i], given[i]), shift);
        assert!(
            len(d) <= 4.0 * EPS * lengths[0],
            "heavy atom {i} moved rigidly"
        );
    }
    for a in 0..3 {
        assert!((lengths[a] / k[a] as f64 - spacing).abs() <= 4.0 * EPS * spacing);
        let lo = given.iter().map(|p| p[a]).fold(f64::INFINITY, f64::min);
        let hi = given.iter().map(|p| p[a]).fold(f64::NEG_INFINITY, f64::max);
        let need = hi - lo + 2.0 * solvation.margin;
        assert!(
            lengths[a] >= need && lengths[a] < need + spacing,
            "axis {a}"
        );
        assert!((0.5 * (lo + hi) + shift[a] - 0.5 * lengths[a]).abs() <= 1e-14 * lengths[a]);
    }
    let shake = s.shake();
    let mut moved = 0.0f64;
    for (&[i, j], &d0) in shake.bonds().iter().zip(shake.lengths()) {
        assert!(
            (len(sub(at[i], at[j])) / d0 - 1.0).abs() <= 2e-12,
            "bond {i}–{j} held"
        );
        let h = if s.binding().elements()[i] == Element::H {
            i
        } else {
            j
        };
        moved = moved.max(len(sub(sub(at[h], given[h]), shift)));
    }
    println!(
        "{} held bonds at their lengths; a hydrogen moved at most {:.4} Å onto its",
        shake.len(),
        moved / ANGSTROM
    );

    check_waters_and_ions(&s, &solvation, true);

    // The chloride's parameters.
    let ion = pantometry_forcefield::Ion::chloride();
    let (x, sigma) = (
        ion.vdw_distance,
        pantometry_forcefield::solvated::CHLORIDE_SIGMA,
    );
    let (x3, s3) = (x * x * x, sigma * sigma * sigma);
    assert!((x3 * x3 / (2.0 * s3 * s3) - 1.0).abs() <= 16.0 * EPS);
    // σ as OpenMM's `amber14/tip3p.xml` prints it, 0.4477656957373345 nm, and `x/2` as AMBER's
    // `frcmod.ionsjc_tip3p` prints Joung and Cheatham's `R_min/2`, 2.513 Å, to its last digit
    // (ambermini, commit f7c421c4, SHA-256 eb7a1a59…8a55: `Cl-      2.513     0.0355910`).
    assert_eq!(sigma, 4.477_656_957_373_345e-10);
    assert!(
        (0.5 * x / ANGSTROM - 2.513).abs() < 5e-4,
        "R_min/2 {}",
        0.5 * x / ANGSTROM
    );
    assert!((kcal(ion.vdw_energy) - 0.035_591_0).abs() < 5e-8);
    assert_eq!(ion.charge, -1.0);
}

/// **An ion is a point charge with a Lennard-Jones term**: one chloride beside one TIP3P water in
/// a 30 Å cube, `r_c` = 12 Å, at five distances along a line from the oxygen. The van der Waals
/// energy is `4ε_ij[(σ_ij/r)¹² − (σ_ij/r)⁶]` with `σ_ij = √(σ_Cl σ_O)` and `ε_ij = √(ε_Cl ε_O)` —
/// UFF's geometric rule, written in Lennard-Jones's terms from Joung and Cheatham's σ and ε and
/// Jorgensen's — to `10⁻¹³` of its two parts' magnitudes; a water hydrogen has no well, so the
/// oxygen's is the only pair. The ion carries −1 and is excluded from nothing.
#[test]
fn an_ion_is_a_charge_and_a_lennard_jones_term() {
    let cell = PeriodicBox::cubic(30.0 * ANGSTROM);
    let p = EwaldParameters::for_accuracy(&cell, 12.0 * ANGSTROM, 1e-6);
    let field =
        PeriodicForceField::tip3p(cell, p, 1).with_ions(&[pantometry_forcefield::Ion::chloride()]);
    assert_eq!(field.ion_atoms(), 1);
    assert_eq!(field.charges()[3], -1.0);
    assert_eq!(field.excluded_pairs(), &[[0, 1], [0, 2], [1, 2]]);
    let sigma = (pantometry_forcefield::solvated::CHLORIDE_SIGMA * water::sigma()).sqrt();
    let eps_cl = pantometry_forcefield::solvated::CHLORIDE_EPSILON_KJ * 1e3 / uff::AVOGADRO;
    let epsilon = (eps_cl * water::epsilon()).sqrt();
    let centre = [15.0 * ANGSTROM; 3];
    let body = water::body_frame();
    let o = body[0];
    for r in [3.0, 3.5, 4.2, 5.0, 7.0] {
        let mut at: Vec<[f64; 3]> = body
            .iter()
            .map(|b| [0, 1, 2].map(|k| centre[k] + b[k] - o[k]))
            .collect();
        // Along +y, away from the hydrogens, which the body frame puts at −y.
        at.push([centre[0], centre[1] + r * ANGSTROM, centre[2]]);
        let e = field.energy(&at);
        let x = sigma / (r * ANGSTROM);
        let x6 = x * x * x * x * x * x;
        let want = 4.0 * epsilon * (x6 * x6 - x6);
        let scale = 4.0 * epsilon * (x6 * x6 + x6);
        println!(
            "r {r} Å: van der Waals {:+.10e} kcal/mol, Lennard-Jones {:+.10e}",
            kcal(e.van_der_waals),
            kcal(want)
        );
        assert!((e.van_der_waals - want).abs() <= 1e-13 * scale, "r {r}");
    }
}

// ---------------------------------------------------------------------------------------------
// A protein in a large box

/// **A protein in a large box is the vacuum protein**: the 5 Å pocket and benzene, +1, alone in
/// cubes of 48, 64 and 80 Å with `r_c = L/2`, past every pair, so that no atom meets an image.
///
/// - **As given**, every atom inside the box: the bonded terms are the vacuum force field's to the
///   bit, the van der Waals to `10⁻¹³`.
/// - **Moved so that the chain crosses all three faces, and each atom wrapped into the box on its
///   own**, so that the walk must rebuild every residue, every peptide bond and benzene: the same
///   bonded terms and van der Waals to `10⁻¹³` of their magnitudes together — the rounding of a
///   coordinate moves every term by its gradient, and a nearly planar inversion's energy is small
///   beside its gradient.
/// - **The electrostatics** (the wrapped positions) are the vacuum's plus the two closed-form terms
///   of a charged cluster with Ewald's neutralising background:
///   `k_e [Q²ξ/(2L) + (2π/3V)(Q Σ q_i r_i² − μ²)]`. The first is a lone charge's Wigner energy; the
///   second is the quadratic part of the periodic potential, `ψ(r) = 1/r + ξ/L + 2πr²/3V + O(r⁴/L⁵)`,
///   summed over pairs (translation-invariant: a shift adds `2Q a·μ + Q²a²` to both parts).
///   **What is left falls as `L⁻⁵`** — the first cubic harmonic, `r⁴` — energy and forces, each
///   within a quarter from one box to the next, W1's band, and the closed-form terms are at least
///   four times what is left; **and the departure from `L⁻⁵` shrinks as the next harmonic's
///   `L⁻⁷` says**, by 0.463 from one pair of boxes to the next. The force of the quadratic term is
///   `−(4πk_e q_i/3V)(Q r_i − μ)`.
///
/// **The shrink is the guard, not the quarter.** Measured departures from `L⁻⁵` are 5.2% and 2.2%,
/// so the ±25% band is loose: the Ewald background scaled by 1.01 passed it (3.542 and 2.414) and
/// was caught by the shrink alone, 1.310. The shrink catches the background off by 1e-3 of itself
/// (0.751) and not by 3e-4 (0.541, inside 0.463 ± 0.15). The background moves no force; the forces'
/// remainder is held to the quarter band and to a quarter of the closed form's force.
#[test]
fn a_protein_in_a_large_box_is_the_vacuum_protein() {
    let b = pocket5();
    let c = b.component();
    let types = b.types();
    let q = b.charges().to_vec();
    let vacuum = b.force_field();
    let n = q.len();
    let q_total: f64 = q.iter().sum();
    assert!((q_total - 1.0).abs() < 1e-12);
    let mut rows = Vec::new();
    for side in [48.0, 64.0, 80.0] {
        let l = side * ANGSTROM;
        let cell = PeriodicBox::cubic(l);
        let p = EwaldParameters::for_accuracy(&cell, 0.5 * l, 1e-10);
        let field = PeriodicForceField::new(c, types, cell, p)
            .expect("the pocket")
            .with_charges(q.clone());
        // As given, moved into the box whole.
        let centre = pantometry_forcefield::complex::centroid(b.positions());
        let inside: Vec<[f64; 3]> = b
            .positions()
            .iter()
            .map(|p| [0, 1, 2].map(|k| p[k] - centre[k] + 0.5 * l))
            .collect();
        let v = vacuum.evaluate(&inside);
        let e = field.evaluate(&inside);
        assert_eq!(e.energy.bond, v.energy.bond);
        assert_eq!(e.energy.angle, v.energy.angle);
        assert_eq!(e.energy.torsion, v.energy.torsion);
        assert_eq!(e.energy.inversion, v.energy.inversion);
        assert!(
            (e.energy.van_der_waals - v.energy.van_der_waals).abs()
                <= 1e-13 * v.energy.van_der_waals.abs()
        );
        // Across every face: the centroid at a corner, each atom wrapped on its own.
        let wrapped: Vec<[f64; 3]> = b
            .positions()
            .iter()
            .map(|p| cell.wrap([0, 1, 2].map(|k| p[k] - centre[k] + 0.3 * ANGSTROM)))
            .collect();
        let unwrapped: Vec<[f64; 3]> = b
            .positions()
            .iter()
            .map(|p| [0, 1, 2].map(|k| p[k] - centre[k] + 0.3 * ANGSTROM))
            .collect();
        let crossing = (0..3)
            .map(|a| unwrapped.iter().filter(|p| p[a] < 0.0).count())
            .collect::<Vec<_>>();
        assert!(crossing.iter().all(|&k| k > 0 && k < n), "{crossing:?}");
        let v = vacuum.evaluate(&unwrapped);
        let e = field.evaluate(&wrapped);
        let terms = [
            ("bond", e.energy.bond, v.energy.bond),
            ("angle", e.energy.angle, v.energy.angle),
            ("torsion", e.energy.torsion, v.energy.torsion),
            ("inversion", e.energy.inversion, v.energy.inversion),
            (
                "van der Waals",
                e.energy.van_der_waals,
                v.energy.van_der_waals,
            ),
        ];
        // Wrapping and the walk move each coordinate by a rounding, which moves every term by
        // its gradient times that: a scale of the terms together, not of each.
        let together: f64 = terms.iter().map(|t| t.2.abs()).sum();
        let mut worst = 0.0f64;
        for (name, a, z) in terms {
            worst = worst.max((a - z).abs() / together);
            assert!(
                (a - z).abs() <= 1e-13 * together,
                "{name}: {a:e} against {z:e}"
            );
        }
        println!(
            "{side} Å, wrapped: the bonded terms and van der Waals off by {worst:.2e} of their sum"
        );
        // The electrostatics' closed-form terms, at the whole positions.
        let volume = cell.volume();
        let mu = [0, 1, 2].map(|a| (0..n).map(|i| q[i] * unwrapped[i][a]).sum::<f64>());
        let second: f64 = (0..n).map(|i| q[i] * dot(unwrapped[i], unwrapped[i])).sum();
        let wigner = COULOMB * q_total * q_total * XI / (2.0 * l);
        let quadratic = 2.0 * std::f64::consts::PI * COULOMB / (3.0 * volume)
            * (q_total * second - dot(mu, mu));
        let left = e.energy.electrostatic - v.energy.electrostatic - wigner - quadratic;
        let strength = 4.0 * std::f64::consts::PI * COULOMB / (3.0 * volume);
        // The forces less the vacuum's and the closed form's: what is left is the `r⁴/L⁵`
        // harmonic's, held below by its rate and by being under a quarter of the closed form's.
        let (mut off, mut closed) = (0.0f64, 0.0f64);
        let vf = &v.forces;
        let fscale = largest_force(vf);
        for i in 0..n {
            for a in 0..3 {
                let term = -strength * q[i] * (q_total * unwrapped[i][a] - mu[a]);
                off = off.max((e.forces[i][a] - vf[i][a] - term).abs());
                closed = closed.max(term.abs());
            }
        }
        println!(
            "{side} Å: Ewald − vacuum {:+.6e} kcal/mol; Wigner {:+.6e}, quadratic {:+.6e}, left \
             {:+.3e}; forces: closed form's largest {:.3e} N, left {:.3e}, of {:.3e}",
            kcal(e.energy.electrostatic - v.energy.electrostatic),
            kcal(wigner),
            kcal(quadratic),
            kcal(left),
            closed,
            off,
            fscale
        );
        rows.push((side, left, wigner.abs() + quadratic.abs(), off, closed));
    }
    let mut departures = Vec::new();
    for w in rows.windows(2) {
        let (a, b) = (w[0], w[1]);
        let r = b.0 / a.0;
        let expected = r * r * r * r * r;
        let energy = a.1 / b.1;
        let force = a.3 / b.3;
        departures.push(energy / expected - 1.0);
        println!(
            "{} → {} Å: left falls by {energy:.3}, forces' by {force:.3}, L⁻⁵ {expected:.3}",
            a.0, b.0
        );
        assert!((energy / expected - 1.0).abs() < 0.25, "energy {energy}");
        assert!((force / expected - 1.0).abs() < 0.25, "forces {force}");
        assert!(a.1.abs() < 0.25 * a.2 && a.3 < 0.25 * a.4);
    }
    // And the departure from `L⁻⁵` is the next harmonic's, `r⁶/L⁷`: for `A/L⁵ + B/L⁷` the ratio's
    // departure is `(B/A)(1/a² − 1/b²)`, so the second's over the first's is
    // `(1/64² − 1/80²)/(1/48² − 1/64²)` = 0.4629, whatever `B/A`. Held to ±0.15: the harmonic after
    // it, `r⁸/L⁹`, is a further `(a/L)²` smaller, and the measured value is printed.
    let predicted =
        (1.0 / (64.0 * 64.0) - 1.0 / (80.0 * 80.0)) / (1.0 / (48.0 * 48.0) - 1.0 / (64.0 * 64.0));
    let shrink = departures[1] / departures[0];
    println!("the departure from L⁻⁵ shrinks by {shrink:.3}, predicted {predicted:.4}");
    assert!((shrink - predicted).abs() < 0.15, "shrink {shrink}");
}

// ---------------------------------------------------------------------------------------------
// The decoupling in the box

/// **The end states are the systems they say**, on the 5 Å pocket with its chloride, at a
/// configuration constrained dynamics left with everything mobile, held bonds sharing atoms:
///
/// - **coupled** (`λ_r` = 0): the solvated field less benzene alone in the same box plus benzene in
///   vacuum;
/// - **decoupled and restrained** (`λ_r` = 1, `λ_e` = `λ_v` = 0): the pocket, the waters and the ion
///   alone in the box — a force field built without benzene — plus benzene in vacuum plus the
///   Boresch restraint's energy, `Boresch::energy`;
///
/// the energy to `10⁻¹²` of the parts and every force, each atom's from the system it belongs to
/// and the restraint's, to `10⁻¹²` of the largest. The mesh and the Ewald parameters are the
/// solvated field's in every system.
///
/// **What this cannot see**: the coupled expectation is built from `s.field()` itself, so an error
/// in assembling the solvated field — a wrong charge, a wrong well, an ion's parameters — is on
/// both sides and cancels; and the decoupled one builds the rest with the same `with_ions`. Those
/// are held elsewhere: the ion by `an_ion_is_a_charge_and_a_lennard_jones_term`, the complex in a
/// box by `a_protein_in_a_large_box_is_the_vacuum_protein`, the box by
/// `the_box_is_neutral_at_its_density_and_clear_of_the_complex`.
#[test]
fn the_end_states_are_the_systems_they_say() {
    let s = small5(0xE5D);
    let b = s.binding();
    let n0 = b.pocket_len();
    let lig = s.ligand();
    let restraint = anchors(b, &s.positions()[..s.complex_atoms()]);
    let d = s.decoupling().with_restraint(restraint);
    let at = moved(
        &s,
        &d,
        Flexibility::Mobile,
        Lambda::new(0.5, 0.5, 0.7),
        6,
        3,
    );
    let n = at.len();
    let cell = s.cell();
    let params = s.field().ewald().parameters();
    let mesh = s.field().ewald().mesh().expect("the mesh").parameters();
    let bnz = benzene(b);
    let vacuum = b.ligand_force_field().evaluate(&at[lig.clone()]);
    let alone = PeriodicForceField::new(&bnz, &b.types()[lig.clone()], cell, params)
        .unwrap()
        .with_charges(b.charges()[lig.clone()].to_vec())
        .with_mesh(mesh)
        .energy(&at[lig.clone()]);
    let whole = s.field().energy(&at);
    let scale = magnitude(&whole);
    let mut f = vec![[0.0; 3]; n];
    let coupled = d.energy_and_forces(&at, Lambda::COUPLED, &mut f);
    let expected = whole.total - alone.total + vacuum.energy.total;
    println!(
        "coupled {coupled:.10e} J against {expected:.10e}; benzene's own images {:+.4e} kcal/mol",
        kcal(alone.total - vacuum.energy.total)
    );
    assert!((coupled - expected).abs() <= 1e-12 * scale);

    // The rest alone: the pocket, the waters and the ion, without benzene.
    let rest = PeriodicForceField::solvated(
        b.pocket_component(),
        &b.types()[..n0],
        s.waters(),
        cell,
        params,
    )
    .unwrap()
    .with_solute_charges(b.charges()[..n0].to_vec())
    .with_ions(&[pantometry_forcefield::Ion::chloride()])
    .with_mesh(mesh);
    let rest_at: Vec<[f64; 3]> = (0..n).filter(|i| !lig.contains(i)).map(|i| at[i]).collect();
    let r = rest.evaluate(&rest_at);
    let decoupled = d.energy_and_forces(&at, Lambda::new(1.0, 0.0, 0.0), &mut f);
    let u_b = restraint.energy(&at);
    let expected = r.energy.total + vacuum.energy.total + u_b;
    println!(
        "decoupled {decoupled:.10e} J against {expected:.10e}; U_B {:.4} kcal/mol",
        kcal(u_b)
    );
    assert!(u_b > 1e-3 * KCAL_PER_MOL, "the restraint is stretched");
    assert!((decoupled - expected).abs() <= 1e-12 * scale);
    let mut want = vec![[0.0; 3]; n];
    let mut k = 0;
    for i in 0..n {
        if lig.contains(&i) {
            want[i] = vacuum.forces[i - lig.start];
        } else {
            want[i] = r.forces[k];
            k += 1;
        }
    }
    restraint.add_forces(&at, &mut want, 1.0);
    let fscale = largest_force(&f);
    for i in 0..n {
        for a in 0..3 {
            assert!(
                (f[i][a] - want[i][a]).abs() <= 1e-12 * fscale,
                "atom {i} axis {a}: {:e} against {:e}",
                f[i][a],
                want[i][a]
            );
        }
    }
}

/// **`∂U/∂λ_r`, `∂U/∂λ_e` and `∂U/∂λ_v` are the energy's derivatives** on the 3 Å pocket in water
/// with the restraint, at configurations constrained dynamics left at four states — the restraint
/// a quarter on, the charges part-way, the van der Waals part-way, and near the decoupled end
/// where the soft core is softest — against central differences of the whole energy, each where its
/// difference stays inside [0, 1], each to the tolerance it earns, and each able to see a part in
/// a thousand. **A quarter, not a half**: at `λ_r = ½` a restraint scaled by `λ_r²` has the
/// derivative `2λ_r U_B = U_B`, which the linear one has, and passed there.
#[test]
fn the_lambda_derivatives_are_the_whole_energys() {
    let s = small3(0xDE1);
    let restraint = anchors(s.binding(), &s.positions()[..s.complex_atoms()]);
    let d = s.decoupling().with_restraint(restraint);
    let mut checked = [0usize; 3];
    for (lambda, steps) in [
        (Lambda::new(0.25, 1.0, 1.0), 8),
        (Lambda::new(1.0, 0.4, 1.0), 8),
        (Lambda::new(1.0, 0.0, 0.5), 8),
        (Lambda::new(1.0, 0.0, 0.1), 12),
    ] {
        let at = moved(&s, &d, Flexibility::Mobile, lambda, steps, 11);
        let c = d.coupling(&at, lambda);
        let scale = magnitude(&s.field().energy(&at));
        for (k, h) in [(0usize, 1e-3), (1, 1e-3), (2, 1e-4)] {
            let at_k = lambda.components()[k];
            if at_k < 2.0 * h || at_k > 1.0 - 2.0 * h {
                continue;
            }
            let energy = |e: f64| {
                let mut l = lambda;
                match k {
                    0 => l.restraint += e,
                    1 => l.electrostatics += e,
                    _ => l.van_der_waals += e,
                }
                let mut f = vec![[0.0; 3]; at.len()];
                (d.energy_and_forces(&at, l, &mut f), scale)
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
            checked[k] += 1;
        }
    }
    assert!(
        checked.iter().all(|&c| c > 0),
        "every λ was checked: {checked:?}"
    );
}

/// **Every state of `couplings` is its own `coupling`, to the bit**, with the restraint, at a
/// configuration dynamics left; the couplings' differences are the whole energy's to `64 ε` of its
/// parts. **And the decoupled reciprocal energy itself is the rest's own whole sum plus `λ_e`
/// times the cross term of three whole sums**, through the mesh at `λ_e` = 1, 0.4 and 0: which a
/// set of checks that are all differences or symmetric in the two sets cannot tell from the rest
/// and the group swapped (the particle-mesh review's lesson).
#[test]
fn every_state_of_couplings_is_its_own_coupling() {
    let s = small3(0xC0B);
    let restraint = anchors(s.binding(), &s.positions()[..s.complex_atoms()]);
    let d = s.decoupling().with_restraint(restraint);
    let at = moved(
        &s,
        &d,
        Flexibility::Mobile,
        Lambda::new(1.0, 0.0, 0.6),
        8,
        13,
    );
    let states = grid();
    let all = d.couplings(&at, &states);
    for (k, l) in states.iter().enumerate() {
        let one = d.coupling(&at, *l);
        assert_eq!(one.energy.to_bits(), all[k].energy.to_bits(), "{l:?}");
        assert_eq!(
            one.gradient.map(f64::to_bits),
            all[k].gradient.map(f64::to_bits)
        );
    }
    let some = d.couplings(&at, &[states[20], states[3], states[20]]);
    assert_eq!(some[0].energy.to_bits(), all[20].energy.to_bits());
    assert_eq!(some[1].energy.to_bits(), all[3].energy.to_bits());
    let scale = magnitude(&s.field().energy(&at));
    let mut f = vec![[0.0; 3]; at.len()];
    let base = d.energy_and_forces(&at, states[0], &mut f);
    for (k, l) in states.iter().enumerate().skip(1) {
        let whole = d.energy_and_forces(&at, *l, &mut f) - base;
        let part = all[k].energy - all[0].energy;
        assert!(
            (whole - part).abs() <= 64.0 * EPS * scale,
            "{l:?}: {whole:e} against {part:e}"
        );
    }
    // The reciprocal energy itself.
    let group = s.ligand_mask();
    let reciprocal = |keep: &dyn Fn(usize) -> bool| {
        let q: Vec<f64> = (0..group.len())
            .map(|i| if keep(i) { s.field().charges()[i] } else { 0.0 })
            .collect();
        s.field()
            .clone()
            .with_charges(q)
            .energy(&at)
            .ewald
            .reciprocal
    };
    let whole = reciprocal(&|_| true);
    let rest = reciprocal(&|i| !group[i]);
    let own = reciprocal(&|i| group[i]);
    for e in [1.0, 0.4, 0.0] {
        let got = d.energy(&at, Lambda::new(1.0, e, 1.0)).ewald.reciprocal;
        let want = rest + e * (whole - rest - own);
        println!(
            "λ_e {e}: reciprocal {:+.10e} against {:+.10e} kcal/mol",
            kcal(got),
            kcal(want)
        );
        assert!((got - want).abs() <= 64.0 * EPS * scale);
    }
    // The cross term is a thousand times what the comparison allows, so a reciprocal energy
    // without it, or with the two sets swapped, is seen.
    assert!(
        (whole - rest - own).abs() > 1e3 * 64.0 * EPS * scale,
        "a cross term to see"
    );
}

/// **Each flexibility freezes what it says**: on the 3 Å pocket, `Frozen` moves benzene and no
/// protein atom, `Zone(r)` exactly `Complex::zone`'s residues, `Mobile` everything; the water and
/// the ion always move. After ten steps every frozen atom has its bits and every mobile one has
/// moved; the held bonds are every bond to a hydrogen of the complex, and the degrees of freedom
/// are `3 N_free − 3 N_water − N_held`, a held bond frozen at both ends counting none.
#[test]
fn each_flexibility_freezes_what_it_says() {
    let s = small3(0xF1E);
    let b = s.binding();
    let n = s.positions().len();
    let lig = s.ligand();
    let shake = s.shake();
    let hydrogen_bonds = b
        .component()
        .bonds()
        .iter()
        .filter(|bond| bond.atoms.iter().any(|&i| b.elements()[i] == Element::H))
        .count();
    assert_eq!(shake.len(), hydrogen_bonds);
    // A radius between the two nearest residues' nearest heavy atoms, so that the zone frees one
    // residue and freezes the rest.
    let crystal = b.crystal_positions();
    let mut nearest: Vec<f64> = b
        .residue_spans()
        .iter()
        .map(|span| {
            span.clone()
                .filter(|&k| b.elements()[k] != Element::H)
                .flat_map(|k| lig.clone().map(move |l| len(sub(crystal[k], crystal[l]))))
                .fold(f64::INFINITY, f64::min)
        })
        .collect();
    nearest.sort_by(f64::total_cmp);
    let radius = 0.5 * (nearest[0] + nearest[1]);
    let zone = Complex::zone(b, radius);
    let freed = zone[..b.pocket_len()].iter().filter(|&&m| m).count();
    assert!(
        freed > 0 && freed < b.pocket_len(),
        "the zone frees some and not all"
    );
    for flexibility in [
        Flexibility::Frozen,
        Flexibility::Zone(radius),
        Flexibility::Mobile,
    ] {
        let mobile = s.mobile(flexibility);
        for i in 0..n {
            let want = if i >= s.complex_atoms() || lig.contains(&i) {
                true
            } else {
                match flexibility {
                    Flexibility::Frozen => false,
                    Flexibility::Zone(_) => zone[i],
                    Flexibility::Mobile => true,
                }
            };
            assert_eq!(mobile[i], want, "{flexibility:?} atom {i}");
        }
        let d = s.decoupling();
        let at = moved(&s, &d, flexibility, Lambda::COUPLED, 10, 7);
        for i in 0..n {
            if mobile[i] {
                assert_ne!(at[i], s.positions()[i], "{flexibility:?}: atom {i} moved");
            } else {
                assert_eq!(
                    at[i],
                    s.positions()[i],
                    "{flexibility:?}: atom {i} kept its bits"
                );
            }
        }
        let md = s.dynamics(flexibility).with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 1e12,
            seed: 1,
        });
        let free = mobile.iter().filter(|&&m| m).count();
        let held = shake
            .bonds()
            .iter()
            .filter(|&&[i, j]| mobile[i] || mobile[j])
            .count();
        assert_eq!(md.degrees_of_freedom(), 3 * free - 3 * s.waters() - held);
        println!(
            "{flexibility:?}: {free} free atoms, {held} held bonds, {} degrees of freedom",
            md.degrees_of_freedom()
        );
    }
}

// ---------------------------------------------------------------------------------------------
// Constrained NVE

/// The 3 Å pocket, benzene and the eight waters nearest benzene in the small box that are 2.5 Å
/// clear of the complex, moved into a box
/// 52 Å across with `r_c` = 26 Å: no pair reaches the cutoff and no image comes inside it, so the
/// potential is smooth and NVE's energy error is the integrator's. The held bonds are the
/// pocket's and benzene's, which share atoms in every methyl.
fn cluster(seed: u64) -> (PeriodicForceField, Vec<[f64; 3]>, MolecularDynamics, usize) {
    let s = small3(seed);
    let b = s.binding();
    let n = s.complex_atoms();
    let centre = pantometry_forcefield::complex::centroid(&s.positions()[s.ligand()]);
    let mut waters: Vec<(f64, [usize; 3])> = s
        .field()
        .rigid_waters()
        .iter()
        .map(|&w| (len(sub(s.positions()[w[0]], centre)), w))
        .collect();
    waters.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Of those, only waters with every atom at least 2.5 Å from every complex atom: the lattice's
    // own placement, 1.4 Å from the complex at closest, is a clash a 2 fs step cannot start from.
    let clear = |w: &[usize; 3]| {
        w.iter().all(|&i| {
            s.positions()[..n]
                .iter()
                .all(|&p| len(sub(s.positions()[i], p)) >= 2.5 * ANGSTROM)
        })
    };
    let kept: Vec<[usize; 3]> = waters.iter().map(|w| w.1).filter(clear).take(8).collect();
    let side = 52.0 * ANGSTROM;
    let big = PeriodicBox::cubic(side);
    let p = EwaldParameters::for_accuracy(&big, 26.0 * ANGSTROM, 1e-5);
    let field = PeriodicForceField::solvated(b.component(), b.types(), kept.len(), big, p)
        .expect("the cluster")
        .with_solute_charges(b.charges().to_vec());
    let shift = sub([0.5 * side; 3], centre);
    let moved = |p: [f64; 3]| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]];
    let mut at: Vec<[f64; 3]> = s.positions()[..n].iter().map(|&p| moved(p)).collect();
    let mut masses = s.masses()[..n].to_vec();
    for w in &kept {
        for &i in w {
            at.push(moved(s.positions()[i]));
            masses.push(s.masses()[i]);
        }
    }
    let span = at
        .iter()
        .flat_map(|p| at.iter().map(move |q| len(sub(*p, *q))))
        .fold(0.0f64, f64::max);
    assert!(span < 26.0 * ANGSTROM - 4.0 * ANGSTROM, "span {span:e}");
    let md = MolecularDynamics::new(masses)
        .with_constraints(pantometry_forcefield::Settle::tip3p(
            field.rigid_waters().to_vec(),
        ))
        .with_bond_constraints(s.shake().clone());
    (field, at, md, kept.len())
}

/// The RMS departure of `kinetic + potential` from its start over `steps` of `dt`, NVE.
fn nve_rms(
    field: &PeriodicForceField,
    start: &[[f64; 3]],
    md: &MolecularDynamics,
    dt: f64,
    steps: usize,
    seed: u64,
) -> f64 {
    let mut at = start.to_vec();
    let mut md = md.clone().thermalised(&at, KELVIN, seed);
    md.prepare(field, &at);
    let e0 = md.kinetic_energy() + md.potential_energy().expect("prepared");
    let mut squares = 0.0;
    for _ in 0..steps {
        md.step(field, &mut at, dt);
        let e = md.kinetic_energy() + md.potential_energy().expect("stepped");
        squares += (e - e0) * (e - e0);
    }
    (squares / steps as f64).sqrt()
}

/// The three ratios of the RMS energy error per halving, 2 → 1, 1 → 0.5 and 0.5 → 0.25 fs, each
/// over 0.1 ps from the same start and velocities, and the RMS error at 2 fs.
fn ratios(seed: u64) -> ([f64; 3], f64) {
    let (field, at, md, _) = cluster(seed);
    let e = [(2.0, 50), (1.0, 100), (0.5, 200), (0.25, 400)]
        .map(|(h, n)| nve_rms(&field, &at, &md, h * FS, n, seed));
    ([e[0] / e[1], e[1] / e[2], e[2] / e[3]], e[0])
}

/// The bands, **measured** over sixteen starts by [`the_h_squared_band_measured`] (release): the
/// ratio from 2 fs 4.620–5.194, from 1 fs 4.131–4.253, from 0.5 fs 4.032–4.062, and the excess
/// over four shrinking by 0.234–0.262 per halving, `(r₃ − 4)/(r₂ − 4)`, which an `h²` integrator
/// with a positive `h⁴` term puts at a quarter and a defect of lower order with a small coefficient
/// pushes up. **2 fs is further from `h²` here than for benzene in water** (W3: 4.229–4.508): the
/// excess there is 0.6–1.2 against W3's 0.2–0.5, the stiffest motion now a protein's bends with a
/// hydrogen in them, which nothing holds. The ratio bands are about three times their measured
/// spreads. **The shrink's band is earned here, 0.20–0.30**, the measured 0.234–0.262 with 3.5
/// times its half-spread either side; it was W3's 0.15–0.35 and `Δq/h` × 0.9995 passed it at
/// 0.3326. **The shrink is the only guard against an error of order `h` with a small coefficient**:
/// `Δq/h` × 0.999 moved none of the three ratios out of its band.
const H_SQUARED_FROM_TWO: std::ops::Range<f64> = 4.2..5.8;
const H_SQUARED_FROM_ONE: std::ops::Range<f64> = 3.95..4.45;
const H_SQUARED_FROM_HALF: std::ops::Range<f64> = 3.99..4.10;
const H_SQUARED_SHRINKS: std::ops::Range<f64> = 0.20..0.30;

/// **Constrained velocity Verlet on a piece of the protein**: the 3 Å pocket's four residues, every
/// bond to a hydrogen held — the methyls' three sharing their carbon — benzene's six, and eight
/// rigid waters, NVE over 0.1 ps from 298 K velocities. The RMS energy error falls by four per
/// halving from 2 fs, and its excess over four shrinks by four per halving, each held to its band.
#[test]
fn constrained_nve_on_the_pocket_falls_as_h_squared() {
    let (r, e2) = ratios(0xA4);
    let shrink = (r[2] - 4.0) / (r[1] - 4.0);
    println!(
        "RMS ΔE over 0.1 ps at 2 fs {:.3e} kcal/mol; ratios {:.4}, {:.4}, {:.4}; the excess \
         shrinks by {shrink:.4}",
        kcal(e2),
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

/// The measurement the bands come from: sixteen starts, each a lattice turned by its own seed and
/// velocities from it. Release: `cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_h_squared_band_measured -- --ignored --nocapture`.
#[test]
#[ignore = "sixteen starts of four NVE runs: run with --release -- --ignored --nocapture"]
fn the_h_squared_band_measured() {
    let mut lo = [f64::INFINITY; 4];
    let mut hi = [f64::NEG_INFINITY; 4];
    for start in 0..16u64 {
        let (r, e2) = ratios(0xA4 + start);
        let shrink = (r[2] - 4.0) / (r[1] - 4.0);
        println!(
            "start {start}: ratios {:.4} {:.4} {:.4}, (r₃ − 4)/(r₂ − 4) {shrink:.4}, RMS at 2 fs \
             {:.3e} kcal/mol",
            r[0],
            r[1],
            r[2],
            kcal(e2)
        );
        for (k, x) in [r[0], r[1], r[2], shrink].into_iter().enumerate() {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    println!(
        "from 2 fs {:.4}–{:.4}, from 1 fs {:.4}–{:.4}, from 0.5 fs {:.4}–{:.4}, excess shrink \
         {:.4}–{:.4}",
        lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], lo[3], hi[3]
    );
}

// ---------------------------------------------------------------------------------------------
// The restraint's release

/// `∫ f` over `[a, b]` by Simpson's rule on `m` (even) intervals.
fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64, m: usize) -> f64 {
    let h = (b - a) / m as f64;
    let mut s = f(a) + f(b);
    for i in 1..m {
        let w = if i % 2 == 1 { 4.0 } else { 2.0 };
        s += w * f(a + i as f64 * h);
    }
    s * h / 3.0
}

/// **The Boresch release is its closed form**, at the anchors 3c's rule chooses on the 3 Å
/// pocket's box: `Boresch::release_free_energy_extended` against `−k_BT ln(8π² V°/Z_r)` with
/// `Z_r` Clark et al.'s eq 6 by quadrature — `r² e^{−βu_r}` over `[0, r₀ + 12σ_r]`, `sin θ e^{−βu_θ}`
/// over `[0, π]`, each dihedral over `(−π, π]`, Simpson's rule on 4000 intervals — to `10⁻⁹`
/// kcal/mol, where the Gaussians' tails beyond the ranges are below `e^{−20}` of each factor here;
/// and eq 7 (`release_free_energy`) apart from it by exactly
/// `k_BT [ln(1 + k_BT/(K_r r₀²)) − k_BT/2K_θA − k_BT/2K_θB]`, the closed form of the two
/// approximations, written out here.
#[test]
fn the_boresch_release_is_its_closed_form() {
    let s = small3(0xB0E);
    let restraint = anchors(s.binding(), &s.positions()[..s.complex_atoms()]);
    let kt = BOLTZMANN.to_si() * KELVIN;
    let [r0, ta, tb, ..] = restraint.reference;
    let k = restraint.force_constants;
    let sigma = |kk: f64| (kt / kk).sqrt();
    // Each factor of Z_r.
    let zr = simpson(
        |r| r * r * (-k[0] * (r - r0) * (r - r0) / (2.0 * kt)).exp(),
        0.0,
        r0 + 12.0 * sigma(k[0]),
        4000,
    );
    let angle = |t0: f64, kk: f64| {
        simpson(
            |t| t.sin() * (-kk * (t - t0) * (t - t0) / (2.0 * kt)).exp(),
            0.0,
            std::f64::consts::PI,
            4000,
        )
    };
    let dihedral = |kk: f64| {
        simpson(
            |p| (-kk * p * p / (2.0 * kt)).exp(),
            -std::f64::consts::PI,
            std::f64::consts::PI,
            4000,
        )
    };
    let z =
        zr * angle(ta, k[1]) * angle(tb, k[2]) * dihedral(k[3]) * dihedral(k[4]) * dihedral(k[5]);
    let pi2 = std::f64::consts::PI * std::f64::consts::PI;
    let quadrature = -kt * (8.0 * pi2 * pantometry_forcefield::boresch::STANDARD_VOLUME / z).ln();
    let extended = restraint.release_free_energy_extended(KELVIN);
    let eq7 = restraint.release_free_energy(KELVIN);
    let gap = kt * ((1.0 + kt / (k[0] * r0 * r0)).ln() - kt / (2.0 * k[1]) - kt / (2.0 * k[2]));
    println!(
        "anchors {:?} {:?}: r₀ {:.3} Å, θ_A {:.1}°, θ_B {:.1}°; release by quadrature {:+.10} \
         kcal/mol, extended eq 7 {:+.10}, eq 7 {:+.10}",
        restraint.receptor,
        restraint.ligand,
        r0 / ANGSTROM,
        ta.to_degrees(),
        tb.to_degrees(),
        kcal(quadrature),
        kcal(extended),
        kcal(eq7)
    );
    // Every angle at least 6σ from 0 and π, so that the tails are what the bound says.
    for t in [ta, tb] {
        assert!(t > 6.0 * sigma(k[1]) && std::f64::consts::PI - t > 6.0 * sigma(k[1]));
    }
    assert!(kcal(quadrature - extended).abs() < 1e-9);
    assert!(kcal(extended - eq7 - gap).abs() < 1e-12);
    assert!(
        kcal(gap).abs() > 1e-3,
        "the two forms differ by more than the check's tolerance"
    );
}

// ---------------------------------------------------------------------------------------------
// The whole complex: the measurements

/// Where the measurements write: `PANTOMETRY_W4_DIR`, or `pantometry-w4` in the system's temporary
/// directory. Never the repository.
fn results_dir() -> std::path::PathBuf {
    let dir = std::env::var_os("PANTOMETRY_W4_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("pantometry-w4"));
    std::fs::create_dir_all(&dir).expect("the results directory");
    dir
}

/// A lock that makes a second run in the same directory refuse, rather than write into the files
/// the first is writing. Removed when dropped; a run that was killed leaves it, and says so.
struct Lock(std::path::PathBuf);

impl Lock {
    fn take(dir: &std::path::Path, name: &str) -> Lock {
        let path = dir.join(format!("{name}.lock"));
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap_or_else(|e| {
                panic!(
                    "{} exists ({e}): another run is writing here, or one was killed — remove it \
                     by hand once no run is",
                    path.display()
                )
            });
        Lock(path)
    }
}

impl Drop for Lock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Appends a line to `name.log` and prints it.
fn log(dir: &std::path::Path, name: &str, line: &str) {
    use std::io::Write;
    println!("{line}");
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join(format!("{name}.log")))
        .expect("the log");
    writeln!(f, "{line}").expect("the log");
}

/// A float as text that parses back to the same bits: `{:e}`, Rust's shortest form that
/// round-trips. Every number the leg writes goes through this, and
/// `positions_round_trip_to_the_bit` holds it to the bit, subnormals, `−0.0` and the extremes
/// included.
fn exact(x: f64) -> String {
    format!("{x:e}")
}

/// FNV-1a, 64-bit, over each float's bits, little-endian: a digest that tells two runs' inputs
/// apart, written here so that nothing new is a dependency.
fn digest(xs: impl IntoIterator<Item = f64>) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for x in xs {
        for b in x.to_bits().to_le_bytes() {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
    }
    h
}

/// [`digest`] of positions, atom by atom, axis by axis.
fn positions_digest(at: &[[f64; 3]]) -> u64 {
    digest(at.iter().flatten().copied())
}

/// Writes `header` (one line), then the positions one atom a line, each number [`exact`], then a
/// trailer `# end N digest`, **to `path` with the extension `partial` and renamed to `path` once
/// complete**: a run killed while writing leaves no file at `path`, and a file cut short anywhere
/// is refused on reading ([`read_positions`]).
fn write_positions(path: &std::path::Path, header: &str, at: &[[f64; 3]]) {
    assert!(!header.contains('\n'), "a header is one line");
    let mut text = format!("{header}\n");
    for p in at {
        text.push_str(&format!(
            "{} {} {}\n",
            exact(p[0]),
            exact(p[1]),
            exact(p[2])
        ));
    }
    text.push_str(&format!(
        "# end {} {:016x}\n",
        at.len(),
        positions_digest(at)
    ));
    let partial = path.with_extension("partial");
    std::fs::write(&partial, text).expect("positions");
    std::fs::rename(&partial, path).expect("rename");
}

/// A positions file's header line and its positions.
type Header = (String, Vec<[f64; 3]>);

/// What [`write_positions`] wrote at `path`: `None` when there is no file; its header and
/// positions when it is whole; and why not when it is not — no trailer, a line that is not three
/// numbers, or a count or digest the positions read back do not match.
fn read_positions(path: &std::path::Path) -> Option<Result<Header, String>> {
    let text = std::fs::read_to_string(path).ok()?;
    Some(parse_positions(&text))
}

/// [`read_positions`]'s reading of a file's text.
fn parse_positions(text: &str) -> Result<Header, String> {
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| "an empty file".to_string())?
        .to_string();
    let mut at = Vec::new();
    for line in lines {
        if let Some(end) = line.strip_prefix("# end ") {
            let mut words = end.split_whitespace();
            let n: Option<usize> = words.next().and_then(|x| x.parse().ok());
            let d = words.next().and_then(|x| u64::from_str_radix(x, 16).ok());
            return match (n, d) {
                (Some(n), _) if n != at.len() => {
                    Err(format!("{} positions where the trailer says {n}", at.len()))
                }
                (Some(_), Some(d)) if d == positions_digest(&at) => Ok((header, at)),
                (Some(_), Some(_)) => Err("the positions' digest is not the trailer's".into()),
                _ => Err(format!("a trailer cut short: {line:?}")),
            };
        }
        let v: Result<Vec<f64>, _> = line.split_whitespace().map(str::parse).collect();
        match v {
            Ok(v) if v.len() == 3 => at.push([v[0], v[1], v[2]]),
            _ => return Err(format!("not a position: {line:?}")),
        }
    }
    Err("no trailer: the file is not whole".into())
}

/// **The leg's hydrogens**: `binding`'s, each system's own, relaxed in vacuum with every heavy
/// atom frozen — 3b's and 3c's `relaxing_hydrogens(20 000, HYDROGEN_TOLERANCE)`, step A-1 —
/// before it goes in water. Phase 1 put the hydrogens in water as 2c-1 placed them, and the
/// crystal's contacts, a backbone amide H 1.50 Å from an Asp Oδ1 among them, were what a released
/// protein let go of at its first step.
fn leg_relaxed(binding: &Binding) -> Binding {
    binding.relaxing_hydrogens(20_000, Binding::HYDROGEN_TOLERANCE)
}

/// The leg's binding: 181L's residues within `cutoff` of benzene, hydrogens relaxed
/// ([`leg_relaxed`]).
fn leg_binding(cutoff: f64) -> Binding {
    leg_relaxed(&Binding::new(the_system(), cutoff).expect("181L"))
}

/// The leg's box: `binding` in water as `solvation` says, **its real-space pairs in a Verlet list**
/// of [`PeriodicForceField::NEIGHBOUR_SKIN`], which the decoupling keeps a list of its own of.
fn leg_box(binding: &Binding, solvation: &Solvation) -> SolvatedComplex {
    SolvatedComplex::new(binding, solvation)
        .expect("the box")
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN)
}

/// How much of the protein the leg moves: all of it, decided after phase 1's measurement.
const LEG_FLEXIBILITY: Flexibility = Flexibility::Mobile;

/// The seed of the leg's preparation ([`SolvatedComplex::prepare`]).
const PREPARATION_SEED: u64 = 0x4_4A7F;

/// The whole of 181L as a binding: every residue (cut at 100 Å, past the farthest heavy atom from
/// benzene), QEq on the protein at +9 and on benzene at 0 — benzene's charges 3c-3's and W3's —
/// and the hydrogens as 2c-1 placed them, **not yet relaxed**.
fn whole_binding() -> Binding {
    let b = Binding::new(the_system(), 100.0 * ANGSTROM).expect("the whole of 181L");
    assert_eq!(b.positions().len(), the_system().component().atoms().len());
    assert!(
        b.system_atoms().iter().enumerate().all(|(k, &i)| k == i),
        "the whole binding is the system, in its order"
    );
    // Benzene's charges are the 3 Å binding's, which W3 used, to the bit.
    assert_eq!(
        b.charges()[b.ligand_range()],
        pocket3().charges()[pocket3().ligand_range()]
    );
    b
}

/// The whole of 181L in W4's box as the leg runs it: [`whole_binding`], its hydrogens relaxed
/// ([`leg_relaxed`]), in [`leg_box`].
fn whole_complex() -> SolvatedComplex {
    leg_box(&leg_relaxed(&whole_binding()), &Solvation::w4())
}

/// The header [`prepared`] writes and reads back: every stage literally (flexibility, time step,
/// steps, friction), the seed and every stage's seed as [`SolvatedComplex::stage_seed`] derives
/// it, the temperature, the atom count, and a digest of the box's positions **as built**, before
/// any stage — which moves with the hydrogens' relaxation, the solvation and the lattice.
fn prepared_header(s: &SolvatedComplex, stages: &[Stage], seed: u64, temperature: f64) -> String {
    let st: Vec<String> = stages
        .iter()
        .map(|st| {
            format!(
                "{:?}/{}/{}/{}",
                st.flexibility,
                exact(st.time_step),
                st.steps,
                exact(st.friction)
            )
        })
        .collect();
    let seeds: Vec<String> = (0..stages.len())
        .map(|k| format!("{:016x}", SolvatedComplex::stage_seed(seed, k)))
        .collect();
    let seeding = format!("seed {seed:#x} stage-seeds {}", seeds.join(","));
    format!(
        "# prepared stages {} {seeding} temperature {} atoms {} built {:016x}",
        st.join(","),
        exact(temperature),
        s.positions().len(),
        positions_digest(s.positions())
    )
}

/// `stages` of [`SolvatedComplex::prepare`] on `potential` from the box as built, at `seed` and
/// `temperature`, each stage logged as it ends, and the end written to `file` in `dir`
/// ([`write_positions`], atomically) under [`prepared_header`]. **Read back instead when a run has
/// written it with the same header**; a file prepared differently, or not whole, is refused by
/// name, never reused.
#[allow(clippy::too_many_arguments)]
fn prepared(
    s: &SolvatedComplex,
    potential: &dyn pantometry_forcefield::Potential,
    stages: &[Stage],
    seed: u64,
    temperature: f64,
    dir: &std::path::Path,
    file: &str,
    name: &str,
) -> Vec<[f64; 3]> {
    let path = dir.join(file);
    let want = prepared_header(s, stages, seed, temperature);
    if let Some(read) = read_positions(&path) {
        let (head, at) = read.unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it to prepare again",
                path.display()
            )
        });
        assert_eq!(
            head,
            want,
            "{} was prepared differently: remove it, and every window_NN.txt run from it, to \
             prepare again",
            path.display()
        );
        assert_eq!(at.len(), s.positions().len());
        return at;
    }
    let t = std::time::Instant::now();
    let mut at = s.positions().to_vec();
    s.prepare(
        potential,
        stages,
        temperature,
        seed,
        &mut at,
        |k, stage, md, _| {
            log(
                dir,
                name,
                &format!(
                    "stage {k}: {} steps of {} fs, {:?}, {} ps⁻¹: {:.0} s in all, T {:.1} K, U \
                     {:.1} kcal/mol",
                    stage.steps,
                    stage.time_step / FS,
                    stage.flexibility,
                    stage.friction / 1e12,
                    t.elapsed().as_secs_f64(),
                    md.temperature(),
                    kcal(md.potential_energy().unwrap())
                ),
            )
        },
    );
    write_positions(&path, &want, &at);
    at
}

/// The candidate states, in path order: the restraint on, the charges off, the van der Waals off
/// — 3c-3's grid, the van der Waals's tail W3's. Every window records its energy at every
/// candidate (the couplings are the cross terms alone), so a window inserted between two that ran
/// needs neither run again.
fn grid() -> Vec<Lambda> {
    let mut g = Vec::new();
    for r in [0.0, 0.05, 0.1, 0.175, 0.25, 0.375, 0.5, 0.75, 1.0] {
        g.push(Lambda::new(r, 1.0, 1.0));
    }
    for e in [0.75, 0.5, 0.25, 0.0] {
        g.push(Lambda::new(1.0, e, 1.0));
    }
    for k in 1..=20 {
        // 0.95, 0.9, …, 0: twenty steps of 0.05, written as integers over 20.
        g.push(Lambda::new(1.0, 0.0, f64::from(20 - k) / 20.0));
    }
    g
}

/// The complex leg's schedule, as indices into [`grid`]: **the restraint on in one interval**,
/// λ_r 0 → 1, which the insertion splits where it overlaps too little; the charges off in two
/// (W3's); the van der Waals off in ten (W3's). Fourteen windows before any is inserted.
const SCHEDULE: [usize; 14] = [0, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32];

/// W3's protocol: 2 fs, 1 ps⁻¹ at 298.15 K, 20 ps discarded, 2000 samples every 100 fs.
const PROTOCOL: Protocol = Protocol {
    time_step: 2e-15,
    temperature: KELVIN,
    friction: 1e12,
    equilibration: 10_000,
    stride: 50,
    samples: 2000,
    seed: 0x4_4A7E,
};

/// An interval whose overlap is below this gets a window at the candidate between its ends: 3c's
/// and W3's rule and threshold ([`pantometry_forcefield::free_energy::refine_schedule`]).
const OVERLAP_THRESHOLD: f64 = 0.1;

/// How many rounds of insertion the leg makes: W3's two.
const INSERTION_ROUNDS: usize = 2;

/// W3's solvent leg, benzene decoupled from TIP3P: +1.570 ± 0.140 kcal/mol (CHANGELOG, W3).
const SOLVENT_LEG: (f64, f64) = (1.570, 0.140);

/// Constrained NVE from `at` with `velocities` under `template` (no bath): `n` steps of `dt`. The
/// RMS departure of the total energy from its start, the drift by a least-squares line per
/// second, ⟨T⟩ and the worst held bond with a mobile atom, relative to its length.
#[allow(clippy::too_many_arguments)]
fn nve(
    s: &SolvatedComplex,
    potential: &dyn pantometry_forcefield::Potential,
    template: &MolecularDynamics,
    mobile: &[bool],
    at: &[[f64; 3]],
    velocities: &[[f64; 3]],
    dt: f64,
    n: usize,
) -> (f64, f64, f64, f64) {
    let mut at = at.to_vec();
    let mut md = template.clone().with_velocities(velocities.to_vec());
    md.prepare(potential, &at);
    let e0 = md.kinetic_energy() + md.potential_energy().unwrap();
    let (mut squares, mut series, mut t_sum) = (0.0, Vec::with_capacity(n), 0.0);
    for _ in 0..n {
        md.step(potential, &mut at, dt);
        let e = md.kinetic_energy() + md.potential_energy().unwrap() - e0;
        squares += e * e;
        series.push(e);
        t_sum += md.temperature();
    }
    // The drift: a least-squares line through the energy against time.
    let m = n as f64;
    let tm = (m + 1.0) / 2.0;
    let em = series.iter().sum::<f64>() / m;
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for (k, e) in series.iter().enumerate() {
        let x = k as f64 + 1.0 - tm;
        sxy += x * (e - em);
        sxx += x * x;
    }
    let worst_bond = s
        .shake()
        .bonds()
        .iter()
        .zip(s.shake().lengths())
        .filter(|(&[i, j], _)| mobile[i] || mobile[j])
        .map(|(&[i, j], &d0)| (len(sub(at[i], at[j])) - d0).abs() / d0)
        .fold(0.0f64, f64::max);
    ((squares / m).sqrt(), sxy / sxx / dt, t_sum / m, worst_bond)
}

/// One flexibility's cost, as a line for the log: the mobile atoms and the degrees of freedom; a
/// protein that moves released from the frozen start ([`SolvatedComplex::release`]); ms per step
/// at 2 fs in a 1 ps⁻¹ bath, the fastest and slowest of three runs of ten steps; then constrained
/// NVE over 0.3 ps at 2 fs and at 1 fs from the same positions and velocities ([`nve`]).
fn flexibility_cost(
    s: &SolvatedComplex,
    d: &PeriodicDecoupling,
    start: &[[f64; 3]],
    flexibility: Flexibility,
    couplings: f64,
) -> String {
    let mobile = s.mobile(flexibility);
    let free_complex = mobile[..s.complex_atoms()].iter().filter(|&&m| m).count();
    let md0 = s.dynamics(flexibility);
    let potential = AtLambda {
        hamiltonian: d,
        lambda: Lambda::COUPLED,
    };
    let mut at = start.to_vec();
    let md = if flexibility != Flexibility::Frozen {
        s.prepare(
            &potential,
            &SolvatedComplex::release(flexibility),
            KELVIN,
            0x5EED,
            &mut at,
            |_, _, _, _| {},
        )
    } else {
        md0.clone().thermalised(&at, KELVIN, 0x5EED)
    };
    let mut md = md.with_bath(Bath::Langevin {
        temperature: KELVIN,
        friction: 1e12,
        seed: 0x5EEF,
    });
    let mut chunks = Vec::new();
    for _ in 0..3 {
        let t = std::time::Instant::now();
        md.run(&potential, &mut at, 2.0 * FS, 10);
        chunks.push(t.elapsed().as_secs_f64() / 10.0);
    }
    let fastest = chunks.iter().copied().fold(f64::INFINITY, f64::min);
    let slowest = chunks.iter().copied().fold(0.0, f64::max);
    let dof = md.degrees_of_freedom();
    let v = md.velocities().to_vec();
    let two = nve(s, &potential, &md0, &mobile, &at, &v, 2.0 * FS, 150);
    let one = nve(s, &potential, &md0, &mobile, &at, &v, 1.0 * FS, 300);
    let windows = SCHEDULE.len() as f64;
    let leg = |step: f64| {
        windows * PROTOCOL.steps_per_window() as f64 * step
            + windows * PROTOCOL.samples as f64 * couplings
    };
    format!(
        "{free_complex} of {} complex atoms mobile, {dof} degrees of freedom; {:.1}–{:.1} ms a \
         step at 2 fs; NVE 0.3 ps: at 2 fs RMS {:.4} kcal/mol, drift {:+.4} kcal/mol/ps, ⟨T⟩ {:.1} \
         K, worst held bond {:.1e}; at 1 fs RMS {:.4}, drift {:+.4}, ⟨T⟩ {:.1}; RMS ratio {:.2}; \
         the leg, {} windows of {} steps: {:.0}–{:.0} h",
        s.complex_atoms(),
        fastest * 1e3,
        slowest * 1e3,
        kcal(two.0),
        kcal(two.1) * 1e-12,
        two.2,
        two.3,
        kcal(one.0),
        kcal(one.1) * 1e-12,
        one.2,
        two.0 / one.0,
        SCHEDULE.len(),
        PROTOCOL.steps_per_window(),
        leg(fastest) / 3600.0,
        leg(slowest) / 3600.0
    )
}

/// **What each way of holding the protein costs** (phase 1's measurement, release, one core, about
/// twenty-five minutes), now on the box the leg runs ([`whole_complex`]: hydrogens relaxed, the
/// neighbour list on): the whole of 181L in W4's box melted with the protein frozen
/// ([`SolvatedComplex::MELT`]); then for the protein frozen, residues beyond 8 and 12 Å of benzene
/// frozen, and nothing frozen, [`flexibility_cost`]; and the evaluation's parts and the couplings'
/// cost. A flexibility whose dynamics fails is logged as failed, with why, and the others still
/// run. `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_cost_of_each_flexibility_measured -- --ignored --nocapture`.
#[test]
#[ignore = "the whole complex in 9 160 waters, about twenty-five minutes with --release: run with --release -- --ignored --nocapture"]
fn the_cost_of_each_flexibility_measured() {
    let dir = results_dir();
    let name = "cost";
    let _lock = Lock::take(&dir, name);
    let t = std::time::Instant::now();
    let s = whole_complex();
    let mesh = s.field().ewald().mesh().unwrap().parameters();
    let l = s.cell().lengths();
    log(
        &dir,
        name,
        &format!(
            "181L in water, built in {:.0} s: {} atoms — {} complex, {} waters, {} Cl- — lattice \
             {:?} ({} waters, {} overlapping removed); box {:.3} × {:.3} × {:.3} Å; α {:.4} Å⁻¹, \
             mesh order {} on {:?}; {} held bonds",
            t.elapsed().as_secs_f64(),
            s.positions().len(),
            s.complex_atoms(),
            s.waters(),
            s.ions(),
            s.lattice(),
            s.lattice().iter().product::<usize>(),
            s.overlapping(),
            l[0] / ANGSTROM,
            l[1] / ANGSTROM,
            l[2] / ANGSTROM,
            s.field().ewald().parameters().alpha * ANGSTROM,
            mesh.order,
            mesh.grid,
            s.shake().len()
        ),
    );
    let q: f64 = s.field().charges().iter().sum();
    log(&dir, name, &format!("net charge {q:+.3e} e"));
    let d = s.decoupling();
    let coupled = AtLambda {
        hamiltonian: &d,
        lambda: Lambda::COUPLED,
    };
    let start = prepared(
        &s,
        &coupled,
        &SolvatedComplex::MELT,
        PREPARATION_SEED,
        KELVIN,
        &dir,
        "start.txt",
        name,
    );

    // The evaluation's parts.
    // The fastest of `k`: the machine's load only ever adds.
    let time = |f: &dyn Fn(), k: usize| {
        let mut best = f64::INFINITY;
        for _ in 0..k {
            let t = std::time::Instant::now();
            f();
            best = best.min(t.elapsed().as_secs_f64());
        }
        best
    };
    let full = time(&|| drop(s.field().evaluate(&start)), 5);
    let no_charges = s
        .field()
        .clone()
        .with_charges(vec![0.0; s.positions().len()]);
    let vdw_only = time(&|| drop(no_charges.evaluate(&start)), 5);
    let bonded = s.field().bonded().unwrap();
    let complex_at = &start[..s.complex_atoms()];
    let bonded_only = time(&|| drop(bonded.evaluate(complex_at)), 5);
    let states = grid();
    let couplings = time(&|| drop(d.couplings(&start, &states)), 5);
    log(
        &dir,
        name,
        &format!(
            "an evaluation {:.1} ms: van der Waals and the bonded terms without charges {:.1}, the \
             protein's bonded terms {:.1}; the couplings at all {} candidates {:.1} ms",
            full * 1e3,
            vdw_only * 1e3,
            bonded_only * 1e3,
            states.len(),
            couplings * 1e3,
        ),
    );

    let mut schemes = vec![("protein frozen", Flexibility::Frozen)];
    for r in [8.0, 12.0] {
        schemes.push(("zone", Flexibility::Zone(r * ANGSTROM)));
    }
    schemes.push(("nothing frozen", Flexibility::Mobile));
    for (label, flexibility) in schemes {
        let outcome = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            flexibility_cost(&s, &d, &start, flexibility, couplings)
        }));
        match outcome {
            Ok(line) => log(&dir, name, &format!("{label} {flexibility:?}: {line}")),
            Err(e) => {
                let why = e
                    .downcast_ref::<String>()
                    .cloned()
                    .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default();
                log(
                    &dir,
                    name,
                    &format!("{label} {flexibility:?}: FAILED: {why}"),
                );
            }
        }
    }
    log(
        &dir,
        name,
        &format!("the measurement took {:.0} s", t.elapsed().as_secs_f64()),
    );
}

// ---------------------------------------------------------------------------------------------
// Phase 2: the leg's pieces

/// The `k` closest contacts of a hydrogen with an atom of the complex at `at` that is neither
/// bonded to it nor bonded to its neighbour — the 1-2 and 1-3 pairs, which van der Waals leaves
/// out — by brute force: `(distance, hydrogen, other)`, nearest first, each pair once.
fn closest_hydrogen_contacts(b: &Binding, at: &[[f64; 3]], k: usize) -> Vec<(f64, usize, usize)> {
    let n = at.len();
    let el = b.elements();
    let mut neighbours = vec![Vec::new(); n];
    for bond in b.component().bonds() {
        let [i, j] = bond.atoms;
        neighbours[i].push(j);
        neighbours[j].push(i);
    }
    let mut out = Vec::new();
    for h in (0..n).filter(|&i| el[i] == Element::H) {
        let mut excluded = vec![h];
        for &j in &neighbours[h] {
            excluded.push(j);
            excluded.extend(neighbours[j].iter().copied());
        }
        for x in 0..n {
            if excluded.contains(&x) || (el[x] == Element::H && x < h) {
                continue;
            }
            out.push((len(sub(at[h], at[x])), h, x));
        }
    }
    out.sort_by(|a, b| a.0.total_cmp(&b.0));
    out.truncate(k);
    out
}

/// A contact for a log: `A:THR54 H … A:ASP47 OD1 1.505 Å`.
fn contact_label(b: &Binding, c: (f64, usize, usize)) -> String {
    let s = the_system();
    let name = |k: usize| {
        let i = b.system_atoms()[k];
        format!("{} {}", s.residue_of(i).label(), s.atom_name(i))
    };
    format!("{} … {} {:.3} Å", name(c.1), name(c.2), c.0 / ANGSTROM)
}

/// **The leg's hydrogens are relaxed in vacuum before it is solvated**, on the 5 Å pocket through
/// the leg's own [`leg_binding`]: each of the three relaxations converged; **the forces on the
/// freed hydrogens, recomputed from the complex's whole force field, are each within the
/// tolerance** — a stationary point in what was freed, which a binding not relaxed is not (the
/// largest is asserted to be ten times the tolerance there); every heavy atom and every hydrogen
/// left held keeps its bits; the closest contact a freed hydrogen makes is longer than at the
/// placement, and the energy lower; and the box built from it starts from those positions.
#[test]
fn the_legs_hydrogens_are_relaxed_in_vacuum_first() {
    let b = leg_binding(5.0 * ANGSTROM);
    let placed = pocket5();
    let r = b
        .hydrogen_relaxation()
        .expect("the leg's binding is relaxed");
    for (system, p) in [
        ("complex", r.complex),
        ("pocket", r.pocket),
        ("ligand", r.ligand),
    ] {
        println!("{system}: {:?} in {} steps", p.status, p.steps);
        assert_eq!(
            p.status,
            pantometry_forcefield::Status::Converged,
            "{system}"
        );
    }
    assert_eq!(r.tolerance, Binding::HYDROGEN_TOLERANCE);
    assert_eq!(r.free, b.free_hydrogens());
    // Which hydrogens are free, by brute force: every ligand hydrogen, and each pocket hydrogen
    // whose bonded heavy atom is within the cutoff of a ligand atom at the crystal positions.
    let crystal = b.crystal_positions();
    let ligand = b.ligand_range();
    let within = |k: usize| {
        ligand
            .clone()
            .any(|l| len(sub(crystal[k], crystal[l])) < b.cutoff())
    };
    let parent = |h: usize| {
        b.component()
            .bonds()
            .iter()
            .find_map(|bond| match bond.atoms {
                [i, j] if i == h => Some(j),
                [i, j] if j == h => Some(i),
                _ => None,
            })
            .expect("a hydrogen is bonded")
    };
    let brute: Vec<bool> = (0..b.positions().len())
        .map(|k| b.elements()[k] == Element::H && (ligand.contains(&k) || within(parent(k))))
        .collect();
    assert_eq!(r.free, brute, "the freed hydrogens are the rule's");
    let freed = brute.iter().filter(|&&f| f).count();
    let hydrogens = b.elements().iter().filter(|&&e| e == Element::H).count();
    println!("{freed} of {hydrogens} hydrogens freed");
    assert!(
        freed > 0 && freed < hydrogens,
        "the rule frees some and holds some"
    );
    // **The leg's 100 Å binding frees every hydrogen**, by the same rule on the system's own
    // positions, without its QEq: every protein hydrogen's heavy atom is within 100 Å of benzene.
    let system = the_system();
    let atoms = system.component().atoms();
    let ligand_atoms = system.atoms_in(pantometry_forcefield::Part::Ligand);
    let mut farthest = 0.0f64;
    for bond in system.component().bonds() {
        for (h, x) in [
            (bond.atoms[0], bond.atoms[1]),
            (bond.atoms[1], bond.atoms[0]),
        ] {
            if atoms[h].element == Element::H && atoms[x].element != Element::H {
                let d = ligand_atoms
                    .iter()
                    .map(|&l| len(sub(atoms[x].at, atoms[l].at)))
                    .fold(f64::INFINITY, f64::min);
                farthest = farthest.max(d);
            }
        }
    }
    println!(
        "a hydrogen's heavy atom is at most {:.2} Å from benzene",
        farthest / ANGSTROM
    );
    assert!(farthest < 100.0 * ANGSTROM);
    let worst = |at: &[[f64; 3]]| {
        let f = b.force_field().evaluate(at).forces;
        (0..at.len())
            .filter(|&k| r.free[k])
            .map(|k| len(f[k]))
            .fold(0.0f64, f64::max)
    };
    let (after, before) = (worst(b.positions()), worst(placed.positions()));
    println!(
        "the largest force on a freed hydrogen: {:.3e} kcal/mol/Å relaxed, {:.3e} as placed",
        after / pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM,
        before / pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM
    );
    assert!(after <= Binding::HYDROGEN_TOLERANCE);
    assert!(before > 10.0 * Binding::HYDROGEN_TOLERANCE);
    let mut moved = 0;
    for k in 0..b.positions().len() {
        if r.free[k] {
            moved += usize::from(b.positions()[k] != placed.positions()[k]);
        } else {
            assert_eq!(
                b.positions()[k],
                placed.positions()[k],
                "atom {k} kept its bits"
            );
        }
    }
    assert!(moved > 0);
    // The contacts a freed hydrogen makes: the held ones cannot move.
    let near = |at: &[[f64; 3]]| -> Vec<(f64, usize, usize)> {
        closest_hydrogen_contacts(&b, at, 100)
            .into_iter()
            .filter(|c| r.free[c.1] || r.free[c.2])
            .take(3)
            .collect()
    };
    let (c0, c1) = (near(placed.positions()), near(b.positions()));
    for (label, c) in [("placed", &c0), ("relaxed", &c1)] {
        let text: Vec<String> = c.iter().map(|&x| contact_label(&b, x)).collect();
        println!("{label}: {}", text.join("; "));
    }
    assert!(c1[0].0 > c0[0].0, "the closest contact is relieved");
    let e = |at: &[[f64; 3]]| b.force_field().energy(at).total;
    assert!(e(b.positions()) < e(placed.positions()));
    let s = leg_box(&b, &small_solvation(3.5, 5.5, 0x4E1));
    let given = b.positions();
    let shift = sub(s.positions()[0], given[0]);
    assert_eq!(
        b.elements()[0],
        Element::N,
        "the first atom is heavy, so not moved onto a bond"
    );
    let edge = s.cell().lengths()[0];
    // Each freed hydrogen in the box is where the relaxation put it, less only SHAKE's move onto
    // its bond's length — at most 0.05 Å, phase 1 measured 0.046 on this pocket — and, for every
    // one the relaxation moved more than 0.1 Å, nearer its relaxed place than its placed one.
    let (mut farthest, mut far_moved) = (0.0f64, 0);
    for k in 0..given.len() {
        let in_box = sub(s.positions()[k], shift);
        if b.elements()[k] != Element::H {
            assert!(len(sub(in_box, given[k])) <= 4.0 * EPS * edge);
        } else if r.free[k] {
            let to_relaxed = len(sub(in_box, given[k]));
            farthest = farthest.max(to_relaxed);
            assert!(to_relaxed <= 0.05 * ANGSTROM, "hydrogen {k} {to_relaxed:e}");
            if len(sub(given[k], placed.positions()[k])) > 0.1 * ANGSTROM {
                far_moved += 1;
                assert!(
                    to_relaxed < len(sub(in_box, placed.positions()[k])),
                    "hydrogen {k}"
                );
            }
        }
    }
    println!(
        "in the box, a freed hydrogen at most {:.4} Å from its relaxed place; {far_moved} moved \
         more than 0.1 Å by the relaxation",
        farthest / ANGSTROM
    );
    assert!(far_moved > 0);
}

/// **W4's solvation is the numbers it says**, each field pinned to a literal: a 12 Å margin,
/// 33.00 waters per nm³, the 1.4 Å clearance and the oxygen's 2.6 Å, `r_c` = 9 Å, δ = 10⁻⁶, the
/// seed, and chloride.
#[test]
fn w4s_solvation_is_its_numbers() {
    let w4 = Solvation::w4();
    assert_eq!(w4.margin, 12.0 * ANGSTROM);
    assert_eq!(w4.density, 33.0e27 * water::molecular_mass());
    assert_eq!(w4.clearance, 1.4 * ANGSTROM);
    assert_eq!(w4.oxygen_clearance, 2.6 * ANGSTROM);
    assert_eq!(w4.cutoff, 9.0 * ANGSTROM);
    assert_eq!(w4.accuracy, 1e-6);
    assert_eq!(w4.seed, 0x4_4A7E);
    assert_eq!(w4.counter_ion, pantometry_forcefield::Ion::chloride());
}

/// **The oxygen's rule at its two edges**, on one lattice water of a 25 Å box and two carbon atoms
/// placed about its oxygen by hand: the angle `a–O–b` at 89° and at 91°, each atom at 2.59 Å or
/// 2.61 Å from the oxygen, either side of W4's 2.6 Å. Only both atoms inside 2.6 Å **and** more
/// than 90° apart removes the water; every hydrogen of it is more than 1.4 Å from both, so the
/// clearance removes none. The atoms are put along the bisector of the angle, in a plane away
/// from the water's hydrogens.
#[test]
fn the_oxygens_rule_is_between_two_atoms_past_ninety_degrees_within_two_point_six() {
    let w4 = Solvation::w4();
    let lattice = WaterBox::lattice_box([8, 8, 8], w4.density, w4.seed);
    let keep = 3 * 64 + 3 * 8 + 3;
    let drop: Vec<bool> = (0..lattice.count()).map(|w| w != keep).collect();
    let one = lattice.without_waters(&drop);
    let water = [one.positions()[0], one.positions()[1], one.positions()[2]];
    let o = water[0];
    // The plane: perpendicular to the water's own plane, through the oxygen, along the H–O–H
    // bisector's opposite, so both carbons are on the side away from the hydrogens.
    let mid = [0, 1, 2].map(|k| 0.5 * (water[1][k] + water[2][k]) - o[k]);
    let away = [-mid[0] / len(mid), -mid[1] / len(mid), -mid[2] / len(mid)];
    let hh = sub(water[1], water[2]);
    let n = [
        mid[1] * hh[2] - mid[2] * hh[1],
        mid[2] * hh[0] - mid[0] * hh[2],
        mid[0] * hh[1] - mid[1] * hh[0],
    ];
    let side = [n[0] / len(n), n[1] / len(n), n[2] / len(n)];
    let elements = [Element::C, Element::C];
    let mut seen = Vec::new();
    for degrees in [89.0f64, 91.0] {
        let half = (0.5 * degrees).to_radians();
        let (s, c) = half.sin_cos();
        let u = [0, 1, 2].map(|k| c * away[k] + s * side[k]);
        let v = [0, 1, 2].map(|k| c * away[k] - s * side[k]);
        for (ra, rb) in [(2.59, 2.59), (2.59, 2.61), (2.61, 2.61)] {
            let a = [0, 1, 2].map(|k| o[k] + ra * ANGSTROM * u[k]);
            let b = [0, 1, 2].map(|k| o[k] + rb * ANGSTROM * v[k]);
            let angle = (dot(sub(a, o), sub(b, o)) / (len(sub(a, o)) * len(sub(b, o))))
                .acos()
                .to_degrees();
            assert!((angle - degrees).abs() < 1e-9, "{angle}");
            let clear = water
                .iter()
                .all(|&p| len(sub(p, a)) > 1.4 * ANGSTROM && len(sub(p, b)) > 1.4 * ANGSTROM);
            assert!(clear, "no water atom within the clearance");
            let removed = w4.removed(&one, &[a, b], &elements)[0];
            println!("{degrees}°, {ra} and {rb} Å: removed {removed}");
            seen.push(removed);
        }
    }
    assert_eq!(seen, [false, false, false, true, false, false]);
}

/// **The preparation draws on other streams than the windows**: no stage of the leg's
/// preparation has the seed of any window the leg could run, at any of the grid's candidates.
#[test]
fn the_preparation_draws_on_other_streams_than_the_windows() {
    let stages = SolvatedComplex::preparation(LEG_FLEXIBILITY).len();
    let prepared: std::collections::BTreeSet<u64> = (0..stages)
        .map(|k| SolvatedComplex::stage_seed(PREPARATION_SEED, k))
        .collect();
    let windows: std::collections::BTreeSet<u64> = (0..grid().len())
        .map(|g| window_seed(PROTOCOL.seed, g))
        .collect();
    assert_eq!(prepared.len(), stages);
    assert!(prepared.is_disjoint(&windows));
}

/// Phase 1's crevice water: lattice water 4663 of W4's box around the whole of 181L, which the
/// melt pushed into the crevice between Arg95 and Trp126 until a 2 fs step turned it past what
/// SETTLE solves.
const CREVICE_WATER: usize = 4663;

/// **The crevice water between Arg95 and Trp126 is removed, and by the rule that is for it**: W4's
/// lattice around the whole of 181L as phase 1 built it — the crystal's positions, hydrogens as
/// placed, the 20 × 21 × 24 lattice — and lattice water [`CREVICE_WATER`], checked by brute force
/// over 27 images to be that water: its oxygen within 2.6 Å of heavy atoms of Arg95 and of Trp126,
/// two of them on opposite sides of it. **The clearance alone keeps it** (every atom of it 1.4 Å
/// or more from every complex atom), and the oxygen's rule, between two heavy atoms, removes it.
/// No QEq is needed: the box's waters are a matter of positions.
#[test]
fn the_crevice_water_between_arg95_and_trp126_is_removed() {
    let system = the_system();
    let atoms = system.component().atoms();
    let at: Vec<[f64; 3]> = atoms.iter().map(|a| a.at).collect();
    let elements: Vec<Element> = atoms.iter().map(|a| a.element).collect();
    let w4 = Solvation::w4();
    let (k, shift) = w4.lattice_for(&at);
    assert_eq!(k, [20, 21, 24], "phase 1's lattice");
    let lattice = WaterBox::lattice_box(k, w4.density, w4.seed);
    let complex: Vec<[f64; 3]> = at
        .iter()
        .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
        .collect();
    let only: Vec<bool> = (0..lattice.count()).map(|w| w != CREVICE_WATER).collect();
    let one = lattice.without_waters(&only);
    let water = &lattice.positions()[3 * CREVICE_WATER..3 * CREVICE_WATER + 3];
    assert_eq!(one.positions(), water);
    let l = lattice.cell().lengths();
    let nearest = |p: [f64; 3], q: [f64; 3]| {
        let mut best = [f64::INFINITY; 3];
        for i in -1..=1 {
            for j in -1..=1 {
                for m in -1..=1 {
                    let image = [
                        q[0] + f64::from(i) * l[0],
                        q[1] + f64::from(j) * l[1],
                        q[2] + f64::from(m) * l[2],
                    ];
                    let v = sub(p, image);
                    if dot(v, v) < dot(best, best) {
                        best = v;
                    }
                }
            }
        }
        best
    };
    let closest = water
        .iter()
        .flat_map(|&p| complex.iter().map(move |&c| len(nearest(p, c))))
        .fold(f64::INFINITY, f64::min);
    let close: Vec<(usize, [f64; 3])> = (0..complex.len())
        .filter(|&i| elements[i] != Element::H)
        .map(|i| (i, nearest(water[0], complex[i])))
        .filter(|(_, v)| len(*v) < w4.oxygen_clearance)
        .collect();
    let residues: std::collections::BTreeSet<String> = close
        .iter()
        .map(|&(i, _)| system.residue_of(i).label())
        .collect();
    let mut most_opposite = 1.0f64;
    for (a, &(_, u)) in close.iter().enumerate() {
        for &(_, v) in &close[..a] {
            most_opposite = most_opposite.min(dot(u, v) / (len(u) * len(v)));
        }
    }
    println!(
        "water {CREVICE_WATER}: {:.3} Å from the complex at closest; {} heavy atoms within {:.3} \
         Å of its oxygen, of {residues:?}, the widest pair at cos {most_opposite:.3}",
        closest / ANGSTROM,
        close.len(),
        w4.oxygen_clearance / ANGSTROM
    );
    assert!(residues.contains("A:ARG95") && residues.contains("A:TRP126"));
    assert!(most_opposite < 0.0, "two of them on opposite sides");
    assert!(closest >= w4.clearance, "the clearance alone keeps it");
    let clearance_only = Solvation {
        oxygen_clearance: 0.0,
        ..w4
    };
    assert_eq!(clearance_only.removed(&one, &complex, &elements), [false]);
    assert_eq!(w4.removed(&one, &complex, &elements), [true], "removed");
}

/// **The leg's box keeps a neighbour list, and so does its decoupling, and neither moves a bit**:
/// the 3 Å pocket through [`leg_box`], its field's list at [`PeriodicForceField::NEIGHBOUR_SKIN`],
/// its decoupling's too; twelve steps of the leg's own dynamics at a state part-way along the
/// path go through the decoupling's list — one evaluation a step and one to start, at least one
/// build — and end on the same bits, positions and velocities, as the same box without a list.
#[test]
fn the_legs_box_keeps_a_neighbour_list() {
    let solvation = small_solvation(3.5, 6.0, 0x4E2);
    let b = pocket3();
    let with = leg_box(b, &solvation);
    let without = SolvatedComplex::new(b, &solvation).unwrap();
    assert!(without.field().neighbour_list().is_none());
    let skin = PeriodicForceField::NEIGHBOUR_SKIN;
    assert_eq!(with.field().neighbour_list().map(|l| l.skin), Some(skin));
    let run = |s: &SolvatedComplex| {
        let d = s.decoupling();
        let mut at = s.positions().to_vec();
        let mut md = s
            .dynamics(LEG_FLEXIBILITY)
            .with_bath(Bath::Langevin {
                temperature: KELVIN,
                friction: 50e12,
                seed: 5,
            })
            .thermalised(&at, KELVIN, 5);
        let potential = AtLambda {
            hamiltonian: &d,
            lambda: Lambda::new(0.5, 0.5, 1.0),
        };
        md.run(&potential, &mut at, 0.5 * FS, 12);
        (at, md.velocities().to_vec(), d.field().neighbour_list())
    };
    let (a, va, list) = run(&with);
    let (b2, vb, none) = run(&without);
    let list = list.expect("the decoupling keeps a list");
    println!("the decoupling's list: {list:?}");
    assert_eq!(list.skin, skin);
    assert_eq!(list.evaluations, 13);
    assert!(list.builds >= 1);
    assert!(none.is_none());
    let bits =
        |x: &[[f64; 3]]| -> Vec<[u64; 3]> { x.iter().map(|p| p.map(f64::to_bits)).collect() };
    assert_eq!(bits(&a), bits(&b2), "positions");
    assert_eq!(bits(&va), bits(&vb), "velocities");
}

/// **The leg's preparation is the frozen melt and then the release**, as the numbers say:
/// [`SolvatedComplex::preparation`] at [`LEG_FLEXIBILITY`] is 1000 steps of 0.5 fs at 50 ps⁻¹ and
/// 500 of 2 fs at 5 ps⁻¹ with the protein frozen, then 400 of 0.5 fs at 50 ps⁻¹ and 250 of 2 fs
/// at 5 ps⁻¹ with it free — 0.5, 1, 0.2 and 0.5 ps; a frozen protein has the melt alone, a zone
/// the release at the zone. **And [`SolvatedComplex::prepare`] runs the stages it is given**: on
/// the 3 Å pocket, four short stages with different steps, time steps and frictions are each the
/// steps they say, the protein keeps its bits through the frozen ones and moves in the others, the
/// last dynamics frees what the last stage frees, and the end is, to the bit, the same stages run
/// by hand — thermalised from the stage's seed where the flexibility changes, the velocities carried
/// where it does not.
#[test]
fn the_preparation_is_the_melt_then_the_release() {
    let stage = |flexibility, time_step, steps, friction| Stage {
        flexibility,
        time_step,
        steps,
        friction,
    };
    let melt = [
        stage(Flexibility::Frozen, 0.5e-15, 1000, 50e12),
        stage(Flexibility::Frozen, 2e-15, 500, 5e12),
    ];
    let release = |f| [stage(f, 0.5e-15, 400, 50e12), stage(f, 2e-15, 250, 5e12)];
    assert_eq!(LEG_FLEXIBILITY, Flexibility::Mobile);
    let mut want = melt.to_vec();
    want.extend(release(Flexibility::Mobile));
    let leg = SolvatedComplex::preparation(LEG_FLEXIBILITY);
    assert_eq!(leg, want);
    for (st, ps) in leg.iter().zip([0.5, 1.0, 0.2, 0.5]) {
        assert!((st.duration() / (ps * 1e-12) - 1.0).abs() < 1e-12);
    }
    assert_eq!(SolvatedComplex::preparation(Flexibility::Frozen), melt);
    let zone = Flexibility::Zone(8.0 * ANGSTROM);
    let mut zoned = melt.to_vec();
    zoned.extend(release(zone));
    assert_eq!(SolvatedComplex::preparation(zone), zoned);

    // Not the 298.15 K everything else here runs at, so that a temperature not passed on is seen.
    const WARM: f64 = 310.0;
    let s = small3(0x9E3);
    let d = s.decoupling();
    let potential = AtLambda {
        hamiltonian: &d,
        lambda: Lambda::COUPLED,
    };
    let short = [
        stage(Flexibility::Frozen, 0.5e-15, 4, 50e12),
        stage(Flexibility::Frozen, 1e-15, 3, 5e12),
        stage(Flexibility::Mobile, 0.5e-15, 3, 20e12),
        stage(Flexibility::Mobile, 1e-15, 2, 5e12),
    ];
    let protein = 0..s.binding().pocket_len();
    let mut at = s.positions().to_vec();
    let mut seen = Vec::new();
    let md = s.prepare(&potential, &short, WARM, 0x51, &mut at, |k, st, md, now| {
        assert_eq!(md.steps(), st.steps as u64, "stage {k}");
        assert_eq!(
            md.frozen(),
            s.dynamics(st.flexibility).frozen(),
            "stage {k}"
        );
        let still = protein.clone().all(|i| now[i] == s.positions()[i]);
        seen.push((k, still));
    });
    assert_eq!(seen, [(0, true), (1, true), (2, false), (3, false)]);
    assert!(md.frozen().iter().all(|&f| !f));
    let Bath::Langevin { temperature, .. } = md.bath() else {
        panic!("a bath")
    };
    assert_eq!(temperature, WARM);
    // By hand.
    let mut by_hand = s.positions().to_vec();
    let bath = |st: &Stage, k: usize| Bath::Langevin {
        temperature: WARM,
        friction: st.friction,
        seed: window_seed(0x51 ^ SolvatedComplex::PREPARATION_STREAM, k),
    };
    let mut m = s
        .dynamics(Flexibility::Frozen)
        .with_bath(bath(&short[0], 0))
        .thermalised(
            &by_hand,
            WARM,
            window_seed(0x51 ^ SolvatedComplex::PREPARATION_STREAM, 0),
        );
    m.run(&potential, &mut by_hand, 0.5e-15, 4);
    let v = m.velocities().to_vec();
    let mut m = s
        .dynamics(Flexibility::Frozen)
        .with_bath(bath(&short[1], 1))
        .with_velocities(v);
    m.run(&potential, &mut by_hand, 1e-15, 3);
    let mut m = s
        .dynamics(Flexibility::Mobile)
        .with_bath(bath(&short[2], 2))
        .thermalised(
            &by_hand,
            WARM,
            window_seed(0x51 ^ SolvatedComplex::PREPARATION_STREAM, 2),
        );
    m.run(&potential, &mut by_hand, 0.5e-15, 3);
    let v = m.velocities().to_vec();
    let mut m = s
        .dynamics(Flexibility::Mobile)
        .with_bath(bath(&short[3], 3))
        .with_velocities(v);
    m.run(&potential, &mut by_hand, 1e-15, 2);
    assert_eq!(at, by_hand, "positions to the bit");
    assert_eq!(md.velocities(), m.velocities(), "velocities to the bit");
}

/// **The leg's schedule starts the restraint in one interval**: [`SCHEDULE`] on [`grid`] is λ_r 0
/// and 1, then λ_e 0.5 and 0, then λ_v 0.9 to 0 by 0.1 — fourteen windows — and the insertion is
/// W3's, below 0.1 and at most twice.
#[test]
fn the_legs_schedule_starts_the_restraint_in_one_interval() {
    let g = grid();
    let states: Vec<[f64; 3]> = SCHEDULE.iter().map(|&k| g[k].components()).collect();
    let mut want = vec![
        [0.0, 1.0, 1.0],
        [1.0, 1.0, 1.0],
        [1.0, 0.5, 1.0],
        [1.0, 0.0, 1.0],
    ];
    for k in (0..=9).rev() {
        want.push([1.0, 0.0, f64::from(k) / 10.0]);
    }
    assert_eq!(states.len(), 14);
    for (got, want) in states.iter().zip(&want) {
        for c in 0..3 {
            assert!((got[c] - want[c]).abs() < 1e-15, "{got:?} against {want:?}");
        }
    }
    assert_eq!((OVERLAP_THRESHOLD, INSERTION_ROUNDS), (0.1, 2));
}

/// Each point's distance to the nearest of a set of atoms, by minimum image, up to `reach`
/// (anything farther reads `reach`): the atoms in cells at least `reach` wide, so that the 27
/// cells about a point hold every atom within `reach` of it.
struct Nearest {
    cell: PeriodicBox,
    reach: f64,
    counts: [usize; 3],
    cells: Vec<Vec<[f64; 3]>>,
}

impl Nearest {
    fn new(cell: PeriodicBox, atoms: &[[f64; 3]], reach: f64) -> Nearest {
        let l = cell.lengths();
        let counts = [0, 1, 2].map(|k| (l[k] / reach).floor() as usize);
        assert!(counts.iter().all(|&c| c >= 3), "three cells an axis");
        let mut cells = vec![Vec::new(); counts[0] * counts[1] * counts[2]];
        for &p in atoms {
            let i = Nearest::index(&cell, counts, p);
            cells[i].push(p);
        }
        Nearest {
            cell,
            reach,
            counts,
            cells,
        }
    }

    fn home(cell: &PeriodicBox, counts: [usize; 3], p: [f64; 3]) -> [usize; 3] {
        let w = cell.wrap(p);
        let l = cell.lengths();
        [0, 1, 2].map(|k| ((w[k] / l[k] * counts[k] as f64) as usize).min(counts[k] - 1))
    }

    fn index(cell: &PeriodicBox, counts: [usize; 3], p: [f64; 3]) -> usize {
        let h = Nearest::home(cell, counts, p);
        (h[0] * counts[1] + h[1]) * counts[2] + h[2]
    }

    fn distance(&self, p: [f64; 3]) -> f64 {
        let c = self.counts;
        let h = Nearest::home(&self.cell, c, p);
        let mut best = self.reach * self.reach;
        for dx in [c[0] - 1, 0, 1] {
            for dy in [c[1] - 1, 0, 1] {
                for dz in [c[2] - 1, 0, 1] {
                    let i = (((h[0] + dx) % c[0]) * c[1] + (h[1] + dy) % c[1]) * c[2]
                        + (h[2] + dz) % c[2];
                    for &q in &self.cells[i] {
                        let v = self.cell.minimum_image(sub(p, q));
                        best = best.min(dot(v, v));
                    }
                }
            }
        }
        best.sqrt()
    }
}

/// Water's number density, nm⁻³, at least each of `shells` from every complex atom in `at`: the
/// oxygens counted, and the volume as the share of a grid of `spacing` that far.
fn far_densities(s: &SolvatedComplex, at: &[[f64; 3]], shells: &[f64], spacing: f64) -> Vec<f64> {
    let cell = s.cell();
    let reach = shells.iter().copied().fold(0.0, f64::max);
    let near = Nearest::new(cell, &at[..s.complex_atoms()], reach);
    let l = cell.lengths();
    let n = [0, 1, 2].map(|k| (l[k] / spacing).round() as usize);
    let mut volume = vec![0usize; shells.len()];
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                let p = [
                    (i as f64 + 0.5) * l[0] / n[0] as f64,
                    (j as f64 + 0.5) * l[1] / n[1] as f64,
                    (k as f64 + 0.5) * l[2] / n[2] as f64,
                ];
                let r = near.distance(p);
                for (v, &t) in volume.iter_mut().zip(shells) {
                    *v += usize::from(r >= t);
                }
            }
        }
    }
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut oxygens = vec![0usize; shells.len()];
    for w in s.field().rigid_waters() {
        let r = near.distance(at[w[0]]);
        for (o, &t) in oxygens.iter_mut().zip(shells) {
            *o += usize::from(r >= t);
        }
    }
    oxygens
        .iter()
        .zip(&volume)
        .map(|(&o, &v)| o as f64 / (v as f64 / total * cell.volume() * 1e27))
        .collect()
}

/// **What the mobile leg costs, and the box it starts from** (release, one core, about twenty
/// minutes; run once before production):
///
/// - **the hydrogens**: the whole of 181L as placed and relaxed in vacuum ([`leg_relaxed`]), the
///   closest hydrogen contacts before and after, and Thr54's amide H against Asp47's Oδ1;
/// - **the box**: the waters each rule removes, phase 1's rule (every oxygen within 2.6 Å of a
///   heavy atom) counted beside it, and [`CREVICE_WATER`]'s fate;
/// - **the preparation** the leg runs ([`SolvatedComplex::preparation`] at [`LEG_FLEXIBILITY`]),
///   stage by stage;
/// - **the density**: then 1 ps more at 2 fs in the last stage's bath, water's number density
///   farther than 6, 8, 10 and 12 Å from every complex atom every 50 fs, against W3's 33.00 nm⁻³;
/// - **the step**: ms a step at 2 fs in the protocol's bath through the list, the fastest and
///   slowest of three runs of ten, an evaluation with and without the list, the couplings at every
///   candidate;
/// - **NVE** at 2 fs over 0.3 ps and at 1 fs over 0.3 ps from the same state: RMS, drift, ⟨T⟩ and
///   the worst held bond;
/// - **the leg**: [`SCHEDULE`]'s windows of [`PROTOCOL`], and each inserted one.
///
/// `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_cost_of_the_mobile_leg_measured -- --ignored --nocapture`.
#[test]
#[ignore = "the whole complex in water, released and timed, about twenty minutes with --release: run with --release -- --ignored --nocapture"]
fn the_cost_of_the_mobile_leg_measured() {
    let dir = results_dir();
    let name = "mobile";
    let _lock = Lock::take(&dir, name);
    let t = std::time::Instant::now();
    let placed = whole_binding();
    log(
        &dir,
        name,
        &format!("the whole binding in {:.0} s", t.elapsed().as_secs_f64()),
    );
    let t1 = std::time::Instant::now();
    let relaxed = leg_relaxed(&placed);
    let r = relaxed.hydrogen_relaxation().unwrap();
    log(
        &dir,
        name,
        &format!(
            "hydrogens relaxed in {:.1} s: complex {:?} in {} steps, protein alone {:?} in {}, \
             benzene alone {:?} in {}",
            t1.elapsed().as_secs_f64(),
            r.complex.status,
            r.complex.steps,
            r.pocket.status,
            r.pocket.steps,
            r.ligand.status,
            r.ligand.steps
        ),
    );
    for (label, b) in [("placed", &placed), ("relaxed", &relaxed)] {
        let c = closest_hydrogen_contacts(b, b.positions(), 6);
        let text: Vec<String> = c.iter().map(|&x| contact_label(b, x)).collect();
        log(&dir, name, &format!("{label}: {}", text.join("; ")));
    }
    let find = |residue: &str, atom: &str| {
        let s = the_system();
        (0..placed.positions().len())
            .find(|&k| s.residue_of(k).label() == residue && s.atom_name(k) == atom)
            .unwrap()
    };
    let (h, o) = (find("A:THR54", "H"), find("A:ASP47", "OD1"));
    log(
        &dir,
        name,
        &format!(
            "Thr54 H … Asp47 OD1: {:.3} Å placed, {:.3} Å relaxed",
            len(sub(placed.positions()[h], placed.positions()[o])) / ANGSTROM,
            len(sub(relaxed.positions()[h], relaxed.positions()[o])) / ANGSTROM
        ),
    );

    let t2 = std::time::Instant::now();
    let w4 = Solvation::w4();
    let s = leg_box(&relaxed, &w4);
    let n = s.complex_atoms();
    let lattice = WaterBox::lattice_box(s.lattice(), w4.density, w4.seed);
    let clearance_only = Solvation {
        oxygen_clearance: 0.0,
        ..w4
    }
    .removed(&lattice, &s.positions()[..n], relaxed.elements());
    let both = w4.removed(&lattice, &s.positions()[..n], relaxed.elements());
    let cell = s.cell();
    let heavy: Vec<[f64; 3]> = (0..n)
        .filter(|&k| relaxed.elements()[k] != Element::H)
        .map(|k| s.positions()[k])
        .collect();
    let phase1 = (0..lattice.count())
        .filter(|&w| {
            clearance_only[w]
                || heavy.iter().any(|&p| {
                    let v = cell.minimum_image(sub(lattice.positions()[3 * w], p));
                    len(v) < w4.oxygen_clearance
                })
        })
        .count();
    let count = |m: &[bool]| m.iter().filter(|&&x| x).count();
    // 18 419 g/mol at 0.73 cm³/g, in waters at 33.00 nm⁻³.
    let volume = 18_419.0 / 6.022_140_76e23 * 0.73 * 1e-6 * 33.0e27;
    log(
        &dir,
        name,
        &format!(
            "the box in {:.0} s: {} atoms, {} waters, {} Cl-; of {} lattice waters the clearance \
             removes {}, the oxygen between two heavy atoms {} more, {} in all (phase 1's rule \
             would remove {}; the protein's volume at 0.73 cm³/g is {volume:.0}); the crevice \
             water {}: clearance {}, either rule {}",
            t2.elapsed().as_secs_f64(),
            s.positions().len(),
            s.waters(),
            s.ions(),
            lattice.count(),
            count(&clearance_only),
            count(&both) - count(&clearance_only),
            s.overlapping(),
            phase1,
            CREVICE_WATER,
            clearance_only[CREVICE_WATER],
            both[CREVICE_WATER]
        ),
    );
    assert_eq!(count(&both), s.overlapping());

    let d = s.decoupling();
    let potential = AtLambda {
        hamiltonian: &d,
        lambda: Lambda::COUPLED,
    };
    let mut at = s.positions().to_vec();
    let t3 = std::time::Instant::now();
    let mut md = s.prepare(
        &potential,
        &SolvatedComplex::preparation(LEG_FLEXIBILITY),
        KELVIN,
        PREPARATION_SEED,
        &mut at,
        |k, st, md, _| {
            log(
                &dir,
                name,
                &format!(
                    "stage {k}: {} steps of {} fs, {:?}, {} ps⁻¹: {:.0} s in all, T {:.1} K, U \
                     {:.1} kcal/mol",
                    st.steps,
                    st.time_step / FS,
                    st.flexibility,
                    st.friction / 1e12,
                    t3.elapsed().as_secs_f64(),
                    md.temperature(),
                    kcal(md.potential_energy().unwrap())
                ),
            )
        },
    );

    // The density.
    let shells = [
        6.0 * ANGSTROM,
        8.0 * ANGSTROM,
        10.0 * ANGSTROM,
        12.0 * ANGSTROM,
    ];
    let mut series: Vec<Vec<f64>> = vec![Vec::new(); shells.len()];
    let t4 = std::time::Instant::now();
    for _ in 0..20 {
        md.run(&potential, &mut at, 2.0 * FS, 25);
        for (s2, x) in series
            .iter_mut()
            .zip(far_densities(&s, &at, &shells, ANGSTROM))
        {
            s2.push(x);
        }
    }
    for (k, x) in series.iter().enumerate() {
        let e = Estimate::of(x);
        let half = |r: std::ops::Range<usize>| x[r.clone()].iter().sum::<f64>() / r.len() as f64;
        log(
            &dir,
            name,
            &format!(
                "water farther than {:.0} Å from the complex: {:.3} ± {:.3} nm⁻³ ({:+.2}% of \
                 33.00), halves {:.3} and {:.3}; n {}, τ {:.2} samples",
                shells[k] / ANGSTROM,
                e.mean,
                e.error,
                100.0 * (e.mean / 33.0 - 1.0),
                half(0..10),
                half(10..20),
                e.samples,
                e.tau
            ),
        );
    }
    log(
        &dir,
        name,
        &format!("the density's 1 ps in {:.0} s", t4.elapsed().as_secs_f64()),
    );

    // The step.
    let mut md = md.with_bath(Bath::Langevin {
        temperature: KELVIN,
        friction: PROTOCOL.friction,
        seed: 0x5EEF,
    });
    let mut chunks = Vec::new();
    for _ in 0..3 {
        let t = std::time::Instant::now();
        md.run(&potential, &mut at, 2.0 * FS, 10);
        chunks.push(t.elapsed().as_secs_f64() / 10.0);
    }
    let fastest = chunks.iter().copied().fold(f64::INFINITY, f64::min);
    let slowest = chunks.iter().copied().fold(0.0, f64::max);
    let time = |f: &dyn Fn(), k: usize| {
        let mut best = f64::INFINITY;
        for _ in 0..k {
            let t = std::time::Instant::now();
            f();
            best = best.min(t.elapsed().as_secs_f64());
        }
        best
    };
    let listed = time(&|| drop(d.field().evaluate(&at)), 5);
    let cells = d.field().clone().without_neighbour_list();
    let unlisted = time(&|| drop(cells.evaluate(&at)), 3);
    let states = grid();
    let couplings = time(&|| drop(d.couplings(&at, &states)), 3);
    log(
        &dir,
        name,
        &format!(
            "{:.1}–{:.1} ms a step at 2 fs, {} degrees of freedom; an evaluation {:.1} ms through \
             the list, {:.1} by the cell list; the couplings at {} candidates {:.1} ms; the list {:?}",
            fastest * 1e3,
            slowest * 1e3,
            md.degrees_of_freedom(),
            listed * 1e3,
            unlisted * 1e3,
            states.len(),
            couplings * 1e3,
            d.field().neighbour_list()
        ),
    );

    // NVE.
    let mobile = s.mobile(LEG_FLEXIBILITY);
    let template = s.dynamics(LEG_FLEXIBILITY);
    let v = md.velocities().to_vec();
    for (dt, steps) in [(2.0 * FS, 150), (1.0 * FS, 300)] {
        let (rms, drift, temperature, bond) =
            nve(&s, &potential, &template, &mobile, &at, &v, dt, steps);
        log(
            &dir,
            name,
            &format!(
                "NVE at {} fs over 0.3 ps: RMS {:.4} kcal/mol, drift {:+.4} kcal/mol/ps, ⟨T⟩ {:.1} \
                 K, worst held bond {:.1e}",
                dt / FS,
                kcal(rms),
                kcal(drift) * 1e-12,
                temperature,
                bond
            ),
        );
    }

    // The leg.
    let per_window = |step: f64| {
        (PROTOCOL.steps_per_window() as f64 * step + PROTOCOL.samples as f64 * couplings) / 3600.0
    };
    let windows = SCHEDULE.len() as f64;
    log(
        &dir,
        name,
        &format!(
            "the leg: {} windows of {} steps ({:.0} ps each), a window {:.1}–{:.1} h, {:.0}–{:.0} h in \
             all before any is inserted, {:.1}–{:.1} h each inserted; the preparation {:.0} min",
            SCHEDULE.len(),
            PROTOCOL.steps_per_window(),
            PROTOCOL.steps_per_window() as f64 * PROTOCOL.time_step / 1e-12,
            per_window(fastest),
            per_window(slowest),
            windows * per_window(fastest),
            windows * per_window(slowest),
            per_window(fastest),
            per_window(slowest),
            t3.elapsed().as_secs_f64() / 60.0
        ),
    );
    log(
        &dir,
        name,
        &format!("the measurement took {:.0} s", t.elapsed().as_secs_f64()),
    );
}

// ---------------------------------------------------------------------------------------------
// The complex leg: written, not run

/// What one window recorded, read back from its file.
#[derive(Clone, Debug)]
struct Recorded {
    grid: usize,
    lambda: Lambda,
    gradients: Vec<[f64; 3]>,
    /// Per sample, the coupling energy at every candidate, joules.
    energies: Vec<Vec<f64>>,
    seconds: f64,
}

impl Recorded {
    /// Reduced `u(to) − u(own)` per sample.
    fn delta(&self, to: usize, kt: f64) -> Vec<f64> {
        self.energies
            .iter()
            .map(|e| (e[to] - e[self.grid]) / kt)
            .collect()
    }

    fn slice(&self, range: std::ops::Range<usize>) -> Recorded {
        Recorded {
            gradients: self.gradients[range.clone()].to_vec(),
            energies: self.energies[range].to_vec(),
            ..self.clone()
        }
    }
}

/// BAR between neighbouring windows, reduced, the delta method's variance.
fn interval_bar(a: &Recorded, b: &Recorded, kt: f64) -> pantometry_forcefield::Bar {
    let reverse: Vec<f64> = b.delta(a.grid, kt).iter().map(|x| -x).collect();
    pantometry_forcefield::free_energy::bar_correlated(&a.delta(b.grid, kt), &reverse)
}

/// BAR along `records`, reduced, by the delta method window by window.
fn chain(records: &[Recorded], kt: f64) -> pantometry_forcefield::FreeEnergy {
    let n = records.len();
    let windows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|j| {
            let prev = if j > 0 {
                records[j].delta(records[j - 1].grid, kt)
            } else {
                Vec::new()
            };
            let next = if j + 1 < n {
                records[j].delta(records[j + 1].grid, kt)
            } else {
                Vec::new()
            };
            (prev, next)
        })
        .collect();
    bennett_chain(&windows)
}

/// TI on the trapezoid along `records`' path, kcal/mol, with its error.
fn trapezoid(records: &[Recorded]) -> (f64, f64) {
    let n = records.len();
    let terms: Vec<Estimate> = (0..n)
        .map(|k| {
            let lo = records[k.saturating_sub(1)].lambda.components();
            let hi = records[(k + 1).min(n - 1)].lambda.components();
            let w = [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]));
            let ys: Vec<f64> = records[k]
                .gradients
                .iter()
                .map(|g| (0..3).map(|c| g[c] * w[c]).sum())
                .collect();
            Estimate::of(&ys)
        })
        .collect();
    (
        kcal(terms.iter().map(|t| t.mean).sum()),
        kcal(terms.iter().map(|t| t.error * t.error).sum::<f64>().sqrt()),
    )
}

/// One window of the leg, as `Windows` would run it at schedule position `index`: the template's
/// dynamics in the protocol's bath, thermalised from `window_seed(protocol.seed, index)`, at
/// `lambda`, from `start`. After the equilibration, every `stride` steps, `on` is given the step
/// and the couplings at `lambda` followed by every one of `states`. Returns where it ended. W3's,
/// on [`sample_frames`].
#[allow(clippy::too_many_arguments)]
fn sample_window(
    h: &PeriodicDecoupling,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    p: Protocol,
    index: usize,
    lambda: Lambda,
    states: &[Lambda],
    mut on: impl FnMut(u64, &[Coupling]),
) -> Vec<[f64; 3]> {
    let mut lambdas = vec![lambda];
    lambdas.extend_from_slice(states);
    sample_frames(h, start, template, p, index, lambda, |s, at| {
        on(s, &h.couplings(at, &lambdas))
    })
}

/// The frames a window samples: the template's dynamics in the protocol's bath, thermalised from
/// `window_seed(protocol.seed, index)`, at `lambda`, from `start`; after the equilibration, every
/// `stride` steps, `on` is given the step and the positions. Returns where it ended. What
/// [`sample_window`] records couplings from, and what the decomposition takes its frames from, so
/// that both are the leg's own dynamics.
fn sample_frames(
    h: &PeriodicDecoupling,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    p: Protocol,
    index: usize,
    lambda: Lambda,
    on: impl FnMut(u64, &[[f64; 3]]),
) -> Vec<[f64; 3]> {
    let potential = AtLambda {
        hamiltonian: h,
        lambda,
    };
    sample_on(&potential, start, template, p, index, on)
}

/// [`sample_frames`] on any `potential`: the template's dynamics in the protocol's bath,
/// thermalised from `window_seed(protocol.seed, index)`, from `start`; after the equilibration,
/// every `stride` steps, `on` is given the step and the positions. What the staged states of
/// variant 4 run on ([`Staged`]), so that a state at `s` = 0 is a window's dynamics to the bit.
fn sample_on(
    potential: &dyn pantometry_forcefield::Potential,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    p: Protocol,
    index: usize,
    mut on: impl FnMut(u64, &[[f64; 3]]),
) -> Vec<[f64; 3]> {
    let seed = window_seed(p.seed, index);
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: p.temperature,
            friction: p.friction,
            seed,
        })
        .thermalised(start, p.temperature, seed);
    let mut at = start.to_vec();
    md.prepare(potential, &at);
    let mut taken = 0;
    // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
    // builds on 1.78.
    #[allow(clippy::manual_is_multiple_of)]
    while taken < p.samples {
        md.step(potential, &mut at, p.time_step);
        let s = md.steps();
        if s > p.equilibration && (s - p.equilibration) % p.stride == 0 {
            on(s, &at);
            taken += 1;
        }
    }
    at
}

/// Window `g`'s header: the candidate and its λ, the protocol, the atom count, the flexibility,
/// the number of candidates, **and digests of the leg's start and of its Boresch restraint**
/// (anchors, references and force constants), so that a window run from another preparation, or
/// another restraint, is refused rather than mixed with these.
fn header(leg: &Leg, g: usize) -> String {
    let p = leg.protocol;
    let l = leg.grid[g];
    let r = leg.restraint;
    let restraint = digest(
        r.receptor
            .iter()
            .chain(&r.ligand)
            .map(|&i| i as f64)
            .chain(r.reference)
            .chain(r.force_constants),
    );
    format!(
        "# candidate {g} lambda {} {} {} protocol {} {} {} {} {} {} {} atoms {} flexibility {:?} \
         states {} start {:016x} restraint {restraint:016x}",
        exact(l.restraint),
        exact(l.electrostatics),
        exact(l.van_der_waals),
        exact(p.time_step),
        exact(p.temperature),
        exact(p.friction),
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        leg.start.len(),
        leg.flexibility,
        leg.grid.len(),
        positions_digest(&leg.start)
    )
}

fn parse(text: &str) -> Recorded {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().expect("a header").split_whitespace().collect();
    let field = |name: &str| head.iter().position(|&w| w == name).expect(name) + 1;
    let grid_index: usize = head[field("candidate")].parse().unwrap();
    let l = field("lambda");
    let lambda = Lambda::new(
        head[l].parse().unwrap(),
        head[l + 1].parse().unwrap(),
        head[l + 2].parse().unwrap(),
    );
    let states: usize = head[field("states")].parse().unwrap();
    let (mut gradients, mut energies, mut seconds) = (Vec::new(), Vec::new(), f64::NAN);
    for line in lines {
        if let Some(rest) = line.strip_prefix("# done seconds ") {
            seconds = rest.trim().parse().unwrap();
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(v.len(), 3 + states, "a sample line");
        gradients.push([v[0], v[1], v[2]]);
        energies.push(v[3..].to_vec());
    }
    assert!(seconds.is_finite(), "a complete window says so");
    Recorded {
        grid: grid_index,
        lambda,
        gradients,
        energies,
        seconds,
    }
}

/// The leg's inputs: its decoupling with the restraint on, the dynamics every window starts from,
/// the prepared start, the protocol and the candidates. The real leg's are [`PROTOCOL`] and
/// [`grid`]; the resume test's are a few steps on a pocket.
struct Leg {
    decoupling: PeriodicDecoupling,
    template: MolecularDynamics,
    start: Vec<[f64; 3]>,
    flexibility: Flexibility,
    restraint: Boresch,
    protocol: Protocol,
    grid: Vec<Lambda>,
    /// The log a window's progress goes to, `NAME.log` in the results directory: `complex` for
    /// the leg. **Not in the header**: where a window logs changes no bit of what it records.
    log: &'static str,
}

/// The leg from a prepared `start`: the Boresch restraint chosen there by 3c's rule, on the box's
/// decoupling, and the leg's flexibility.
fn leg_from(
    s: &SolvatedComplex,
    start: Vec<[f64; 3]>,
    protocol: Protocol,
    grid: Vec<Lambda>,
) -> Leg {
    let restraint = anchors(s.binding(), &start[..s.complex_atoms()]);
    Leg {
        decoupling: s.decoupling().with_restraint(restraint),
        template: s.dynamics(LEG_FLEXIBILITY),
        start,
        flexibility: LEG_FLEXIBILITY,
        restraint,
        protocol,
        grid,
        log: "complex",
    }
}

/// Window `g`'s record: **from its file if a run completed it under the same header**
/// ([`header`]), or run now — written sample by sample to `window_NN.partial`, which is begun
/// afresh whatever a killed run left there, and renamed to `window_NN.txt` when complete. A
/// finished window is never run again; a partial one is never read. W3's.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn window(leg: &Leg, dir: &std::path::Path, g: usize) -> Recorded {
    use std::io::Write;
    let path = dir.join(format!("window_{g:02}.txt"));
    let grid = &leg.grid;
    let p = leg.protocol;
    let want = header(leg, g);
    if let Ok(text) = std::fs::read_to_string(&path) {
        assert_eq!(
            text.lines().next().unwrap_or(""),
            want,
            "{} was run differently",
            path.display()
        );
        let r = parse(&text);
        assert_eq!(r.energies.len(), p.samples);
        return r;
    }
    let t = std::time::Instant::now();
    let partial = path.with_extension("partial");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
    writeln!(file, "{want}").unwrap();
    let total = p.steps_per_window();
    log(
        dir,
        leg.log,
        &format!(
            "window at candidate {g} λ = {:?} begins: {total} steps",
            grid[g]
        ),
    );
    let mut taken = 0usize;
    sample_window(
        &leg.decoupling,
        &leg.start,
        &leg.template,
        p,
        g,
        grid[g],
        grid,
        |s, c| {
            let mut line = format!(
                "{s} {} {} {}",
                exact(c[0].gradient[0]),
                exact(c[0].gradient[1]),
                exact(c[0].gradient[2])
            );
            for x in &c[1..] {
                line.push(' ');
                line.push_str(&exact(x.energy));
            }
            writeln!(file, "{line}").unwrap();
            taken += 1;
            if taken % 100 == 0 {
                file.flush().unwrap();
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    leg.log,
                    &format!(
                        "candidate {g}: {taken} of {} samples, step {s} of {total}, {secs:.0} s, \
                         {:.1} ms per step",
                        p.samples,
                        secs * 1e3 / s as f64
                    ),
                );
            }
        },
    );
    let seconds = t.elapsed().as_secs_f64();
    writeln!(file, "# done seconds {}", exact(seconds)).unwrap();
    file.flush().unwrap();
    drop(file);
    std::fs::rename(&partial, &path).expect("rename");
    log(
        dir,
        leg.log,
        &format!("window at candidate {g} done in {seconds:.0} s"),
    );
    parse(&std::fs::read_to_string(&path).unwrap())
}

/// A fresh directory under the system's temporary one for a test that writes files, named by the
/// process and `tag`.
#[cfg(not(target_family = "wasm"))]
fn fresh_dir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("pantometry-w4-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("a temporary directory");
    d
}

/// **Positions round-trip to the bit**: forty-two numbers — the smallest subnormal and the
/// largest, the smallest normal, `±0.0`, `±f64::MAX`, ε, a third, a tenth and thirty more from
/// arbitrary bits (no NaN, no infinity) — through [`exact`] alone and through
/// [`write_positions`] and [`read_positions`], header and all, `to_bits` equal, `−0.0` kept apart
/// from `0.0`; nothing is left at the `.partial` path; and **the file cut short at every byte**
/// short of its last newline is refused, whatever it then parses as.
#[cfg(not(target_family = "wasm"))]
#[test]
fn positions_round_trip_to_the_bit() {
    let mut xs = vec![
        f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::MIN_POSITIVE,
        -0.0,
        0.0,
        f64::MAX,
        -f64::MAX,
        f64::EPSILON,
        1.0 / 3.0,
        0.1,
        -2.5e-10,
        1e-300,
    ];
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    while xs.len() < 42 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let x = f64::from_bits(state);
        if x.is_finite() {
            xs.push(x);
        }
    }
    for &x in &xs {
        let back: f64 = exact(x).parse().unwrap();
        assert_eq!(back.to_bits(), x.to_bits(), "{x:e}");
    }
    let at: Vec<[f64; 3]> = xs.chunks(3).map(|c| [c[0], c[1], c[2]]).collect();
    let dir = fresh_dir("round-trip");
    let path = dir.join("positions.txt");
    write_positions(&path, "# a header", &at);
    assert!(!path.with_extension("partial").exists());
    let (head, back) = read_positions(&path).expect("written").expect("whole");
    assert_eq!(head, "# a header");
    let bits =
        |v: &[[f64; 3]]| -> Vec<[u64; 3]> { v.iter().map(|p| p.map(f64::to_bits)).collect() };
    assert_eq!(bits(&back), bits(&at));
    let text = std::fs::read_to_string(&path).unwrap();
    for cut in (0..text.len() - 1).filter(|&k| text.is_char_boundary(k)) {
        assert!(
            parse_positions(&text[..cut]).is_err(),
            "cut at byte {cut} of {} read as whole",
            text.len()
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// What a panic said.
#[cfg(not(target_family = "wasm"))]
fn refusal(f: impl FnOnce()) -> String {
    let e = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).expect_err("refused");
    e.downcast_ref::<String>()
        .cloned()
        .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
        .unwrap_or_default()
}

/// **A killed leg resumes to the same bits**, through the leg's own [`prepared`], [`leg_from`] and
/// [`window`] on the 3 Å pocket, two stages of preparation and windows of eight 0.5 fs steps and
/// three samples, at candidates 0 and 8, in temporary directories:
///
/// - **uninterrupted**: preparation and both windows;
/// - **killed and resumed**: preparation and window 0, then what a run killed during window 8
///   leaves — its `.partial` file, a header and one sample — then the leg again. Window 0 is
///   read, not run (its file's bytes, the run's `seconds` among them, unchanged), `prepared.txt`
///   is read, the partial is begun afresh and renamed, and **every gradient and energy of both
///   records is the uninterrupted run's to the bit**.
///
/// And each guard refuses by name: the preparation asked for with another seed, a window asked for
/// from a start one ulp away, and a `prepared.txt` cut short.
#[cfg(not(target_family = "wasm"))]
#[test]
fn a_killed_leg_resumes_to_the_same_bits() {
    // The restraint's nine candidates: windows 0 and 8 are its ends.
    let restraint_grid = || grid()[..9].to_vec();
    let s = small3(0x2E5);
    let coupled = s.decoupling();
    let potential = AtLambda {
        hamiltonian: &coupled,
        lambda: Lambda::COUPLED,
    };
    let stages = [
        Stage {
            flexibility: Flexibility::Frozen,
            time_step: 0.5e-15,
            steps: 3,
            friction: 50e12,
        },
        Stage {
            flexibility: LEG_FLEXIBILITY,
            time_step: 0.5e-15,
            steps: 2,
            friction: 50e12,
        },
    ];
    let protocol = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: 0x2E6,
    };
    let start = |dir: &std::path::Path, seed: u64| {
        prepared(
            &s,
            &potential,
            &stages,
            seed,
            KELVIN,
            dir,
            "prepared.txt",
            "complex",
        )
    };
    let run = |dir: &std::path::Path, schedule: &[usize]| -> Vec<Recorded> {
        let leg = leg_from(&s, start(dir, PREPARATION_SEED), protocol, restraint_grid());
        schedule.iter().map(|&g| window(&leg, dir, g)).collect()
    };
    let whole = fresh_dir("whole");
    let full = run(&whole, &[0, 8]);
    let killed = fresh_dir("killed");
    run(&killed, &[0]);
    let bytes = |p: std::path::PathBuf| std::fs::read(p).expect("a file");
    let (finished, prepared_bytes) = (
        bytes(killed.join("window_00.txt")),
        bytes(killed.join("prepared.txt")),
    );
    let text = std::fs::read_to_string(whole.join("window_08.txt")).unwrap();
    let partial: String = text.lines().take(2).map(|l| format!("{l}\n")).collect();
    std::fs::write(killed.join("window_08.partial"), partial).unwrap();
    let resumed = run(&killed, &[0, 8]);
    // Compared whole, not printed: a difference would dump the file.
    assert!(
        bytes(killed.join("window_00.txt")) == finished,
        "window 0 read, not run"
    );
    assert!(bytes(killed.join("prepared.txt")) == prepared_bytes);
    assert!(!killed.join("window_08.partial").exists());
    let f64s = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
    assert_eq!(full.len(), resumed.len());
    for (a, b) in full.iter().zip(&resumed) {
        assert_eq!(a.grid, b.grid);
        assert_eq!(
            a.lambda.components().map(f64::to_bits),
            b.lambda.components().map(f64::to_bits)
        );
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
        "resumed: {} records, {} numbers each the uninterrupted run's",
        resumed.len(),
        resumed[0].gradients.len() * 3 + resumed[0].energies.len() * restraint_grid().len()
    );

    let why = refusal(|| drop(start(&killed, PREPARATION_SEED + 1)));
    assert!(why.contains("prepared differently"), "{why}");
    let mut moved = start(&killed, PREPARATION_SEED);
    moved[0][0] = f64::from_bits(moved[0][0].to_bits() + 1);
    let leg = leg_from(&s, moved, protocol, restraint_grid());
    let why = refusal(|| drop(window(&leg, &killed, 0)));
    assert!(why.contains("was run differently"), "{why}");
    let cut = fresh_dir("cut");
    std::fs::write(
        cut.join("prepared.txt"),
        &prepared_bytes[..prepared_bytes.len() - 30],
    )
    .unwrap();
    let why = refusal(|| drop(start(&cut, PREPARATION_SEED)));
    assert!(why.contains("not whole"), "{why}");
    for d in [whole, killed, cut] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// The leg's prepared start, from [`whole_complex`]: [`prepared`] with the coupled box, the
/// preparation at [`LEG_FLEXIBILITY`], [`PREPARATION_SEED`] and 298.15 K, kept in `prepared.txt` in
/// `dir` and read back, header checked, when a run has written it. The one place the leg's start
/// is asked for: the leg, the refined window and the decomposition all call it.
fn production_start(s: &SolvatedComplex, dir: &std::path::Path) -> Vec<[f64; 3]> {
    let coupled = s.decoupling();
    prepared(
        s,
        &AtLambda {
            hamiltonian: &coupled,
            lambda: Lambda::COUPLED,
        },
        &SolvatedComplex::preparation(LEG_FLEXIBILITY),
        PREPARATION_SEED,
        KELVIN,
        dir,
        "prepared.txt",
        "complex",
    )
}

/// The leg as [`the_complex_leg_measured`] runs it from `start`: [`leg_from`] with [`PROTOCOL`]
/// and [`grid`]. A window of any leg this returns has the header the leg's own would have.
fn production_leg(s: &SolvatedComplex, start: Vec<[f64; 3]>) -> Leg {
    leg_from(s, start, PROTOCOL, grid())
}

/// **Benzene decoupled from T4 lysozyme L99A in TIP3P, the complex leg, and the binding free
/// energy it closes with W3's solvent leg.** **Run once**, at `db8a90f`, from 2026-10-07 to
/// 2026-10-13, exit 0 after 539 298 s on one core, with no window inserted: the complex leg BAR
/// +18.362 ± 0.258 kcal/mol, and `ΔG°_bind` = +1.570 − (+18.362) − (−7.749) = **−9.04 ± 0.29**
/// against experiment's −5.19 ± 0.16 (the CHANGELOG has the intervals). Nothing here asserts it.
/// Its cost was measured first (`the_cost_of_the_mobile_leg_measured`). The whole of 181L as
/// [`whole_complex`] builds it — the hydrogens relaxed in vacuum first, the neighbour list on — **the protein free**
/// ([`LEG_FLEXIBILITY`]), prepared by [`SolvatedComplex::preparation`]: the frozen melt and then
/// the release, 0.2 ps at 0.5 fs and 0.5 ps at 2 fs, before any window (kept in `prepared.txt`
/// and read back if there). The Boresch restraint chosen by 3c's rule at that start, [`SCHEDULE`]
/// on [`grid`] with W3's [`PROTOCOL`], and a window added at the candidate between any two whose
/// overlap is below [`OVERLAP_THRESHOLD`], [`INSERTION_ROUNDS`] rounds at most
/// ([`pantometry_forcefield::free_energy::refine_schedule`]). Each window is written to
/// `window_NN.txt` in [`results_dir`] as it runs and read back if there, so a stopped run resumes;
/// progress is in `complex.log`. Printed: each insertion, each interval, the three segments, BAR
/// and TI, the halves, the release, and `ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release` with
/// W3's solvent leg. `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_complex_leg_measured -- --ignored --nocapture`.
#[test]
#[ignore = "the complex leg: six days on one core with --release, run once (2026-10-13)"]
fn the_complex_leg_measured() {
    let flexibility = LEG_FLEXIBILITY;
    let dir = results_dir();
    let lock = Lock::take(&dir, "complex");
    log(
        &dir,
        "complex",
        &format!(
            "the leg begins; its lock is {} — remove it by hand only once no run is writing here",
            lock.0.display()
        ),
    );
    let t = std::time::Instant::now();
    let s = whole_complex();
    let leg = production_leg(&s, production_start(&s, &dir));
    let restraint = leg.restraint;
    let kt = BOLTZMANN.to_si() * KELVIN;
    let release = kcal(restraint.release_free_energy(KELVIN));
    log(
        &dir,
        "complex",
        &format!(
            "{} atoms, flexibility {flexibility:?}; Boresch: receptor {:?}, ligand {:?}; r₀ {:.3} \
             Å, θ_A {:.1}°, θ_B {:.1}°; release to 1 M {release:+.3} kcal/mol (extended {:+.3})",
            leg.start.len(),
            restraint.receptor,
            restraint.ligand,
            restraint.reference[0] / ANGSTROM,
            restraint.reference[1].to_degrees(),
            restraint.reference[2].to_degrees(),
            kcal(restraint.release_free_energy_extended(KELVIN))
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
            "complex",
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
    let to_kcal = |x: f64| kcal(kt * x);
    for k in 0..records.len() - 1 {
        let bar = interval_bar(&records[k], &records[k + 1], kt);
        log(
            &dir,
            "complex",
            &format!(
                "| {k} | {:?} | {:+.3} ± {:.3} | {:.3} |",
                records[k].lambda,
                to_kcal(bar.delta),
                to_kcal(bar.error()),
                bar.overlap
            ),
        );
    }
    let last = |p: &dyn Fn(&Lambda) -> bool| records.iter().rposition(|r| p(&r.lambda)).unwrap();
    let on = last(&|l| l.electrostatics == 1.0 && l.van_der_waals == 1.0);
    let off = last(&|l| l.van_der_waals == 1.0);
    for (name, range) in [
        ("restraint on", 0..on + 1),
        ("charges off", on..off + 1),
        ("van der Waals off", off..records.len()),
    ] {
        let c = chain(&records[range.clone()], kt);
        let (ti, te) = trapezoid(&records[range]);
        log(
            &dir,
            "complex",
            &format!(
                "{name}: BAR {:+.3} ± {:.3}, TI {ti:+.3} ± {te:.3} kcal/mol",
                to_kcal(c.value),
                to_kcal(c.error)
            ),
        );
    }
    let total = chain(&records, kt);
    let (bar, bar_error) = (to_kcal(total.value), to_kcal(total.error));
    let (ti, ti_error) = trapezoid(&records);
    let half = PROTOCOL.samples / 2;
    let first: Vec<Recorded> = records.iter().map(|r| r.slice(0..half)).collect();
    let second: Vec<Recorded> = records
        .iter()
        .map(|r| r.slice(half..PROTOCOL.samples))
        .collect();
    let (h1, h2) = (chain(&first, kt), chain(&second, kt));
    let bind = SOLVENT_LEG.0 - bar - release;
    let error = (SOLVENT_LEG.1 * SOLVENT_LEG.1 + bar_error * bar_error).sqrt();
    log(
        &dir,
        "complex",
        &format!(
            "complex leg: BAR {bar:+.3} ± {bar_error:.3}, TI {ti:+.3} ± {ti_error:.3} kcal/mol; \
             halves {:+.3} and {:+.3}; ΔG°_bind = {:+.3} − ({bar:+.3}) − ({release:+.3}) = \
             {bind:+.3} ± {error:.3} kcal/mol; experiment −5.19 ± 0.16; {:.0} s of windows, {:.0} \
             s in all",
            to_kcal(h1.value),
            to_kcal(h2.value),
            SOLVENT_LEG.0,
            records.iter().map(|r| r.seconds).sum::<f64>(),
            t.elapsed().as_secs_f64()
        ),
    );
}

// ---------------------------------------------------------------------------------------------
// Where the 3.85 kcal/mol comes from: benzene's interaction with what surrounds it, taken apart

/// What an atom around benzene is, for the decomposition: a protein atom by its UFF type, an atom
/// of a water, or an ion.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Kind {
    Protein(UffType),
    Water,
    Ion,
}

impl Kind {
    fn label(self) -> String {
        match self {
            Kind::Protein(t) => format!("protein {}", t.label()),
            Kind::Water => "water".into(),
            Kind::Ion => "ions".into(),
        }
    }

    /// Its class, an index into [`CLASSES`].
    fn class(self) -> usize {
        match self {
            Kind::Protein(_) => 0,
            Kind::Water => 1,
            Kind::Ion => 2,
        }
    }
}

/// The classes the electrostatics are split into. **The reciprocal sum is not a sum over pairs**,
/// so the cross Coulomb energy is split by zeroing every charge outside a class: the cross term is
/// linear in the rest's charges, its real-space pairs and its background alike, so the classes'
/// sum to the whole's to rounding.
const CLASSES: [&str; 3] = ["protein", "water", "ions"];

/// **Benzene in a box, set up to be taken apart**: every atom's kind and its van der Waals `[x, D]`
/// — **rebuilt here** from the binding's UFF types, TIP3P's oxygen and the chloride, not read from
/// the force field, so that the identity with the field's own coupling checks the parameters as
/// well as the sums — and the decoupling the coupling is read from, with one more per class with
/// every charge outside it zeroed.
struct Surroundings {
    cell: PeriodicBox,
    cutoff: f64,
    ligand: std::ops::Range<usize>,
    /// Each atom's kind; `None` for benzene's.
    kinds: Vec<Option<Kind>>,
    vdw: Vec<[f64; 2]>,
    /// Benzene's carbons: its ring.
    ring: Vec<usize>,
    /// The protein's heavy atoms.
    heavy: Vec<bool>,
    whole: PeriodicDecoupling,
    /// [`CLASSES`]' decouplings, `None` for a class with no atom.
    classes: Vec<Option<PeriodicDecoupling>>,
}

/// Benzene's van der Waals with one kind of atom: the pairs inside the cutoff, the long-range
/// correction's share, and how many pairs.
#[derive(Clone, Copy, Debug, Default)]
struct Share {
    pairs: f64,
    tail: f64,
    count: usize,
}

/// One frame's van der Waals, kind by kind, with `Σ |term|` and the number of terms for the
/// rounding bound.
struct VanDerWaals {
    by_kind: std::collections::BTreeMap<Kind, Share>,
    magnitude: f64,
    terms: usize,
}

impl VanDerWaals {
    fn total(&self) -> f64 {
        self.by_kind.values().map(|s| s.pairs + s.tail).sum()
    }

    fn class(&self, c: usize) -> f64 {
        self.by_kind
            .iter()
            .filter(|(k, _)| k.class() == c)
            .map(|(_, s)| s.pairs + s.tail)
            .sum()
    }
}

/// The decomposition held to the field: each sum's discrepancy over its bound, so that a ratio
/// above 1 is a failure.
#[derive(Clone, Copy, Debug)]
struct Identity {
    /// `Σ` kinds against the couplings' `U(λ_v = 1) − U(λ_v = 0)`, both at `λ_e` = 0.
    van_der_waals: f64,
    /// `Σ` classes against the couplings' cross Coulomb energy.
    coulomb: f64,
    /// Both against the whole evaluation's `U(coupled) − U(decoupled)`, when asked for.
    evaluation: Option<f64>,
}

impl Identity {
    fn worst(&self) -> f64 {
        self.van_der_waals
            .max(self.coulomb)
            .max(self.evaluation.unwrap_or(0.0))
    }
}

impl Surroundings {
    fn new(
        whole: PeriodicDecoupling,
        ligand: std::ops::Range<usize>,
        kinds: Vec<Option<Kind>>,
        vdw: Vec<[f64; 2]>,
        ring: Vec<usize>,
        heavy: Vec<bool>,
    ) -> Surroundings {
        let n = kinds.len();
        assert_eq!(whole.group().len(), n);
        assert_eq!(vdw.len(), n);
        for i in 0..n {
            assert_eq!(
                whole.group()[i],
                kinds[i].is_none(),
                "atom {i}: benzene or not"
            );
            assert_eq!(
                whole.group()[i],
                ligand.contains(&i),
                "atom {i}: the ligand"
            );
        }
        // No pair across is excluded, so every pair inside the cutoff counts.
        assert!(whole
            .field()
            .excluded_pairs()
            .iter()
            .all(|&[i, j]| whole.group()[i] == whole.group()[j]));
        let field = whole.field();
        let (cell, cutoff) = (field.cell(), field.ewald().parameters().cutoff);
        let classes = (0..CLASSES.len())
            .map(|c| {
                if !kinds.iter().any(|k| k.is_some_and(|k| k.class() == c)) {
                    return None;
                }
                let q: Vec<f64> = (0..n)
                    .map(|i| match kinds[i] {
                        Some(k) if k.class() != c => 0.0,
                        _ => field.charges()[i],
                    })
                    .collect();
                Some(
                    PeriodicDecoupling::new(&field.clone().with_charges(q), whole.group())
                        .expect("benzene decouples"),
                )
            })
            .collect();
        Surroundings {
            cell,
            cutoff,
            ligand,
            kinds,
            vdw,
            ring,
            heavy,
            whole,
            classes,
        }
    }

    /// Benzene's `[x, D]` as the field has them.
    fn uff_benzene(&self) -> Vec<[f64; 2]> {
        self.vdw[self.ligand.clone()].to_vec()
    }

    /// The minimum-image separation of atoms `a < b` from wrapped positions `w`, as the field's
    /// cross terms take it — a comparison per axis, exact — and its square.
    fn separation(&self, w: &[[f64; 3]], a: usize, b: usize) -> ([f64; 3], f64) {
        let l = self.cell.lengths();
        let mut d = sub(w[a], w[b]);
        for x in 0..3 {
            if d[x] > 0.5 * l[x] {
                d[x] -= l[x];
            } else if d[x] < -0.5 * l[x] {
                d[x] += l[x];
            }
        }
        (d, dot(d, d))
    }

    /// Benzene's van der Waals with everything else at `at`, kind by kind, **with benzene's own
    /// `[x, D]` replaced by `benzene`** and every other atom's as built, combined by UFF's
    /// geometric rule: every pair inside the cutoff, and each kind's share of the long-range
    /// correction, `2 (2π/V) [S₆ᵏ S₆ᴮ/(9 r_c⁹) − 2 S₃ᵏ S₃ᴮ/(3 r_c³)]` with `S_m = Σ √D x^m`.
    fn van_der_waals(&self, at: &[[f64; 3]], benzene: &[[f64; 2]]) -> VanDerWaals {
        assert_eq!(benzene.len(), self.ligand.len());
        let param = |k: usize| {
            if self.ligand.contains(&k) {
                benzene[k - self.ligand.start]
            } else {
                self.vdw[k]
            }
        };
        let w: Vec<[f64; 3]> = at.iter().map(|p| self.cell.wrap(*p)).collect();
        let rc2 = self.cutoff * self.cutoff;
        let mut by_kind: std::collections::BTreeMap<Kind, Share> = Default::default();
        let (mut magnitude, mut terms) = (0.0, 0usize);
        for g in self.ligand.clone() {
            for (j, kind) in self.kinds.iter().enumerate() {
                let Some(kind) = *kind else { continue };
                let (a, b) = (g.min(j), g.max(j));
                let (_, r2) = self.separation(&w, a, b);
                if r2 >= rc2 {
                    continue;
                }
                let (pa, pb) = (param(a), param(b));
                let pair = Pair {
                    atoms: [a, b],
                    distance: (pa[0] * pb[0]).sqrt(),
                    well: (pa[1] * pb[1]).sqrt(),
                };
                let u = pair.energy(r2.sqrt());
                let share = by_kind.entry(kind).or_default();
                share.pairs += u;
                share.count += 1;
                magnitude += u.abs();
                terms += 1;
            }
        }
        let sums = |p: [f64; 2]| {
            let x3 = p[0] * p[0] * p[0];
            [p[1].sqrt() * x3 * x3, p[1].sqrt() * x3]
        };
        let mut group = [0.0; 2];
        for p in benzene {
            let s = sums(*p);
            group[0] += s[0];
            group[1] += s[1];
        }
        let mut rest: std::collections::BTreeMap<Kind, [f64; 2]> = Default::default();
        for (j, kind) in self.kinds.iter().enumerate() {
            if let Some(kind) = kind {
                let s = sums(self.vdw[j]);
                let e = rest.entry(*kind).or_default();
                e[0] += s[0];
                e[1] += s[1];
            }
        }
        let rc3 = self.cutoff * self.cutoff * self.cutoff;
        let rc9 = rc3 * rc3 * rc3;
        let volume = self.cell.volume();
        for (kind, s) in rest {
            let tail = 2.0 * 2.0 * std::f64::consts::PI / volume
                * (s[0] * group[0] / (9.0 * rc9) - 2.0 * s[1] * group[1] / (3.0 * rc3));
            by_kind.entry(kind).or_default().tail = tail;
            magnitude += tail.abs();
            terms += 1;
        }
        VanDerWaals {
            by_kind,
            magnitude,
            terms,
        }
    }

    /// Benzene's cross Coulomb energy with each class at `at`, joules: each class's decoupling's
    /// own, at full strength.
    fn coulomb(&self, at: &[[f64; 3]]) -> [f64; 3] {
        let mut c = [0.0; 3];
        for (k, d) in self.classes.iter().enumerate() {
            if let Some(d) = d {
                c[k] = d.couplings(at, &[Lambda::COUPLED])[0].gradient[1];
            }
        }
        c
    }

    /// The field's own: benzene's whole van der Waals with the rest, `U(0, 0, 1) − U(0, 0, 0)` of
    /// the couplings, and the cross Coulomb energy.
    fn coupling(&self, at: &[[f64; 3]]) -> (f64, f64) {
        let c = self.whole.couplings(
            at,
            &[Lambda::new(0.0, 0.0, 1.0), Lambda::new(0.0, 0.0, 0.0)],
        );
        (c[0].energy - c[1].energy, c[0].gradient[1])
    }

    /// **The parts held to the whole**: the kinds' van der Waals against the couplings' to the
    /// summation bound `2 (n + 1) ε Σ|term|` — the two sums add the same `n` terms in different
    /// orders, the long-range correction split by kind among them — and the classes' Coulomb to
    /// `64 ε` of the whole evaluation's parts (W1's and W3's bound for the cross terms through the
    /// mesh, measured against the evaluation at the frame given). With `evaluation`, both against
    /// the whole evaluation's `U(coupled) − U(decoupled)` as well, to `64 ε` of both evaluations'
    /// parts: another path through the field, the pair search and the mesh's potentials.
    fn identity(
        &self,
        at: &[[f64; 3]],
        vdw: &VanDerWaals,
        coulomb: [f64; 3],
        scale: f64,
        evaluation: bool,
    ) -> Identity {
        let (whole_vdw, whole_coulomb) = self.coupling(at);
        let vdw_bound = 2.0 * (vdw.terms as f64 + 1.0) * EPS * vdw.magnitude;
        let parts_coulomb: f64 = coulomb.iter().sum();
        let c_bound = 64.0 * EPS * scale;
        let evaluation = evaluation.then(|| {
            let on = self.whole.energy(at, Lambda::COUPLED);
            let off = self.whole.energy(at, Lambda::new(0.0, 0.0, 0.0));
            let bound = 64.0 * EPS * (magnitude(&on) + magnitude(&off)) + vdw_bound;
            ((on.total - off.total) - (vdw.total() + parts_coulomb)).abs() / bound
        });
        Identity {
            van_der_waals: (vdw.total() - whole_vdw).abs() / vdw_bound,
            coulomb: (parts_coulomb - whole_coulomb).abs() / c_bound,
            evaluation,
        }
    }

    /// The scale of the Coulomb bound at `at`: the whole evaluation's parts.
    fn scale(&self, at: &[[f64; 3]]) -> f64 {
        magnitude(&self.whole.energy(at, Lambda::COUPLED))
    }

    /// The atoms of each kind within `radius` of a ring carbon, and the protein's heavy atoms
    /// among them.
    fn contacts(
        &self,
        at: &[[f64; 3]],
        radius: f64,
    ) -> (std::collections::BTreeMap<Kind, usize>, usize) {
        let w: Vec<[f64; 3]> = at.iter().map(|p| self.cell.wrap(*p)).collect();
        let mut near: std::collections::BTreeMap<Kind, usize> = Default::default();
        let mut heavy = 0;
        for (j, kind) in self.kinds.iter().enumerate() {
            let Some(kind) = *kind else { continue };
            let close = self.ring.iter().any(|&g| {
                let (a, b) = (g.min(j), g.max(j));
                self.separation(&w, a, b).1 < radius * radius
            });
            if close {
                *near.entry(kind).or_default() += 1;
                heavy += usize::from(self.heavy[j]);
            }
        }
        (near, heavy)
    }
}

/// The ligand's carbons in the complex's order.
fn ring_of(elements: &[Element], ligand: std::ops::Range<usize>) -> Vec<usize> {
    ligand.filter(|&i| elements[i] == Element::C).collect()
}

/// A solvated complex set up to be taken apart: the protein by its UFF types, the waters, the
/// ions.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn complex_surroundings(s: &SolvatedComplex) -> Surroundings {
    let n = s.positions().len();
    let nc = s.complex_atoms();
    let first_ion = n - s.ions();
    let b = s.binding();
    let ligand = s.ligand();
    let ion = Ion::chloride();
    let oxygen = water::oxygen_vdw();
    assert!(
        s.field().charges()[first_ion..]
            .iter()
            .all(|&q| q == ion.charge),
        "the ions are chlorides"
    );
    let mut kinds = Vec::with_capacity(n);
    let mut vdw = Vec::with_capacity(n);
    for i in 0..n {
        if i < nc {
            let t = b.types()[i];
            let p = t.parameters();
            vdw.push([p.vdw_distance, p.vdw_energy]);
            kinds.push((!ligand.contains(&i)).then_some(Kind::Protein(t)));
        } else if i < first_ion {
            vdw.push(if (i - nc) % 3 == 0 {
                oxygen
            } else {
                [0.0, 0.0]
            });
            kinds.push(Some(Kind::Water));
        } else {
            vdw.push([ion.vdw_distance, ion.vdw_energy]);
            kinds.push(Some(Kind::Ion));
        }
    }
    let heavy = (0..n)
        .map(|i| i < nc && !ligand.contains(&i) && b.elements()[i] != Element::H)
        .collect();
    Surroundings::new(
        s.decoupling(),
        ligand.clone(),
        kinds,
        vdw,
        ring_of(b.elements(), ligand),
        heavy,
    )
}

/// W3's benzene in its water, set up the same way.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn water_surroundings(b: &hydration::Benzene, w: &hydration::Solvated) -> Surroundings {
    let nb = b.at.len();
    let n = w.at.len();
    let oxygen = water::oxygen_vdw();
    let mut kinds = Vec::with_capacity(n);
    let mut vdw = Vec::with_capacity(n);
    for (i, t) in b.types.iter().enumerate() {
        let p = t.parameters();
        vdw.push([p.vdw_distance, p.vdw_energy]);
        kinds.push(None);
        assert!(i < nb);
    }
    for i in nb..n {
        vdw.push(if (i - nb) % 3 == 0 {
            oxygen
        } else {
            [0.0, 0.0]
        });
        kinds.push(Some(Kind::Water));
    }
    let elements: Vec<Element> = b.component.atoms().iter().map(|a| a.element).collect();
    Surroundings::new(
        w.decoupling(),
        0..nb,
        kinds,
        vdw,
        ring_of(&elements, 0..nb),
        vec![false; n],
    )
}

/// **GAFF's benzene van der Waals**, `[x, D]` in UFF's form: `x = 2 R*/2`, `D = ε`, for `ca` on
/// each carbon and `ha` on each hydrogen. Read from `gaff.dat` (GAFF 1.7, Nov 2013) in ambermini
/// (github.com/choderalab/ambermini, commit `f7c421c46dca9f0e`, `share/amber/dat/leap/parm/gaff.dat`,
/// SHA-256 `95e3c655…44d8433`), its `MOD4 RE` block: `ca 1.9080 0.0860 OPLS` and
/// `ha 1.4590 0.0150 Spellmeyer`; and the same in OpenMM's `gaff.xml` (openmmforcefields, commit
/// `3f73b6ff`, from gaff-1.81.dat): `ca` σ 0.3399669508 nm, ε 0.359824 kJ/mol, and `ha` σ
/// 0.2599642460 nm, ε 0.06276 kJ/mol, which are the same numbers through `R*/2 = 2^⅙ σ/2` and
/// 4.184. **Combined here by UFF's geometric rule with everything else, which GAFF is not**: GAFF
/// combines by Lorentz–Berthelot against GAFF and AMBER atoms, so this is a hybrid, benzene's
/// size and depth GAFF's and the mixing and the environment UFF's and TIP3P's.
fn gaff_benzene(elements: &[Element]) -> Vec<[f64; 2]> {
    elements
        .iter()
        .map(|e| match e {
            Element::C => [2.0 * 1.9080 * ANGSTROM, 0.0860 * KCAL_PER_MOL],
            Element::H => [2.0 * 1.4590 * ANGSTROM, 0.0150 * KCAL_PER_MOL],
            other => panic!("benzene has no {other:?}"),
        })
        .collect()
}

/// [`Surroundings`] at `at` held to the field: the kinds' van der Waals to the couplings' and the
/// classes' Coulomb to the cross Coulomb energy, both to the whole evaluation's `U(coupled) −
/// U(decoupled)` too, each within its bound ([`Surroundings::identity`]). Every kind carries a share
/// of the long-range correction and has pairs exactly when its energy from them is not zero; water
/// has pairs; every class with an atom has a Coulomb energy of its own; and GAFF's benzene
/// ([`gaff_benzene`]) moves every kind's pair energy. Returns how many protein types have pairs.
fn surroundings_sum_to_the_coupling(name: &str, sur: &Surroundings, at: &[[f64; 3]]) -> usize {
    let vdw = sur.van_der_waals(at, &sur.uff_benzene());
    let coulomb = sur.coulomb(at);
    let id = sur.identity(at, &vdw, coulomb, sur.scale(at), true);
    let (whole_vdw, whole_c) = sur.coupling(at);
    println!(
        "{name}: van der Waals {:+.4} kcal/mol in {} terms, {:+.2e} of it off the coupling's \
         ({:.1e} of the bound); Coulomb {:+.4} ({:.1e} of the bound); against the evaluation \
         {:.1e} of the bound",
        kcal(whole_vdw),
        vdw.terms,
        kcal(vdw.total() - whole_vdw),
        id.van_der_waals,
        kcal(whole_c),
        id.coulomb,
        id.evaluation.unwrap()
    );
    for (k, share) in &vdw.by_kind {
        println!(
            "  {:<13} {:+.5} kcal/mol, {} pairs, tail {:+.2e}",
            k.label(),
            kcal(share.pairs + share.tail),
            share.count,
            kcal(share.tail)
        );
        assert!(share.tail != 0.0, "{name}: {k:?}'s tail");
        assert_eq!(share.count == 0, share.pairs == 0.0, "{name}: {k:?}");
    }
    for c in 0..CLASSES.len() {
        if sur.classes[c].is_some() {
            println!(
                "  {:<13} Coulomb {:+.5} kcal/mol",
                CLASSES[c],
                kcal(coulomb[c])
            );
            assert_ne!(coulomb[c], 0.0, "{name}: {}", CLASSES[c]);
        }
    }
    assert!(
        vdw.by_kind[&Kind::Water].count > 0,
        "{name}: water within the cutoff"
    );
    assert!(id.worst() <= 1.0, "{name}: {id:?}");
    let gaff = sur.van_der_waals(at, &gaff_benzene(&ring_and_hydrogens(sur)));
    for (k, share) in vdw.by_kind.iter().filter(|(_, s)| s.count > 0) {
        assert_ne!(
            share.pairs, gaff.by_kind[k].pairs,
            "{name}: {k:?} under GAFF"
        );
    }
    vdw.by_kind
        .iter()
        .filter(|(k, s)| matches!(k, Kind::Protein(_)) && s.count > 0)
        .count()
}

/// **Benzene's surroundings sum to its coupling**: on the 3 Å pocket in water as built, by its five
/// protein types or more, and on W3's small water box, [`surroundings_sum_to_the_coupling`]. The
/// pocket is neutral and has no ion: [`an_ions_share_of_benzenes_coupling_is_its_own`] is that
/// branch.
#[test]
fn benzenes_surroundings_sum_to_its_coupling() {
    let s = small3(0x5EC0);
    let sur = complex_surroundings(&s);
    assert!(sur.classes[2].is_none(), "the 3 Å pocket is neutral");
    let types = surroundings_sum_to_the_coupling("3 Å pocket", &sur, s.positions());
    assert!(types >= 5, "{types} protein types");
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC2, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    surroundings_sum_to_the_coupling("W3's small box", &water_surroundings(&b, &w), &w.at);
}

/// **An ion's share of benzene's coupling is its own**: the 5 Å pocket with its chloride, as built,
/// [`surroundings_sum_to_the_coupling`], the ion a kind and a class of its own.
#[test]
fn an_ions_share_of_benzenes_coupling_is_its_own() {
    let s = small5(0x5EC1);
    let sur = complex_surroundings(&s);
    assert!(sur.classes[2].is_some(), "the 5 Å pocket has its chloride");
    let types = surroundings_sum_to_the_coupling("5 Å pocket", &sur, s.positions());
    assert!(types >= 5, "{types} protein types");
}

/// Benzene's elements as a surroundings' ligand has them: carbons on its ring, hydrogens else.
fn ring_and_hydrogens(sur: &Surroundings) -> Vec<Element> {
    sur.ligand
        .clone()
        .map(|i| {
            if sur.ring.contains(&i) {
                Element::C
            } else {
                Element::H
            }
        })
        .collect()
}

/// The decomposition's frames: W3's and the leg's 2 fs in a 1 ps⁻¹ bath at 298.15 K, 4 ps
/// discarded and then 320 frames 50 fs apart, 16 ps; seeded apart from every window.
const FRAMES: Protocol = Protocol {
    time_step: 2e-15,
    temperature: KELVIN,
    friction: 1e12,
    equilibration: 2_000,
    stride: 25,
    samples: 320,
    seed: 0x4_4DEC,
};

/// A protein atom is in contact with benzene within this of a ring carbon.
const CONTACT: f64 = 4.5 * ANGSTROM;

/// One frame, taken apart.
struct Taken {
    uff: VanDerWaals,
    gaff: VanDerWaals,
    coulomb: [f64; 3],
    near: std::collections::BTreeMap<Kind, usize>,
    heavy: usize,
}

/// `frames` frames of `sur`'s box from `start` through [`sample_frames`] on `h` at `lambda` with
/// `template`, each taken apart under UFF's benzene and GAFF's, and held to the field's coupling;
/// every twentieth frame also against the whole evaluation. Returns the frames, the worst
/// identity, and the milliseconds a step.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn take_apart(
    sur: &Surroundings,
    h: &PeriodicDecoupling,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    lambda: Lambda,
    name: &str,
) -> (Vec<Taken>, Identity, f64) {
    let uff = sur.uff_benzene();
    let gaff = gaff_benzene(&ring_and_hydrogens(sur));
    let mut taken: Vec<Taken> = Vec::new();
    let mut worst = Identity {
        van_der_waals: 0.0,
        coulomb: 0.0,
        evaluation: Some(0.0),
    };
    let mut scale = sur.scale(start);
    let t = std::time::Instant::now();
    let mut last = (0u64, 0.0f64);
    sample_frames(h, start, template, FRAMES, 0, lambda, |step, at| {
        let k = taken.len();
        let full = k % 20 == 0;
        if full {
            scale = sur.scale(at);
        }
        let v = sur.van_der_waals(at, &uff);
        let g = sur.van_der_waals(at, &gaff);
        let c = sur.coulomb(at);
        let id = sur.identity(at, &v, c, scale, full);
        worst.van_der_waals = worst.van_der_waals.max(id.van_der_waals);
        worst.coulomb = worst.coulomb.max(id.coulomb);
        if let (Some(w), Some(e)) = (&mut worst.evaluation, id.evaluation) {
            *w = w.max(e);
        }
        let (near, heavy) = sur.contacts(at, CONTACT);
        if (k + 1) % 40 == 0 {
            let secs = t.elapsed().as_secs_f64();
            println!(
                "{name}: frame {} at step {step}, {secs:.0} s, {:.1} ms a step since the last \
                 report; vdW {:+.3} kcal/mol, {heavy} heavy atoms within 4.5 Å",
                k + 1,
                (secs - last.1) * 1e3 / (step - last.0) as f64,
                kcal(v.total())
            );
            last = (step, secs);
        }
        taken.push(Taken {
            uff: v,
            gaff: g,
            coulomb: c,
            near,
            heavy,
        });
    });
    let per_step = t.elapsed().as_secs_f64() * 1e3 / FRAMES.steps_per_window() as f64;
    (taken, worst, per_step)
}

/// `Estimate::of` in kcal/mol from joules.
fn estimate_kcal(xs: impl Iterator<Item = f64>) -> Estimate {
    let v: Vec<f64> = xs.map(kcal).collect();
    Estimate::of(&v)
}

fn pm(e: &Estimate) -> String {
    format!("{:+.3} ± {:.3} (τ {:.1})", e.mean, e.error, e.tau)
}

/// The free energy of turning benzene's van der Waals into GAFF's on `taken`, by Zwanzig's
/// `ΔG = −kT ln ⟨exp(−βΔU)⟩`, kcal/mol, with what says whether to trust it: the error of the mean
/// weight carried through the logarithm, the spread of `βΔU`, the second cumulant `⟨ΔU⟩ − β σ²/2`,
/// the effective sample size `(Σw)²/Σw²` of the weights, and the same on each half.
struct Zwanzig {
    mean: Estimate,
    exp: f64,
    exp_error: f64,
    spread: f64,
    cumulant: f64,
    ess: f64,
    halves: [f64; 2],
}

fn zwanzig(du: &[f64], kt: f64) -> Zwanzig {
    let exp_of = |x: &[f64]| -> (f64, f64, f64) {
        let b: Vec<f64> = x.iter().map(|u| u / kt).collect();
        let m = b.iter().copied().fold(f64::INFINITY, f64::min);
        let w: Vec<f64> = b.iter().map(|y| (-(y - m)).exp()).collect();
        let e = Estimate::of(&w);
        let sum: f64 = w.iter().sum();
        let ess = sum * sum / w.iter().map(|v| v * v).sum::<f64>();
        (
            kcal(kt * (m - e.mean.ln())),
            kcal(kt * e.error / e.mean),
            ess,
        )
    };
    let (exp, exp_error, ess) = exp_of(du);
    let half = du.len() / 2;
    let mean = estimate_kcal(du.iter().copied());
    let n = du.len() as f64;
    let mu = du.iter().sum::<f64>() / n;
    let var = du.iter().map(|u| (u - mu) * (u - mu)).sum::<f64>() / n;
    Zwanzig {
        mean,
        exp,
        exp_error,
        spread: var.sqrt() / kt,
        cumulant: kcal(mu - var / (2.0 * kt)),
        ess,
        halves: [exp_of(&du[..half]).0, exp_of(&du[half..]).0],
    }
}

/// **Where the binding's 3.85 kcal/mol comes from** (release, one core; 10 000 steps of the
/// complex, about seventy minutes at 300 ms a step and three hours at the 1 000 its one run had):
/// frames of benzene bound and in water, its interaction with what surrounds it taken apart, and
/// the same frames with GAFF's benzene van der Waals.
///
/// - **Frames.** The leg's own: [`production_start`] reads `prepared.txt` back through its header
///   check (refused if absent: nothing is written to the results directory), and
///   [`sample_frames`] runs the leg's dynamics and Hamiltonian at window 0's state — coupled, the
///   restraint off — from it, [`FRAMES`]; then W3's box ([`hydration::w3_solvated`]) melted as W3
///   melts it ([`hydration::w3_melt`]), through the same [`sample_frames`] and [`FRAMES`].
/// - **Taken apart** ([`Surroundings`]): benzene's van der Waals with each protein UFF type, the
///   waters and the ions, its cross Coulomb with each class, and the protein's heavy atoms within
///   4.5 Å of its ring; means with [`Estimate`]'s autocorrelation-aware errors. Window 0's 2000
///   samples give the whole van der Waals and Coulomb over 200 ps beside them.
/// - **GAFF's benzene** ([`gaff_benzene`], a hybrid): `⟨ΔU⟩` bound and in water,
///   `ΔΔG_bind ≈ ⟨ΔU⟩_complex − ⟨ΔU⟩_water` to first order with `ΔU = U_GAFF − U_UFF`, and
///   Zwanzig's estimate with its spread and effective sample size.
/// - **Every frame held to the field**: the parts' sums against the coupling, and every twentieth
///   against the whole evaluation ([`Surroundings::identity`]), asserted after everything is
///   printed.
///
/// `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p where_the_binding_error_comes_from_measured -- --ignored --nocapture`.
#[test]
#[ignore = "frames of the whole complex and of W3's box, about seventy minutes with --release: run with --release -- --ignored --nocapture"]
fn where_the_binding_error_comes_from_measured() {
    let dir = results_dir();
    for f in ["prepared.txt", "window_00.txt"] {
        assert!(
            dir.join(f).exists(),
            "{} is not there: this measurement reads the leg's results and writes none",
            dir.join(f).display()
        );
    }
    let t = std::time::Instant::now();
    let kt = BOLTZMANN.to_si() * KELVIN;
    let s = whole_complex();
    let leg = production_leg(&s, production_start(&s, &dir));
    println!(
        "the leg's box and start read back in {:.0} s",
        t.elapsed().as_secs_f64()
    );

    // Window 0 at the same state: the whole van der Waals (candidate 12, λ = (1, 0, 1), less the
    // restraint) and the cross Coulomb energy, over its 2000 samples.
    let g = grid();
    assert_eq!(g[0], Lambda::new(0.0, 1.0, 1.0));
    assert_eq!(g[12], Lambda::new(1.0, 0.0, 1.0));
    let w0 = window(&leg, &dir, 0);
    let w0_vdw = estimate_kcal(
        w0.energies
            .iter()
            .zip(&w0.gradients)
            .map(|(e, d)| e[12] - d[0]),
    );
    let w0_c = estimate_kcal(w0.gradients.iter().map(|d| d[1]));

    let bound = complex_surroundings(&s);
    let (complex, complex_id, complex_ms) = take_apart(
        &bound,
        &leg.decoupling,
        &leg.start,
        &leg.template,
        g[0],
        "complex",
    );
    if complex_ms > 400.0 {
        println!("**{complex_ms:.0} ms a step on the complex, above 400: the run was slowed**");
    }

    let benzene = hydration::Benzene::new(Some(pocket3().ligand_force_field().charges().to_vec()));
    assert_eq!(
        &benzene.types[..],
        &pocket3().types()[pocket3().ligand_range()]
    );
    let w3 = hydration::w3_solvated(&benzene);
    let w3_decoupling = w3.decoupling();
    let w3_template = w3.dynamics(&benzene);
    let mut w3_start = w3.at.clone();
    let tm = std::time::Instant::now();
    hydration::w3_melt(&w3_decoupling, &w3_template, KELVIN, &mut w3_start);
    println!(
        "W3's box, {} waters, melted in {:.0} s",
        w3.waters,
        tm.elapsed().as_secs_f64()
    );
    let free = water_surroundings(&benzene, &w3);
    let (water, water_id, water_ms) = take_apart(
        &free,
        &w3_decoupling,
        &w3_start,
        &w3_template,
        Lambda::COUPLED,
        "water",
    );

    // The tables.
    println!(
        "\n{} frames each, {} fs apart after {} ps; {complex_ms:.0} ms a step bound, {water_ms:.1} \
         in water",
        FRAMES.samples,
        FRAMES.stride as f64 * FRAMES.time_step / FS,
        FRAMES.equilibration as f64 * FRAMES.time_step * 1e12
    );
    println!(
        "\n| benzene's van der Waals with | bound (kcal/mol) | in water | atoms within 4.5 Å |"
    );
    println!("| --- | --- | --- | --- |");
    let kinds: std::collections::BTreeSet<Kind> = complex
        .iter()
        .chain(&water)
        .flat_map(|f| f.uff.by_kind.keys().copied())
        .collect();
    let of = |frames: &[Taken], f: &dyn Fn(&VanDerWaals) -> f64| {
        estimate_kcal(frames.iter().map(|x| f(&x.uff)))
    };
    let near_of = |frames: &[Taken], k: Kind| {
        let v: Vec<f64> = frames
            .iter()
            .map(|x| x.near.get(&k).copied().unwrap_or(0) as f64)
            .collect();
        Estimate::of(&v)
    };
    for &k in &kinds {
        let share = |v: &VanDerWaals| v.by_kind.get(&k).map_or(0.0, |s| s.pairs + s.tail);
        let b = of(&complex, &share);
        let n = near_of(&complex, k);
        let w = (k == Kind::Water).then(|| of(&water, &share));
        println!(
            "| {} | {} | {} | {:.1} ± {:.1} |",
            k.label(),
            pm(&b),
            w.map_or("".into(), |w| pm(&w)),
            n.mean,
            n.error
        );
    }
    let total_b = estimate_kcal(complex.iter().map(|x| x.uff.total()));
    let total_w = estimate_kcal(water.iter().map(|x| x.uff.total()));
    println!("| **all** | {} | {} | |", pm(&total_b), pm(&total_w));
    for c in 0..CLASSES.len() {
        let b = estimate_kcal(complex.iter().map(|x| x.uff.class(c)));
        println!("| {} (class) | {} | | |", CLASSES[c], pm(&b));
    }
    let heavy = Estimate::of(&complex.iter().map(|x| x.heavy as f64).collect::<Vec<_>>());
    println!(
        "\nprotein heavy atoms within 4.5 Å of a ring carbon: {:.2} ± {:.2} (τ {:.1}), from {} to {}",
        heavy.mean,
        heavy.error,
        heavy.tau,
        complex.iter().map(|x| x.heavy).min().unwrap(),
        complex.iter().map(|x| x.heavy).max().unwrap()
    );
    println!(
        "window 0 over 200 ps: van der Waals {}, Coulomb {}",
        pm(&w0_vdw),
        pm(&w0_c)
    );

    println!("\n| benzene's cross Coulomb with | bound (kcal/mol) | in water |");
    println!("| --- | --- | --- |");
    for c in 0..CLASSES.len() {
        let b = estimate_kcal(complex.iter().map(|x| x.coulomb[c]));
        let w = (c == 1).then(|| estimate_kcal(water.iter().map(|x| x.coulomb[c])));
        println!(
            "| {} | {} | {} |",
            CLASSES[c],
            pm(&b),
            w.map_or("".into(), |w| pm(&w))
        );
    }
    let cb = estimate_kcal(complex.iter().map(|x| x.coulomb.iter().sum()));
    let cw = estimate_kcal(water.iter().map(|x| x.coulomb.iter().sum()));
    println!("| **all** | {} | {} |", pm(&cb), pm(&cw));

    // The pairs, UFF's against the hybrid's.
    println!("\n| benzene atom with | UFF x_ij (Å), D_ij (kcal/mol) | hybrid x_ij, D_ij |");
    println!("| --- | --- | --- |");
    let gaff = gaff_benzene(&[Element::C, Element::H]);
    let uff_c = UffType::CR.parameters();
    let uff_h = UffType::H.parameters();
    let mut seen = std::collections::BTreeSet::new();
    for (j, kind) in bound.kinds.iter().enumerate() {
        let Some(kind) = *kind else { continue };
        let p = bound.vdw[j];
        if p[1] == 0.0 || !seen.insert(kind) {
            continue;
        }
        for (name, u, h) in [
            ("C", [uff_c.vdw_distance, uff_c.vdw_energy], gaff[0]),
            ("H", [uff_h.vdw_distance, uff_h.vdw_energy], gaff[1]),
        ] {
            println!(
                "| {name} with {} | {:.4}, {:.4} | {:.4}, {:.4} |",
                kind.label(),
                (u[0] * p[0]).sqrt() / ANGSTROM,
                kcal((u[1] * p[1]).sqrt()),
                (h[0] * p[0]).sqrt() / ANGSTROM,
                kcal((h[1] * p[1]).sqrt())
            );
        }
    }

    // GAFF's benzene on the same frames.
    let du = |frames: &[Taken]| -> Vec<f64> {
        frames
            .iter()
            .map(|x| x.gaff.total() - x.uff.total())
            .collect()
    };
    let (zb, zw) = (zwanzig(&du(&complex), kt), zwanzig(&du(&water), kt));
    println!("\n| ΔU = U_hybrid − U_UFF | bound | in water |");
    println!("| --- | --- | --- |");
    for c in 0..2 {
        let b = estimate_kcal(complex.iter().map(|x| x.gaff.class(c) - x.uff.class(c)));
        let w =
            (c == 1).then(|| estimate_kcal(water.iter().map(|x| x.gaff.class(c) - x.uff.class(c))));
        println!(
            "| ⟨ΔU⟩ with {} | {} | {} |",
            CLASSES[c],
            pm(&b),
            w.map_or("".into(), |w| pm(&w))
        );
    }
    for (name, f) in [
        (
            "⟨ΔU⟩",
            &(|z: &Zwanzig| pm(&z.mean)) as &dyn Fn(&Zwanzig) -> String,
        ),
        ("σ(βΔU)", &|z: &Zwanzig| format!("{:.3}", z.spread)),
        ("⟨ΔU⟩ − βσ²/2", &|z: &Zwanzig| {
            format!("{:+.3}", z.cumulant)
        }),
        ("Zwanzig −kT ln⟨e^(−βΔU)⟩", &|z: &Zwanzig| {
            format!("{:+.3} ± {:.3}", z.exp, z.exp_error)
        }),
        ("Zwanzig, halves", &|z: &Zwanzig| {
            format!("{:+.3}, {:+.3}", z.halves[0], z.halves[1])
        }),
        ("effective samples", &|z: &Zwanzig| {
            format!("{:.1} of {}", z.ess, FRAMES.samples)
        }),
    ] {
        println!("| {name} | {} | {} |", f(&zb), f(&zw));
    }
    let first = zb.mean.mean - zw.mean.mean;
    let first_error = (zb.mean.error * zb.mean.error + zw.mean.error * zw.mean.error).sqrt();
    let exp = zb.exp - zw.exp;
    let exp_error = (zb.exp_error * zb.exp_error + zw.exp_error * zw.exp_error).sqrt();
    println!(
        "\nΔΔG_bind = ΔG_bind(hybrid) − ΔG_bind(UFF) = ΔG_swap(bound) − ΔG_swap(water): first order \
         ⟨ΔU⟩_bound − ⟨ΔU⟩_water = {first:+.3} ± {first_error:.3}; second cumulant {:+.3}; \
         Zwanzig {exp:+.3} ± {exp_error:.3} kcal/mol. The binding is {:.2} too strong against \
         experiment; a positive ΔΔG is GAFF's benzene binding less.",
        zb.cumulant - zw.cumulant,
        -9.04 - -5.19
    );

    println!(
        "\nthe identity, worst over the frames (1 is the bound): complex van der Waals {:.3}, \
         Coulomb {:.3}, against the evaluation {:.3}; water {:.3}, {:.3}, {:.3}. {:.0} s in all",
        complex_id.van_der_waals,
        complex_id.coulomb,
        complex_id.evaluation.unwrap(),
        water_id.van_der_waals,
        water_id.coulomb,
        water_id.evaluation.unwrap(),
        t.elapsed().as_secs_f64()
    );
    assert!(complex_id.worst() <= 1.0, "{complex_id:?}");
    assert!(water_id.worst() <= 1.0, "{water_id:?}");
}

// ---------------------------------------------------------------------------------------------
// The rest of the 3.79 kcal/mol: benzene's cross van der Waals reparametrised, on saved frames

/// How two atoms' `[x, D]` combine into their pair's.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Mixing {
    /// UFF's, eqs 21b and 22: `x_ij = √(x_i x_j)` and `D_ij = √(D_i D_j)`. The field's.
    Geometric,
    /// Lorentz–Berthelot, AMBER's and GAFF's: `σ_ij = (σ_i + σ_j)/2` and `ε_ij = √(ε_i ε_j)`. With
    /// `x = 2^⅙ σ`, linear in σ, that is `x_ij = (x_i + x_j)/2` — AMBER's `R*_ij = R*_i/2 + R*_j/2`
    /// — and `D_ij = √(D_i D_j)`.
    LorentzBerthelot,
}

impl Mixing {
    /// The pair's `[x_ij, D_ij]`.
    fn pair(self, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        match self {
            Mixing::Geometric => [(a[0] * b[0]).sqrt(), (a[1] * b[1]).sqrt()],
            Mixing::LorentzBerthelot => [0.5 * (a[0] + b[0]), (a[1] * b[1]).sqrt()],
        }
    }
}

/// **Benzene's cross van der Waals under one model**: every atom's `[x, D]`, benzene's and the
/// rest's, and the rule that combines a pair. Only the pairs across are in it; benzene's own pairs
/// stay UFF's in every variant, so they cancel from every `ΔU` and from every `ΔΔG`.
#[derive(Clone, Debug)]
struct CrossModel {
    vdw: Vec<[f64; 2]>,
    mixing: Mixing,
}

/// The variants, in the order [`variants`] builds them. **0 is 5794556's hybrid**, measured again
/// on the saved frames, which are its frames to the bit if nothing it ran on has changed: its
/// `ΔΔG` was +1.711 ± 0.056.
const VARIANTS: [&str; 4] = [
    "0: GAFF benzene, geometric (5794556's hybrid)",
    "1: UFF benzene, Lorentz–Berthelot",
    "2: GAFF benzene, Lorentz–Berthelot",
    "3: GAFF benzene, ff14SB protein, Lorentz–Berthelot",
];

/// What each atom around benzene is called in the decomposition: its class (a protein atom's
/// ff14SB type, `water`, `ions`) and its residue (`A:LEU84`, `water`, `ions`), and a protein atom's
/// ff14SB `[x, D]`. `None` for benzene's atoms.
#[derive(Clone, Debug)]
struct Labels {
    classes: Vec<String>,
    class: Vec<Option<usize>>,
    residues: Vec<String>,
    residue: Vec<Option<usize>>,
    amber: Vec<Option<[f64; 2]>>,
}

impl Labels {
    fn new() -> Labels {
        Labels {
            classes: Vec::new(),
            class: Vec::new(),
            residues: Vec::new(),
            residue: Vec::new(),
            amber: Vec::new(),
        }
    }

    /// The next atom: benzene's when `class` is `None`.
    fn push(&mut self, class: Option<(&str, &str)>, amber: Option<[f64; 2]>) {
        let index = |names: &mut Vec<String>, name: &str| {
            names.iter().position(|n| n == name).unwrap_or_else(|| {
                names.push(name.to_string());
                names.len() - 1
            })
        };
        match class {
            None => {
                self.class.push(None);
                self.residue.push(None);
            }
            Some((c, r)) => {
                let c = index(&mut self.classes, c);
                let r = index(&mut self.residues, r);
                self.class.push(Some(c));
                self.residue.push(Some(r));
            }
        }
        self.amber.push(amber);
    }

    /// The same labels with every protein atom's ff14SB `[x, D]` replaced by `vdw`'s: what makes
    /// variant 3's protein the current one.
    fn amber_as(&self, vdw: &[[f64; 2]]) -> Labels {
        let mut l = self.clone();
        for (i, a) in l.amber.iter_mut().enumerate() {
            if a.is_some() {
                *a = Some(vdw[i]);
            }
        }
        l
    }
}

/// An ff14SB type's `[x, D]` in UFF's form, as [`gaff_benzene`] writes GAFF's: `x = 2 R*/2`,
/// `D = ε` ([`amber`]).
fn amber_vdw(t: &str) -> [f64; 2] {
    let (r, e) = amber::lennard_jones(t);
    [2.0 * r * ANGSTROM, e * KCAL_PER_MOL]
}

/// [`Labels`] for a solvated complex taken apart as [`complex_surroundings`] takes it: each protein
/// atom by its ff14SB type from `system`'s residue templates ([`amber::atom_type`]), through the
/// binding's map back to the system.
fn complex_labels(s: &SolvatedComplex, system: &System, sur: &Surroundings) -> Labels {
    let map = s.binding().system_atoms();
    let mut labels = Labels::new();
    for (i, kind) in sur.kinds.iter().enumerate() {
        match kind {
            None => labels.push(None, None),
            Some(Kind::Protein(_)) => {
                let t = amber::atom_type(system, map[i]);
                let residue = system.residue_of(map[i]).label();
                labels.push(Some((t, &residue)), Some(amber_vdw(t)));
            }
            Some(Kind::Water) => labels.push(Some(("water", "water")), None),
            Some(Kind::Ion) => labels.push(Some(("ions", "ions")), None),
        }
    }
    labels
}

/// [`Labels`] for W3's benzene in water: every atom around it water's.
fn water_labels(sur: &Surroundings) -> Labels {
    let mut labels = Labels::new();
    for kind in &sur.kinds {
        match kind {
            None => labels.push(None, None),
            Some(Kind::Water) => labels.push(Some(("water", "water")), None),
            Some(other) => panic!("W3's box has a {other:?}"),
        }
    }
    labels
}

/// The model the field is: every `[x, D]` as [`Surroundings`] rebuilt them, combined geometrically.
fn current_model(sur: &Surroundings) -> CrossModel {
    CrossModel {
        vdw: sur.vdw.clone(),
        mixing: Mixing::Geometric,
    }
}

/// A variant: [`current_model`] with **benzene's `[x, D]` replaced by `benzene`** when given,
/// **every protein atom's by its ff14SB `[x, D]`** from `labels` when `amber`, and combined by
/// `mixing`. Water stays TIP3P's (hydrogens without a term) and the ions Joung and Cheatham's.
fn cross_model(
    sur: &Surroundings,
    labels: &Labels,
    benzene: Option<&[[f64; 2]]>,
    amber: bool,
    mixing: Mixing,
) -> CrossModel {
    let mut vdw = sur.vdw.clone();
    if let Some(b) = benzene {
        assert_eq!(b.len(), sur.ligand.len());
        for (k, i) in sur.ligand.clone().enumerate() {
            vdw[i] = b[k];
        }
    }
    if amber {
        for (i, kind) in sur.kinds.iter().enumerate() {
            if let Some(Kind::Protein(_)) = kind {
                vdw[i] = labels.amber[i].expect("a protein atom has an ff14SB type");
            }
        }
    }
    CrossModel { vdw, mixing }
}

/// The four variants of [`VARIANTS`], GAFF's benzene [`gaff_benzene`]'s.
fn variants(sur: &Surroundings, labels: &Labels) -> Vec<CrossModel> {
    let gaff = gaff_benzene(&ring_and_hydrogens(sur));
    vec![
        cross_model(sur, labels, Some(&gaff), false, Mixing::Geometric),
        cross_model(sur, labels, None, false, Mixing::LorentzBerthelot),
        cross_model(sur, labels, Some(&gaff), false, Mixing::LorentzBerthelot),
        cross_model(sur, labels, Some(&gaff), true, Mixing::LorentzBerthelot),
    ]
}

/// A sum of benzene's cross van der Waals: the pairs inside the cutoff, the long-range correction,
/// and `Σ|term|` and the number of terms for its rounding.
#[derive(Clone, Copy, Debug, Default)]
struct CrossSum {
    pairs: f64,
    tail: f64,
    magnitude: f64,
    terms: usize,
}

impl CrossSum {
    fn total(&self) -> f64 {
        self.pairs + self.tail
    }
}

/// **The long-range correction** of benzene's van der Waals with the atoms `select` keeps, under
/// `m`: `(4π/V) Σ_g Σ_j D_gj [x_gj¹²/(9 r_c⁹) − 2 x_gj⁶/(3 r_c³)]`, the integral of the pair beyond
/// the cutoff over a uniform rest, with the rest's atoms counted by their `[x, D]` so that a box of
/// waters is one term per benzene atom. Under the geometric rule it is [`Surroundings::van_der_waals`]'s
/// `S₆ S₆ᴮ`, `S₃ S₃ᴮ` form, which Lorentz–Berthelot does not factor into. It does not move with the
/// atoms: computed once per model.
fn cross_tail(sur: &Surroundings, m: &CrossModel, select: impl Fn(usize) -> bool) -> CrossSum {
    let mut groups: std::collections::BTreeMap<[u64; 2], f64> = Default::default();
    for (j, kind) in sur.kinds.iter().enumerate() {
        if kind.is_some() && select(j) {
            *groups.entry(m.vdw[j].map(f64::to_bits)).or_default() += 1.0;
        }
    }
    let rc3 = sur.cutoff * sur.cutoff * sur.cutoff;
    let rc9 = rc3 * rc3 * rc3;
    let scale = 4.0 * std::f64::consts::PI / sur.cell.volume();
    let mut s = CrossSum::default();
    for g in sur.ligand.clone() {
        for (bits, count) in &groups {
            let p = m.mixing.pair(m.vdw[g], bits.map(f64::from_bits));
            let x3 = p[0] * p[0] * p[0];
            let x6 = x3 * x3;
            let t = count * scale * p[1] * (x6 * x6 / (9.0 * rc9) - 2.0 * x6 / (3.0 * rc3));
            s.tail += t;
            s.magnitude += t.abs();
            s.terms += 1;
        }
    }
    s
}

/// **Benzene's whole cross van der Waals under `m` at `at`**: every pair across inside the cutoff,
/// by its own loop, and `tail`, [`cross_tail`] over every atom.
fn cross_sum(sur: &Surroundings, at: &[[f64; 3]], m: &CrossModel, tail: CrossSum) -> CrossSum {
    let w: Vec<[f64; 3]> = at.iter().map(|p| sur.cell.wrap(*p)).collect();
    let rc2 = sur.cutoff * sur.cutoff;
    let mut s = tail;
    for g in sur.ligand.clone() {
        for (j, kind) in sur.kinds.iter().enumerate() {
            if kind.is_none() {
                continue;
            }
            let (a, b) = (g.min(j), g.max(j));
            let (_, r2) = sur.separation(&w, a, b);
            if r2 >= rc2 {
                continue;
            }
            let p = m.mixing.pair(m.vdw[a], m.vdw[b]);
            let u = Pair {
                atoms: [a, b],
                distance: p[0],
                well: p[1],
            }
            .energy(r2.sqrt());
            s.pairs += u;
            s.magnitude += u.abs();
            s.terms += 1;
        }
    }
    s
}

/// `ΔU` from one model to another at one frame, **pair by pair**: each pair's `u_to − u_from` added
/// to the atom of the rest it is with, and each class's long-range correction's difference.
#[derive(Clone, Debug)]
struct Delta {
    per_atom: Vec<f64>,
    tail: Vec<f64>,
    magnitude: f64,
    terms: usize,
}

impl Delta {
    fn total(&self) -> f64 {
        self.per_atom.iter().sum::<f64>() + self.tail.iter().sum::<f64>()
    }

    /// By class, pairs and long-range correction.
    fn by_class(&self, labels: &Labels) -> Vec<f64> {
        let mut c = self.tail.clone();
        for (j, u) in self.per_atom.iter().enumerate() {
            if let Some(k) = labels.class[j] {
                c[k] += u;
            }
        }
        c
    }

    /// By residue, pairs only: the long-range correction belongs to no residue.
    fn by_residue(&self, labels: &Labels) -> Vec<f64> {
        let mut r = vec![0.0; labels.residues.len()];
        for (j, u) in self.per_atom.iter().enumerate() {
            if let Some(k) = labels.residue[j] {
                r[k] += u;
            }
        }
        r
    }
}

/// Each class's [`cross_tail`] under `m`.
fn class_tails(sur: &Surroundings, labels: &Labels, m: &CrossModel) -> Vec<CrossSum> {
    (0..labels.classes.len())
        .map(|c| cross_tail(sur, m, |j| labels.class[j] == Some(c)))
        .collect()
}

/// **`ΔU` from `base` to each of `to` at `at`, pair by pair** ([`Delta`]): one search for the
/// pairs inside the cutoff, `base`'s energy once per pair and each model's beside it; the
/// long-range corrections from each model's [`class_tails`].
fn cross_deltas(
    sur: &Surroundings,
    labels: &Labels,
    at: &[[f64; 3]],
    base: (&CrossModel, &[CrossSum]),
    to: &[(&CrossModel, &[CrossSum])],
) -> Vec<Delta> {
    let n = sur.kinds.len();
    let mut deltas: Vec<Delta> = to
        .iter()
        .map(|(_, tails)| {
            let mut d = Delta {
                per_atom: vec![0.0; n],
                tail: vec![0.0; labels.classes.len()],
                magnitude: 0.0,
                terms: 0,
            };
            for c in 0..labels.classes.len() {
                d.tail[c] = tails[c].tail - base.1[c].tail;
                d.magnitude += tails[c].magnitude + base.1[c].magnitude;
                d.terms += tails[c].terms;
            }
            d
        })
        .collect();
    let w: Vec<[f64; 3]> = at.iter().map(|p| sur.cell.wrap(*p)).collect();
    let rc2 = sur.cutoff * sur.cutoff;
    let energy = |m: &CrossModel, a: usize, b: usize, r: f64| {
        let p = m.mixing.pair(m.vdw[a], m.vdw[b]);
        Pair {
            atoms: [a, b],
            distance: p[0],
            well: p[1],
        }
        .energy(r)
    };
    for g in sur.ligand.clone() {
        for (j, kind) in sur.kinds.iter().enumerate() {
            if kind.is_none() {
                continue;
            }
            let (a, b) = (g.min(j), g.max(j));
            let (_, r2) = sur.separation(&w, a, b);
            if r2 >= rc2 {
                continue;
            }
            let r = r2.sqrt();
            let u0 = energy(base.0, a, b, r);
            for (d, (m, _)) in deltas.iter_mut().zip(to) {
                let u = energy(m, a, b, r);
                d.per_atom[j] += u - u0;
                d.magnitude += u.abs() + u0.abs();
                d.terms += 1;
            }
        }
    }
    deltas
}

/// The summation bound for comparing two sums of `terms` terms whose magnitudes sum to
/// `magnitude`: `2 (n + 1) ε Σ|term|`, [`Surroundings::identity`]'s.
fn summation_bound(terms: usize, magnitude: f64) -> f64 {
    2.0 * (terms as f64 + 1.0) * EPS * magnitude
}

/// The variants at one frame: each one's `ΔU` pair by pair, held to the difference of the two
/// whole sums ([`cross_sum`]), and variant 3's `ΔU` by class and by residue.
struct FrameDeltas {
    du: Vec<f64>,
    u: f64,
    u_variants: Vec<f64>,
    classes: Vec<f64>,
    residues: Vec<f64>,
    /// `|ΔU − (U_v − U_current)|` over its bound, the worst variant's.
    against_sums: f64,
    /// Variant 3's classes' sum, and its residues' and corrections', against its `ΔU`, over the
    /// bound.
    parts: f64,
}

/// Every model's long-range corrections, which do not move with the atoms: the whole one and the
/// classes', the current model's first.
struct Tails {
    whole: Vec<CrossSum>,
    classes: Vec<Vec<CrossSum>>,
}

fn tails(sur: &Surroundings, labels: &Labels, models: &[&CrossModel]) -> Tails {
    Tails {
        whole: models
            .iter()
            .map(|m| cross_tail(sur, m, |_| true))
            .collect(),
        classes: models.iter().map(|m| class_tails(sur, labels, m)).collect(),
    }
}

/// [`FrameDeltas`] at `at`, `models` the current model and then the variants, `tails` theirs.
fn frame_deltas(
    sur: &Surroundings,
    labels: &Labels,
    at: &[[f64; 3]],
    models: &[&CrossModel],
    tails: &Tails,
) -> FrameDeltas {
    let to: Vec<(&CrossModel, &[CrossSum])> = (1..models.len())
        .map(|k| (models[k], &tails.classes[k][..]))
        .collect();
    let deltas = cross_deltas(sur, labels, at, (models[0], &tails.classes[0]), &to);
    let sums: Vec<CrossSum> = models
        .iter()
        .zip(&tails.whole)
        .map(|(m, t)| cross_sum(sur, at, m, *t))
        .collect();
    let mut against_sums: f64 = 0.0;
    for (k, d) in deltas.iter().enumerate() {
        let (s, s0) = (&sums[k + 1], &sums[0]);
        let bound = summation_bound(
            d.terms.max(s.terms + s0.terms),
            d.magnitude + s.magnitude + s0.magnitude,
        );
        against_sums = against_sums.max((d.total() - (s.total() - s0.total())).abs() / bound);
    }
    let last = deltas.last().expect("a variant");
    let classes = last.by_class(labels);
    let residues = last.by_residue(labels);
    let bound = summation_bound(last.terms, last.magnitude);
    let by_class: f64 = classes.iter().sum();
    let by_residue: f64 = residues.iter().sum::<f64>() + last.tail.iter().sum::<f64>();
    let parts = (by_class - last.total())
        .abs()
        .max((by_residue - last.total()).abs())
        / bound;
    FrameDeltas {
        du: deltas.iter().map(Delta::total).collect(),
        u: sums[0].total(),
        u_variants: sums[1..].iter().map(CrossSum::total).collect(),
        classes,
        residues,
        against_sums,
        parts,
    }
}

/// The frames' file's first line: `header`, then how many frames of how many atoms.
fn frames_line(header: &str, frames: usize, atoms: usize) -> String {
    assert!(!header.contains('\n'), "a header is one line");
    format!("{header} frames {frames} atoms {atoms}")
}

/// FNV-1a, 64-bit, [`digest`]'s, carried over `bytes` from `h`.
fn fnv(mut h: u64, bytes: &[u8]) -> u64 {
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// FNV-1a's offset basis, [`digest`]'s start.
const FNV_BASIS: u64 = 0xcbf2_9ce4_8422_2325;

/// **Frames written as they are sampled**, in binary: [`frames_line`] and a newline, then each
/// frame as its step (`u64`) and every position's three coordinates (`f64` bits), little-endian,
/// then a trailer `# end N digest` and a newline, the digest FNV-1a over every byte between the
/// header and the trailer. **Written to the path with the extension `partial`, begun afresh, and
/// renamed once complete**: a phase killed while sampling leaves no frames file.
struct FrameWriter {
    file: std::io::BufWriter<std::fs::File>,
    partial: std::path::PathBuf,
    path: std::path::PathBuf,
    frames: usize,
    atoms: usize,
    written: usize,
    hash: u64,
}

impl FrameWriter {
    fn create(path: &std::path::Path, header: &str, frames: usize, atoms: usize) -> FrameWriter {
        use std::io::Write;
        let partial = path.with_extension("partial");
        let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("frames"));
        writeln!(file, "{}", frames_line(header, frames, atoms)).expect("frames");
        FrameWriter {
            file,
            partial,
            path: path.to_path_buf(),
            frames,
            atoms,
            written: 0,
            hash: FNV_BASIS,
        }
    }

    fn push(&mut self, step: u64, at: &[[f64; 3]]) {
        use std::io::Write;
        assert_eq!(at.len(), self.atoms, "a frame of the atoms the header says");
        assert!(
            self.written < self.frames,
            "more frames than the header says"
        );
        let mut bytes = Vec::with_capacity(8 + 24 * at.len());
        bytes.extend_from_slice(&step.to_le_bytes());
        for p in at {
            for x in p {
                bytes.extend_from_slice(&x.to_bits().to_le_bytes());
            }
        }
        self.hash = fnv(self.hash, &bytes);
        self.file.write_all(&bytes).expect("frames");
        self.written += 1;
    }

    fn finish(self) {
        use std::io::Write;
        let FrameWriter {
            mut file,
            partial,
            path,
            frames,
            written,
            hash,
            ..
        } = self;
        assert_eq!(written, frames, "fewer frames than the header says");
        writeln!(file, "# end {frames} {hash:016x}").expect("frames");
        file.flush().expect("frames");
        drop(file);
        std::fs::rename(&partial, &path).expect("rename");
    }
}

/// Frames read back: the first line, and each frame's step and positions.
struct Frames {
    line: String,
    steps: Vec<u64>,
    at: Vec<Vec<[f64; 3]>>,
}

/// What [`FrameWriter`] wrote at `path`: `None` when there is no file, the frames when it is whole,
/// and why not when it is not.
fn read_frames(path: &std::path::Path) -> Option<Result<Frames, String>> {
    let file = std::fs::File::open(path).ok()?;
    Some(parse_frames(&mut std::io::BufReader::new(file)))
}

/// [`read_frames`]'s reading: refused for a first line that does not end in `frames N atoms M`, a
/// frame cut short, or anything after the last frame that is not exactly the trailer with the
/// count and the digest of the bytes read.
fn parse_frames(r: &mut impl std::io::BufRead) -> Result<Frames, String> {
    let mut first = Vec::new();
    r.read_until(b'\n', &mut first).map_err(|e| e.to_string())?;
    if first.pop() != Some(b'\n') {
        return Err("no whole first line".into());
    }
    let line = String::from_utf8(first).map_err(|_| "a first line that is not text".to_string())?;
    let words: Vec<&str> = line.split_whitespace().collect();
    let k = words.len();
    if k < 4 || words[k - 4] != "frames" || words[k - 2] != "atoms" {
        return Err(format!("a first line without its counts: {line:?}"));
    }
    let frames: usize = words[k - 3]
        .parse()
        .map_err(|_| "a frame count".to_string())?;
    let atoms: usize = words[k - 1]
        .parse()
        .map_err(|_| "an atom count".to_string())?;
    let mut buf = vec![0u8; 8 + 24 * atoms];
    let (mut steps, mut at) = (Vec::with_capacity(frames), Vec::with_capacity(frames));
    let mut hash = FNV_BASIS;
    let word = |b: &[u8]| u64::from_le_bytes(b.try_into().expect("eight bytes"));
    for f in 0..frames {
        r.read_exact(&mut buf)
            .map_err(|_| format!("cut short in frame {f} of {frames}"))?;
        hash = fnv(hash, &buf);
        steps.push(word(&buf[..8]));
        at.push(
            (0..atoms)
                .map(|a| [0, 1, 2].map(|x| f64::from_bits(word(&buf[8 + 24 * a + 8 * x..][..8]))))
                .collect::<Vec<[f64; 3]>>(),
        );
    }
    let mut rest = Vec::new();
    r.read_to_end(&mut rest).map_err(|e| e.to_string())?;
    let want = format!("# end {frames} {hash:016x}\n");
    if rest != want.as_bytes() {
        return Err(format!(
            "the trailer is not {want:?}: {:?}",
            String::from_utf8_lossy(&rest)
        ));
    }
    Ok(Frames { line, steps, at })
}

/// [`digest`] of a Boresch restraint: its anchors, references and force constants, as [`header`]
/// digests them.
fn restraint_digest(r: &Boresch) -> u64 {
    digest(
        r.receptor
            .iter()
            .chain(&r.ligand)
            .map(|&i| i as f64)
            .chain(r.reference)
            .chain(r.force_constants),
    )
}

/// One phase of the apportioning: a box taken apart ([`Surroundings`], [`Labels`]), the dynamics
/// its frames are sampled with, and what its frames are run from — `start`, called only when the
/// frames have to be run, and `identity`, which says what `start` is and goes into the header.
struct Phase<'a> {
    name: &'static str,
    sur: &'a Surroundings,
    labels: &'a Labels,
    decoupling: &'a PeriodicDecoupling,
    template: &'a MolecularDynamics,
    lambda: Lambda,
    identity: String,
    start: &'a dyn Fn() -> Vec<[f64; 3]>,
}

/// A phase's frames' header: its name, the protocol literally, the state, and its identity.
fn phase_header(phase: &Phase, p: Protocol) -> String {
    let l = phase.lambda;
    format!(
        "# frames {} protocol {} {} {} {} {} {} {:#x} index 0 lambda {} {} {} {}",
        phase.name,
        exact(p.time_step),
        exact(p.temperature),
        exact(p.friction),
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        exact(l.restraint),
        exact(l.electrostatics),
        exact(l.van_der_waals),
        phase.identity
    )
}

/// **A phase's frames: read back from `NAME.frames` in `dir` when a run completed them under the
/// same header**, refused by name when that file is not whole or was run differently, or run now
/// through [`sample_frames`] at `p` — window 0's seed, 5794556's — streamed to `NAME.partial` and
/// renamed, then read back from the file, so that what is measured is what was saved.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn phase_frames(dir: &std::path::Path, phase: &Phase, p: Protocol) -> Frames {
    let path = dir.join(format!("{}.frames", phase.name));
    let header = phase_header(phase, p);
    let atoms = phase.sur.kinds.len();
    let want = frames_line(&header, p.samples, atoms);
    if let Some(read) = read_frames(&path) {
        let f = read.unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it to run the phase again",
                path.display()
            )
        });
        assert_eq!(
            f.line,
            want,
            "{} was run differently: remove it to run the phase again",
            path.display()
        );
        assert!(
            f.steps.windows(2).all(|s| s[0] < s[1]),
            "{}: steps that do not rise",
            path.display()
        );
        log(
            dir,
            "variants",
            &format!("{}: {} frames read back", phase.name, f.at.len()),
        );
        return f;
    }
    let start = (phase.start)();
    assert_eq!(start.len(), atoms);
    let t = std::time::Instant::now();
    log(
        dir,
        "variants",
        &format!(
            "{}: {} steps begin, {} frames",
            phase.name,
            p.steps_per_window(),
            p.samples
        ),
    );
    let mut writer = FrameWriter::create(&path, &header, p.samples, atoms);
    let mut last = (0u64, 0.0f64);
    sample_frames(
        phase.decoupling,
        &start,
        phase.template,
        p,
        0,
        phase.lambda,
        |step, at| {
            writer.push(step, at);
            if writer.written % 40 == 0 {
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    "variants",
                    &format!(
                        "{}: frame {} at step {step}, {secs:.0} s, {:.1} ms a step since the last \
                         report",
                        phase.name,
                        writer.written,
                        (secs - last.1) * 1e3 / (step - last.0) as f64
                    ),
                );
                last = (step, secs);
            }
        },
    );
    writer.finish();
    log(
        dir,
        "variants",
        &format!(
            "{}: frames written in {:.0} s",
            phase.name,
            t.elapsed().as_secs_f64()
        ),
    );
    read_frames(&path).expect("just written").expect("whole")
}

/// A phase's frames, reparametrised: per frame, each variant's `ΔU`, the current model's whole
/// cross van der Waals and each variant's, variant 3's by class and by residue; and the worst of
/// each identity over the frames, over its bound.
#[derive(Clone, Debug, PartialEq)]
struct Reparametrised {
    du: Vec<Vec<f64>>,
    u: Vec<f64>,
    u_variants: Vec<Vec<f64>>,
    classes: Vec<Vec<f64>>,
    residues: Vec<Vec<f64>>,
    /// Each variant's `ΔU` against the difference of the two whole sums.
    against_sums: f64,
    /// Variant 3's classes, and its residues and corrections, against its `ΔU`.
    parts: f64,
    /// The current model's whole sum against the field's own coupling, every `every`th frame.
    against_field: f64,
    /// Variant 0's whole sum against [`Surroundings::van_der_waals`] with GAFF's benzene, 5794556's
    /// path, every frame.
    against_hybrid: f64,
}

/// **Every frame of a phase reparametrised** ([`frame_deltas`]), the current model's whole sum held
/// to the field's coupling every `every`th frame, and variant 0's to 5794556's own sum on every one.
fn reparametrise(phase: &Phase, frames: &Frames, every: usize) -> Reparametrised {
    let (sur, labels) = (phase.sur, phase.labels);
    let current = current_model(sur);
    let vs = variants(sur, labels);
    let mut models = vec![&current];
    models.extend(vs.iter());
    let tails = tails(sur, labels, &models);
    let gaff = gaff_benzene(&ring_and_hydrogens(sur));
    let mut r = Reparametrised {
        du: Vec::new(),
        u: Vec::new(),
        u_variants: Vec::new(),
        classes: Vec::new(),
        residues: Vec::new(),
        against_sums: 0.0,
        parts: 0.0,
        against_field: 0.0,
        against_hybrid: 0.0,
    };
    for (k, at) in frames.at.iter().enumerate() {
        let f = frame_deltas(sur, labels, at, &models, &tails);
        let whole = cross_sum(sur, at, &current, tails.whole[0]);
        if k % every == 0 {
            let (field, _) = sur.coupling(at);
            let bound = summation_bound(whole.terms, whole.magnitude);
            r.against_field = r.against_field.max((whole.total() - field).abs() / bound);
        }
        let hybrid = sur.van_der_waals(at, &gaff);
        let v0 = cross_sum(sur, at, &vs[0], tails.whole[1]);
        let bound = summation_bound(hybrid.terms.max(v0.terms), hybrid.magnitude + v0.magnitude);
        r.against_hybrid = r
            .against_hybrid
            .max((hybrid.total() - v0.total()).abs() / bound);
        r.against_sums = r.against_sums.max(f.against_sums);
        r.parts = r.parts.max(f.parts);
        r.du.push(f.du);
        r.u.push(f.u);
        r.u_variants.push(f.u_variants);
        r.classes.push(f.classes);
        r.residues.push(f.residues);
    }
    r
}

/// **The apportioning**: each phase's frames ([`phase_frames`], read back or run and saved), then
/// every frame reparametrised ([`reparametrise`]). The frames are the expensive part and each
/// phase's are saved whole once run, so a run killed in the water phase reruns only that.
fn apportion(
    dir: &std::path::Path,
    p: Protocol,
    phases: [&Phase; 2],
    every: usize,
) -> [Reparametrised; 2] {
    phases.map(|phase| {
        let frames = phase_frames(dir, phase, p);
        let t = std::time::Instant::now();
        let r = reparametrise(phase, &frames, every);
        log(
            dir,
            "variants",
            &format!(
                "{}: {} frames reparametrised in {:.0} s",
                phase.name,
                frames.at.len(),
                t.elapsed().as_secs_f64()
            ),
        );
        r
    })
}

/// `ΔG°_bind` with W4's fifteen windows (fd2881e), kcal/mol, and experiment's.
const W4_BINDING: (f64, f64) = (-8.975, 0.279);
const EXPERIMENT: (f64, f64) = (-5.19, 0.16);

/// σ(βΔU) above which one step of perturbation is not trusted here.
const ONE_STEP_SPREAD: f64 = 1.5;

/// The tables: each variant's `ΔU` bound and in water, Zwanzig's estimate, the second cumulant
/// and the spread, the effective samples and the halves, and `ΔΔG_bind`; then variant 3's by class
/// and by residue; and the pairs.
fn report_apportioned(kt: f64, bound: &Reparametrised, water: &Reparametrised, labels: &Labels) {
    let n = bound.du.len();
    let mut summary = Vec::new();
    for (v, name) in VARIANTS.iter().enumerate() {
        let du_b: Vec<f64> = bound.du.iter().map(|d| d[v]).collect();
        let du_w: Vec<f64> = water.du.iter().map(|d| d[v]).collect();
        let (zb, zw) = (zwanzig(&du_b, kt), zwanzig(&du_w, kt));
        let ub = estimate_kcal(bound.u.iter().copied());
        let uw = estimate_kcal(water.u.iter().copied());
        let vb = estimate_kcal(bound.u_variants.iter().map(|u| u[v]));
        let vw = estimate_kcal(water.u_variants.iter().map(|u| u[v]));
        println!("\n**Variant {name}**\n\n| | bound | in water |\n| --- | --- | --- |");
        println!("| ⟨U⟩, current | {} | {} |", pm(&ub), pm(&uw));
        println!("| ⟨U⟩, the variant's | {} | {} |", pm(&vb), pm(&vw));
        for (row, f) in [
            (
                "⟨ΔU⟩",
                &(|z: &Zwanzig, _: f64| pm(&z.mean)) as &dyn Fn(&Zwanzig, f64) -> String,
            ),
            ("⟨ΔU⟩ over ⟨U⟩", &|z: &Zwanzig, u: f64| {
                format!("{:.1}%", 100.0 * z.mean.mean / u.abs())
            }),
            ("σ(βΔU)", &|z: &Zwanzig, _: f64| {
                format!("{:.3}", z.spread)
            }),
            ("⟨ΔU⟩ − βσ²/2", &|z: &Zwanzig, _: f64| {
                format!("{:+.3}", z.cumulant)
            }),
            (
                "Zwanzig −kT ln⟨e^(−βΔU)⟩",
                &|z: &Zwanzig, _: f64| format!("{:+.3} ± {:.3}", z.exp, z.exp_error),
            ),
            ("Zwanzig, halves", &|z: &Zwanzig, _: f64| {
                format!("{:+.3}, {:+.3}", z.halves[0], z.halves[1])
            }),
            ("effective samples", &|z: &Zwanzig, _: f64| {
                format!("{:.1} of {n}", z.ess)
            }),
        ] {
            println!("| {row} | {} | {} |", f(&zb, ub.mean), f(&zw, uw.mean));
        }
        let first = zb.mean.mean - zw.mean.mean;
        let first_error = (zb.mean.error * zb.mean.error + zw.mean.error * zw.mean.error).sqrt();
        let exp = zb.exp - zw.exp;
        let exp_error = (zb.exp_error * zb.exp_error + zw.exp_error * zw.exp_error).sqrt();
        let halves = [zb.halves[0] - zw.halves[0], zb.halves[1] - zw.halves[1]];
        let reliable = zb.spread <= ONE_STEP_SPREAD && zw.spread <= ONE_STEP_SPREAD;
        println!(
            "\nΔΔG_bind = ΔG_swap(bound) − ΔG_swap(water): first order {first:+.3} ± {first_error:.3}; \
             second cumulant {:+.3}; Zwanzig {exp:+.3} ± {exp_error:.3} (halves {:+.3}, {:+.3}) \
             kcal/mol; ΔG°_bind {:+.3} + ({exp:+.3}) = {:+.3} against experiment's {:+.2}.{}",
            zb.cumulant - zw.cumulant,
            halves[0],
            halves[1],
            W4_BINDING.0,
            W4_BINDING.0 + exp,
            EXPERIMENT.0,
            if reliable {
                String::new()
            } else {
                format!(
                    " **σ(βΔU) above {ONE_STEP_SPREAD}: one step of perturbation is not to be \
                     trusted here; it wants staged perturbation through intermediate parameters, \
                     or a run of the variant itself.**"
                )
            }
        );
        summary.push((
            name,
            first,
            zb.cumulant - zw.cumulant,
            exp,
            exp_error,
            reliable,
        ));
    }

    println!(
        "\n| variant | ΔΔG_bind, first order | second cumulant | Zwanzig | ΔG°_bind | σ(βΔU) ≤ \
         {ONE_STEP_SPREAD} |\n| --- | --- | --- | --- | --- | --- |"
    );
    for (name, first, cumulant, exp, error, reliable) in &summary {
        println!(
            "| {name} | {first:+.3} | {cumulant:+.3} | {exp:+.3} ± {error:.3} | {:+.3} | {} |",
            W4_BINDING.0 + exp,
            if *reliable { "yes" } else { "**no**" }
        );
    }
    println!(
        "\nsteps, as differences of the Zwanzig estimates: Lorentz–Berthelot alone {:+.3}; then \
         GAFF's benzene {:+.3}; then ff14SB's protein {:+.3}. GAFF's benzene under the geometric \
         rule, 5794556's +1.711, is {:+.3} on these frames. The binding is {:+.2} too strong.",
        summary[1].3,
        summary[2].3 - summary[1].3,
        summary[3].3 - summary[2].3,
        summary[0].3,
        W4_BINDING.0 - EXPERIMENT.0
    );

    // Variant 3 by class and by residue, to first order.
    let du_w3 = estimate_kcal(water.du.iter().map(|d| d[3]));
    let mut rows: Vec<(String, Estimate)> = labels
        .classes
        .iter()
        .enumerate()
        .map(|(c, name)| {
            (
                name.clone(),
                estimate_kcal(bound.classes.iter().map(|x| x[c])),
            )
        })
        .collect();
    rows.sort_by(|a, b| b.1.mean.abs().total_cmp(&a.1.mean.abs()));
    println!(
        "\n**Variant 3 by class**, first order: ΔΔG_bind ≈ Σ_c ⟨ΔU_c⟩_bound − ⟨ΔU⟩_water, the water \
         class's in-water ⟨ΔU⟩ {} beside it.\n\n| class | ⟨ΔU_c⟩ bound (kcal/mol) |\n| --- | --- |",
        pm(&du_w3)
    );
    let mut sum = 0.0;
    for (name, e) in &rows {
        sum += e.mean;
        println!("| {name} | {} |", pm(e));
    }
    println!(
        "| **all** | {sum:+.3} |\n\nfirst order {sum:+.3} − ({:+.3}) = {:+.3}",
        du_w3.mean,
        sum - du_w3.mean
    );
    let mut res: Vec<(String, Estimate)> = labels
        .residues
        .iter()
        .enumerate()
        .map(|(r, name)| {
            (
                name.clone(),
                estimate_kcal(bound.residues.iter().map(|x| x[r])),
            )
        })
        .collect();
    res.sort_by(|a, b| b.1.mean.abs().total_cmp(&a.1.mean.abs()));
    println!(
        "\n**Variant 3 by residue**, pairs only:\n\n| residue | ⟨ΔU_r⟩ bound |\n| --- | --- |"
    );
    for (name, e) in res.iter().take(15) {
        println!("| {name} | {} |", pm(e));
    }
    let rest: f64 = res.iter().skip(15).map(|(_, e)| e.mean).sum();
    let corrections = sum - res.iter().map(|(_, e)| e.mean).sum::<f64>();
    println!(
        "| {} more | {rest:+.3} |\n| long-range corrections | {corrections:+.3} |",
        res.len().saturating_sub(15)
    );
}

/// **The rest of W4's 3.79 kcal/mol, apportioned** (release, one core; the frames are the cost:
/// 10 000 steps of the complex, about fifty minutes at 300 ms a step and two and a half hours at
/// 5794556's throttled 710, and a few minutes of W3's box): benzene's cross van der Waals
/// reparametrised by Zwanzig's one-step perturbation on frames sampled at UFF's, bound and in
/// water, for each of [`VARIANTS`]:
///
/// 1. **Lorentz–Berthelot alone**: benzene's UFF `[x, D]` mixed with everything by
///    `x_ij = (x_i + x_j)/2`, `D_ij = √(D_i D_j)` ([`Mixing`]);
/// 2. **GAFF's benzene with Lorentz–Berthelot**, the environment UFF's and TIP3P's;
/// 3. **an AMBER-like protein**: every protein atom's Lennard-Jones ff14SB's by its type from its
///    residue's template ([`amber`]), GAFF's benzene, Lorentz–Berthelot, TIP3P's water and Joung and
///    Cheatham's chloride as they are;
///
/// and 0, 5794556's hybrid again, whose +1.711 these frames should give. Each with `⟨ΔU⟩`,
/// Zwanzig's estimate, the second cumulant, σ(βΔU), the effective samples and the halves, and
/// `ΔΔG_bind = ΔG_swap(bound) − ΔG_swap(water)`; a variant whose σ(βΔU) is above
/// [`ONE_STEP_SPREAD`] on either side is flagged. Variant 3 is then taken apart by ff14SB type and
/// by residue, to first order.
///
/// - **The frames are saved** in `frames/` under `PANTOMETRY_W4_DIR` — `complex.frames` and
///   `water.frames`, binary, every bit ([`FrameWriter`]), about 230 MB and 13 MB — so that a later
///   variant needs no new dynamics. Each phase's is read back when a run completed it under the
///   same header and refused by name when it is not whole or was run differently; a phase killed
///   while sampling is begun afresh. The lock is `frames/variants.lock`, the log
///   `frames/variants.log`. **Nothing is written beside the leg's files**: `prepared.txt` is read
///   back through its header check and refused if it is not there.
/// - **The protocol is 5794556's**: the leg's box, start and Hamiltonian at window 0's state,
///   [`FRAMES`] — 4 ps discarded, 320 frames 50 fs apart — through [`sample_frames`] with window 0's
///   seed; then W3's box melted as W3 melts it, the same way.
/// - **Every frame held**: each variant's `ΔU` pair by pair against the difference of its two
///   whole sums, variant 3's classes and residues against its `ΔU`, variant 0's whole sum against
///   5794556's own ([`Surroundings::van_der_waals`]), and every twentieth frame the current model's
///   against the field's coupling — each to its summation bound, asserted after everything is
///   printed.
///
/// **Run once**, at `2e09cea`, exit 0 on 2026-10-17 after 3 669 s. Zwanzig's ΔΔG_bind: 0
/// +1.711 ± 0.056 (5794556's, reproduced), 1 −0.059 ± 0.014, 2 +1.698 ± 0.047 and 3
/// +3.365 ± 0.073, which puts `ΔG°_bind` at −5.61 against experiment's −5.19. σ(βΔU) is at most
/// 0.84. Nothing here asserts it; the CHANGELOG has the tables.
///
/// `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_rest_of_the_binding_error_apportioned_measured -- --ignored
/// --nocapture`.
#[test]
#[ignore = "frames of the whole complex and of W3's box, saved; about an hour with --release the first time, run once (2026-10-17)"]
fn the_rest_of_the_binding_error_apportioned_measured() {
    let dir = results_dir();
    assert!(
        dir.join("prepared.txt").exists(),
        "{} is not there: this measurement reads the leg's preparation and writes nothing beside \
         it",
        dir.join("prepared.txt").display()
    );
    let frames = dir.join("frames");
    std::fs::create_dir_all(&frames).expect("the frames directory");
    let _lock = Lock::take(&frames, "variants");
    let t = std::time::Instant::now();
    let kt = BOLTZMANN.to_si() * KELVIN;
    let s = whole_complex();
    let leg = production_leg(&s, production_start(&s, &dir));
    log(
        &frames,
        "variants",
        &format!(
            "the leg's box and start read back in {:.0} s",
            t.elapsed().as_secs_f64()
        ),
    );
    let g = grid();
    assert_eq!(g[0], Lambda::new(0.0, 1.0, 1.0));
    let bound = complex_surroundings(&s);
    let bound_labels = complex_labels(&s, the_system(), &bound);
    let complex_start = || leg.start.clone();
    let complex = Phase {
        name: "complex",
        sur: &bound,
        labels: &bound_labels,
        decoupling: &leg.decoupling,
        template: &leg.template,
        lambda: g[0],
        identity: format!(
            "start {:016x} restraint {:016x} flexibility {:?}",
            positions_digest(&leg.start),
            restraint_digest(&leg.restraint),
            leg.flexibility
        ),
        start: &complex_start,
    };

    let benzene = hydration::Benzene::new(Some(pocket3().ligand_force_field().charges().to_vec()));
    assert_eq!(
        &benzene.types[..],
        &pocket3().types()[pocket3().ligand_range()]
    );
    let w3 = hydration::w3_solvated(&benzene);
    let w3_decoupling = w3.decoupling();
    let w3_template = w3.dynamics(&benzene);
    let free = water_surroundings(&benzene, &w3);
    let free_labels = water_labels(&free);
    let water_start = || {
        let mut at = w3.at.clone();
        let tm = std::time::Instant::now();
        hydration::w3_melt(&w3_decoupling, &w3_template, KELVIN, &mut at);
        log(
            &frames,
            "variants",
            &format!(
                "W3's box, {} waters, melted in {:.0} s",
                w3.waters,
                tm.elapsed().as_secs_f64()
            ),
        );
        at
    };
    let water = Phase {
        name: "water",
        sur: &free,
        labels: &free_labels,
        decoupling: &w3_decoupling,
        template: &w3_template,
        lambda: Lambda::COUPLED,
        identity: format!(
            "built {:016x} melted w3_melt {}",
            positions_digest(&w3.at),
            exact(KELVIN)
        ),
        start: &water_start,
    };

    let [b, w] = apportion(&frames, FRAMES, [&complex, &water], 20);
    report_apportioned(kt, &b, &w, &bound_labels);

    // The pairs: benzene's carbon and hydrogen with each partner, UFF's and the variants'.
    let uff_c = UffType::CR.parameters();
    let uff_h = UffType::H.parameters();
    let gaff = gaff_benzene(&[Element::C, Element::H]);
    println!(
        "\n| benzene atom with | current x_ij (Å), D_ij (kcal/mol) | 1 | 2 | 3 |\n| --- | --- | --- \
         | --- | --- |"
    );
    let mut seen = std::collections::BTreeSet::new();
    for (j, kind) in bound.kinds.iter().enumerate() {
        let Some(kind) = *kind else { continue };
        let class = bound_labels.class[j].map(|c| bound_labels.classes[c].clone());
        if bound.vdw[j][1] == 0.0 || !seen.insert((kind, class.clone())) {
            continue;
        }
        let env3 = bound_labels.amber[j].unwrap_or(bound.vdw[j]);
        for (name, u, h) in [
            ("C", [uff_c.vdw_distance, uff_c.vdw_energy], gaff[0]),
            ("H", [uff_h.vdw_distance, uff_h.vdw_energy], gaff[1]),
        ] {
            let show = |p: [f64; 2]| format!("{:.4}, {:.4}", p[0] / ANGSTROM, kcal(p[1]));
            println!(
                "| {name} with {} ({}) | {} | {} | {} | {} |",
                kind.label(),
                class.clone().unwrap_or_default(),
                show(Mixing::Geometric.pair(u, bound.vdw[j])),
                show(Mixing::LorentzBerthelot.pair(u, bound.vdw[j])),
                show(Mixing::LorentzBerthelot.pair(h, bound.vdw[j])),
                show(Mixing::LorentzBerthelot.pair(h, env3))
            );
        }
    }

    println!(
        "\nthe identities, worst over the frames (1 is the bound): ΔU against the two sums {:.3} \
         bound, {:.3} in water; variant 3's parts {:.3}, {:.3}; the field's coupling {:.3}, \
         {:.3}; 5794556's sum {:.3}, {:.3}. {:.0} s in all",
        b.against_sums,
        w.against_sums,
        b.parts,
        w.parts,
        b.against_field,
        w.against_field,
        b.against_hybrid,
        w.against_hybrid,
        t.elapsed().as_secs_f64()
    );
    for r in [&b, &w] {
        assert!(r.against_sums <= 1.0, "{}", r.against_sums);
        assert!(r.parts <= 1.0, "{}", r.parts);
        assert!(r.against_field <= 1.0, "{}", r.against_field);
        assert!(r.against_hybrid <= 1.0, "{}", r.against_hybrid);
    }
}

/// The 3 Å pocket in water and W3's small box, each taken apart with its labels.
fn small_phases() -> [(Surroundings, Labels, Vec<[f64; 3]>); 2] {
    let s = small3(0x5EC3);
    let sur = complex_surroundings(&s);
    let labels = complex_labels(&s, the_system(), &sur);
    let at = s.positions().to_vec();
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC4, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    let free = water_surroundings(&b, &w);
    let free_labels = water_labels(&free);
    [(sur, labels, at), (free, free_labels, w.at.clone())]
}

/// **A variant at the current parameters changes nothing, exactly**: on the 3 Å pocket in water and
/// on W3's small box, each variant's recipe ([`cross_model`]) with its parameters set to the current
/// ones — benzene's UFF `[x, D]` where GAFF's would go, the protein's UFF `[x, D]` where ff14SB's
/// would ([`Labels::amber_as`]), the geometric rule — gives `ΔU` = 0 to the bit at every atom and
/// in every class's long-range correction; and every variant against itself rebuilt, Lorentz–Berthelot
/// included, the same.
#[test]
fn a_variant_at_the_current_parameters_changes_nothing() {
    for (sur, labels, at) in small_phases() {
        let current = current_model(&sur);
        let current_tails = class_tails(&sur, &labels, &current);
        let as_current = labels.amber_as(&sur.vdw);
        let uff = sur.uff_benzene();
        let recipes = [
            (Some(&uff[..]), false),
            (None, false),
            (Some(&uff[..]), false),
            (Some(&uff[..]), true),
        ];
        let at_current: Vec<CrossModel> = recipes
            .iter()
            .map(|&(b, amber)| cross_model(&sur, &as_current, b, amber, Mixing::Geometric))
            .collect();
        let vs = variants(&sur, &labels);
        let again = variants(&sur, &labels);
        let ts: Vec<Vec<CrossSum>> = at_current
            .iter()
            .chain(&vs)
            .chain(&again)
            .map(|m| class_tails(&sur, &labels, m))
            .collect();
        let to: Vec<(&CrossModel, &[CrossSum])> = at_current
            .iter()
            .enumerate()
            .map(|(k, m)| (m, &ts[k][..]))
            .collect();
        let mut deltas = cross_deltas(&sur, &labels, &at, (&current, &current_tails), &to);
        for k in 0..vs.len() {
            deltas.extend(cross_deltas(
                &sur,
                &labels,
                &at,
                (&vs[k], &ts[4 + k]),
                &[(&again[k], &ts[8 + k])],
            ));
        }
        for (k, d) in deltas.iter().enumerate() {
            assert!(d.terms > 0, "{k}: pairs");
            for (j, u) in d.per_atom.iter().enumerate() {
                assert_eq!(u.to_bits(), 0, "{k}: atom {j}");
            }
            for (c, u) in d.tail.iter().enumerate() {
                assert_eq!(u.to_bits(), 0, "{k}: class {}", labels.classes[c]);
            }
            assert_eq!(d.total().to_bits(), 0, "{k}");
        }
    }
}

/// **Lorentz–Berthelot is its formula for one pair**, written by hand: GAFF's `ca` with ff14SB's
/// `CT` is `R*_ij = 1.9080 + 1.9080 = 3.8160` Å and `ε_ij = √(0.0860 × 0.1094)` kcal/mol, GAFF's
/// `ha` with ff14SB's `O` 3.1202 Å and `√(0.0150 × 0.2100)`; UFF's `C_R` with TIP3P's oxygen is
/// `x_ij = 2^⅙ (σ_i + σ_j)/2` with `σ = x/2^⅙`; each within a few roundings, either order, and the
/// geometric rule's distance not the same.
#[test]
fn lorentz_berthelot_is_its_formula_for_one_pair() {
    let near = |a: f64, b: f64, what: &str| {
        assert!(
            (a - b).abs() <= 4.0 * EPS * b.abs(),
            "{what}: {a:e} against {b:e}, {:.2} ε",
            (a - b).abs() / (EPS * b.abs())
        );
    };
    let lb = Mixing::LorentzBerthelot;
    let gaff = gaff_benzene(&[Element::C, Element::H]);
    for (benzene, partner, r, e, what) in [
        (
            gaff[0],
            amber_vdw("CT"),
            3.8160,
            (0.0860f64 * 0.1094).sqrt(),
            "ca–CT",
        ),
        (
            gaff[1],
            amber_vdw("O"),
            3.1202,
            (0.0150f64 * 0.2100).sqrt(),
            "ha–O",
        ),
    ] {
        for p in [lb.pair(benzene, partner), lb.pair(partner, benzene)] {
            near(p[0], r * ANGSTROM, what);
            near(p[1], e * KCAL_PER_MOL, what);
        }
    }
    // The two rules agree for ca–CT, whose radii are equal, and not for ha–O.
    let g = Mixing::Geometric.pair(gaff[1], amber_vdw("O"))[0];
    assert!(
        (g - 3.1202 * ANGSTROM).abs() > 1e-3 * ANGSTROM,
        "ha–O: {} Å geometric",
        g / ANGSTROM
    );
    // UFF's aromatic carbon with TIP3P's oxygen, through σ.
    let two_sixth = 1.122_462_048_309_373; // 2^⅙
    let c = UffType::CR.parameters();
    let c = [c.vdw_distance, c.vdw_energy];
    let o = water::oxygen_vdw();
    let sigma = (c[0] / two_sixth + o[0] / two_sixth) / 2.0;
    for p in [lb.pair(c, o), lb.pair(o, c)] {
        near(p[0], two_sixth * sigma, "C_R–O by σ");
        near(p[1], (c[1] * o[1]).sqrt(), "C_R–O depth");
    }
    let g = Mixing::Geometric.pair(c, o)[0];
    assert!(
        (g - two_sixth * sigma).abs() > 1e-3 * ANGSTROM,
        "C_R–O: {} Å geometric against {} Å arithmetic",
        g / ANGSTROM,
        two_sixth * sigma / ANGSTROM
    );
}

/// **The ff14SB types are the templates'**: for a handful of 181L's residues — Leu84, Phe153,
/// Met102, Ala99, Val111, Gly107's backbone, His31 as `HIP`, Met1 at the N terminus, Lys162 at the C
/// terminus, Ser117, Trp126, Arg95 and Asp70 — atom by atom against literal expectations read off
/// `amino12.lib`, `aminont12.lib` and `aminoct12.lib`; the templates those residues take; every
/// protein atom of 181L typed, each type's first letter its element's (`2C` and `3C` carbons), thirty
/// types in all; and the Lennard-Jones lines a few types take, equivalences included.
#[test]
fn the_ff14sb_types_are_the_templates() {
    // A residue's number, its name, the template it takes, and some of its atoms' types.
    type Expected<'a> = (i32, &'a str, &'a str, &'a [(&'a str, &'a str)]);
    let s = the_system();
    let atom = |number: i32, name: &str| -> usize {
        (0..s.component().atoms().len())
            .find(|&i| {
                let r = s.residue_of(i);
                r.part == pantometry_forcefield::Part::Protein
                    && r.number == number
                    && s.atom_name(i) == name
            })
            .unwrap_or_else(|| panic!("no {name} in residue {number}"))
    };
    let expect: [Expected; 13] = [
        (
            84,
            "LEU",
            "LEU",
            &[
                ("N", "N"),
                ("H", "H"),
                ("CA", "CX"),
                ("HA", "H1"),
                ("CB", "2C"),
                ("HB2", "HC"),
                ("CG", "3C"),
                ("HG", "HC"),
                ("CD1", "CT"),
                ("HD11", "HC"),
                ("CD2", "CT"),
                ("C", "C"),
                ("O", "O"),
            ],
        ),
        (
            153,
            "PHE",
            "PHE",
            &[
                ("CB", "CT"),
                ("CG", "CA"),
                ("CD1", "CA"),
                ("CZ", "CA"),
                ("HD1", "HA"),
                ("HZ", "HA"),
            ],
        ),
        (
            102,
            "MET",
            "MET",
            &[
                ("CB", "2C"),
                ("CG", "2C"),
                ("HG2", "H1"),
                ("SD", "S"),
                ("CE", "CT"),
                ("HE1", "H1"),
                ("HB2", "HC"),
            ],
        ),
        (
            99,
            "ALA",
            "ALA",
            &[("CA", "CX"), ("CB", "CT"), ("HB1", "HC")],
        ),
        (
            111,
            "VAL",
            "VAL",
            &[("CB", "3C"), ("HB", "HC"), ("CG1", "CT"), ("HG11", "HC")],
        ),
        (
            107,
            "GLY",
            "GLY",
            &[
                ("N", "N"),
                ("H", "H"),
                ("CA", "CX"),
                ("HA2", "H1"),
                ("HA3", "H1"),
                ("C", "C"),
                ("O", "O"),
            ],
        ),
        (
            31,
            "HIS",
            "HIP",
            &[
                ("CG", "CC"),
                ("ND1", "NA"),
                ("HD1", "H"),
                ("CE1", "CR"),
                ("HE1", "H5"),
                ("NE2", "NA"),
                ("CD2", "CW"),
                ("HD2", "H4"),
            ],
        ),
        (
            1,
            "MET",
            "NMET",
            &[
                ("N", "N3"),
                ("H", "H"),
                ("H2", "H"),
                ("H3", "H"),
                ("CA", "CX"),
                ("HA", "HP"),
            ],
        ),
        (
            162,
            "LYS",
            "CLYS",
            &[
                ("C", "C"),
                ("O", "O2"),
                ("OXT", "O2"),
                ("CE", "C8"),
                ("HE2", "HP"),
                ("NZ", "N3"),
                ("HZ1", "H"),
            ],
        ),
        (
            117,
            "SER",
            "SER",
            &[("OG", "OH"), ("HG", "HO"), ("HB2", "H1")],
        ),
        (
            126,
            "TRP",
            "TRP",
            &[
                ("CG", "C*"),
                ("CD1", "CW"),
                ("CD2", "CB"),
                ("NE1", "NA"),
                ("CE2", "CN"),
                ("HE1", "H"),
                ("HD1", "H4"),
            ],
        ),
        (
            95,
            "ARG",
            "ARG",
            &[("NE", "N2"), ("CZ", "CA"), ("HD2", "H1")],
        ),
        (
            70,
            "ASP",
            "ASP",
            &[("CG", "CO"), ("OD1", "O2"), ("HB2", "HC")],
        ),
    ];
    for (number, name, template, atoms) in expect {
        let first = atom(number, atoms[0].0);
        assert_eq!(s.residue_of(first).name, name, "residue {number}");
        assert_eq!(
            amber::template_name(s.residue_of(first)),
            template,
            "residue {number}"
        );
        for (a, t) in atoms {
            assert_eq!(
                amber::atom_type(s, atom(number, a)),
                *t,
                "{name}{number}:{a}"
            );
        }
    }
    let mut types = std::collections::BTreeSet::new();
    for i in 0..s.component().atoms().len() {
        if s.part(i) != pantometry_forcefield::Part::Protein {
            continue;
        }
        let t = amber::atom_type(s, i);
        amber::lennard_jones(t);
        let element = match t.trim_start_matches(&['2', '3'][..]).as_bytes()[0] {
            b'H' => Element::H,
            b'C' => Element::C,
            b'N' => Element::N,
            b'O' => Element::O,
            b'S' => Element::S,
            other => panic!("type {t}: {}", other as char),
        };
        assert_eq!(s.component().atoms()[i].element, element, "atom {i}: {t}");
        types.insert(t);
    }
    assert_eq!(types.len(), 30, "{types:?}");
    for (t, line, r, e) in [
        ("CT", "CT", 1.9080, 0.1094),
        ("3C", "3C", 1.9080, 0.1094),
        ("CA", "C*", 1.9080, 0.0860),
        ("CO", "CO", 1.9080, 0.0860),
        ("NA", "N", 1.8240, 0.1700),
        ("HO", "HO", 0.0, 0.0),
        ("HP", "HP", 1.1000, 0.0157),
        ("S", "S", 2.0000, 0.2500),
    ] {
        assert_eq!(amber::lennard_jones_line(t), line, "{t}");
        assert_eq!(amber::lennard_jones(t), (r, e), "{t}");
    }
}

/// **Each variant's change is two whole sums apart**: on the 3 Å pocket in water as built and on
/// W3's small box, each variant's `ΔU` pair by pair ([`cross_deltas`]) against `U_variant −
/// U_current` from two whole sums of their own ([`cross_sum`]), to the summation bound; variant 3's
/// classes, and its residues with the corrections, against its `ΔU`; the current model's whole sum
/// against the field's own coupling, and variant 0's against 5794556's [`Surroundings::van_der_waals`];
/// every variant's `ΔU` not zero; and in water, which has no protein, variant 3's `ΔU` variant 2's to
/// the bit, where on the pocket it is not.
#[test]
fn each_variants_change_is_two_whole_sums_apart() {
    for (k, (sur, labels, at)) in small_phases().into_iter().enumerate() {
        let current = current_model(&sur);
        let vs = variants(&sur, &labels);
        let mut models = vec![&current];
        models.extend(vs.iter());
        let ts = tails(&sur, &labels, &models);
        let f = frame_deltas(&sur, &labels, &at, &models, &ts);
        let whole = cross_sum(&sur, &at, &current, ts.whole[0]);
        let (field, _) = sur.coupling(&at);
        let to_field =
            (whole.total() - field).abs() / summation_bound(whole.terms, whole.magnitude);
        let hybrid = sur.van_der_waals(&at, &gaff_benzene(&ring_and_hydrogens(&sur)));
        let v0 = cross_sum(&sur, &at, &vs[0], ts.whole[1]);
        let to_hybrid = (hybrid.total() - v0.total()).abs()
            / summation_bound(hybrid.terms.max(v0.terms), hybrid.magnitude + v0.magnitude);
        println!(
            "{k}: ΔU {:?} kcal/mol; against the sums {:.2e}, parts {:.2e}, the field {:.2e}, \
             5794556's {:.2e} of their bounds",
            f.du.iter().map(|&u| kcal(u)).collect::<Vec<_>>(),
            f.against_sums,
            f.parts,
            to_field,
            to_hybrid
        );
        assert!(f.against_sums <= 1.0, "{k}: {}", f.against_sums);
        assert!(f.parts <= 1.0, "{k}: {}", f.parts);
        assert!(to_field <= 1.0, "{k}: {to_field}");
        assert!(to_hybrid <= 1.0, "{k}: {to_hybrid}");
        for (v, u) in f.du.iter().enumerate() {
            assert!(u.abs() > 1e-4 * KCAL_PER_MOL, "{k}: variant {v}, {u:e}");
        }
        if k == 0 {
            assert!(
                (f.du[3] - f.du[2]).abs() > 1e-4 * KCAL_PER_MOL,
                "the pocket: ff14SB's protein moves ΔU"
            );
        } else {
            assert_eq!(f.du[3].to_bits(), f.du[2].to_bits(), "water has no protein");
        }
    }
}

/// **Frames round-trip to the bit**: three frames of fourteen atoms, the forty-two numbers of
/// [`positions_round_trip_to_the_bit`] turned, subnormals, `±0.0` and the extremes among them, and
/// steps up to `u64::MAX`, through [`FrameWriter`] and [`read_frames`]: first line, steps and every
/// position's `to_bits` equal, nothing left at the `.partial` path; **the file cut short at every
/// byte** is refused; one bit flipped in any byte of a frame is refused; and a writer given fewer
/// frames than it promised refuses to finish.
#[cfg(not(target_family = "wasm"))]
#[test]
fn frames_round_trip_to_the_bit() {
    let mut xs = vec![
        f64::from_bits(1),
        f64::from_bits(0x000f_ffff_ffff_ffff),
        f64::MIN_POSITIVE,
        -0.0,
        0.0,
        f64::MAX,
        -f64::MAX,
        f64::EPSILON,
        1.0 / 3.0,
        0.1,
        -2.5e-10,
        1e-300,
    ];
    let mut state: u64 = 0x9E37_79B9_7F4A_7C15;
    while xs.len() < 42 {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        let x = f64::from_bits(state);
        if x.is_finite() {
            xs.push(x);
        }
    }
    let frames: Vec<Vec<[f64; 3]>> = (0..3)
        .map(|f| {
            (0..14)
                .map(|a| [0, 1, 2].map(|x| xs[(3 * a + x + 5 * f) % xs.len()]))
                .collect()
        })
        .collect();
    let steps = [0, 25, u64::MAX];
    let dir = fresh_dir("frames");
    let path = dir.join("test.frames");
    let mut w = FrameWriter::create(&path, "# frames a header", 3, 14);
    for (s, f) in steps.iter().zip(&frames) {
        w.push(*s, f);
    }
    w.finish();
    assert!(!path.with_extension("partial").exists());
    let back = read_frames(&path).expect("written").expect("whole");
    assert_eq!(back.line, frames_line("# frames a header", 3, 14));
    assert_eq!(back.steps, steps);
    let bits = |v: &[Vec<[f64; 3]>]| -> Vec<Vec<[u64; 3]>> {
        v.iter()
            .map(|f| f.iter().map(|p| p.map(f64::to_bits)).collect())
            .collect()
    };
    assert_eq!(bits(&back.at), bits(&frames));
    let bytes = std::fs::read(&path).unwrap();
    for cut in 0..bytes.len() {
        assert!(
            parse_frames(&mut &bytes[..cut]).is_err(),
            "cut at byte {cut} of {} read as whole",
            bytes.len()
        );
    }
    let start = bytes.iter().position(|&b| b == b'\n').unwrap() + 1;
    let end = bytes.len() - format!("# end 3 {:016x}\n", 0u64).len();
    for k in start..end {
        let mut flipped = bytes.clone();
        flipped[k] ^= 0x10;
        assert!(
            parse_frames(&mut &flipped[..]).is_err(),
            "a bit flipped at byte {k} read as whole"
        );
    }
    let short = refusal(|| {
        let mut w = FrameWriter::create(&dir.join("short.frames"), "# short", 2, 14);
        w.push(0, &frames[0]);
        w.finish();
    });
    assert!(short.contains("fewer frames"), "{short}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// **An apportioning resumes to the same bits**, through the measurement's own [`apportion`] on the
/// 3 Å pocket in water and W3's small box, three frames two 0.5 fs steps apart, in a temporary
/// directory: run, with a stale `complex.partial` left as a killed run leaves it; then run again,
/// which reads both phases' frames back, their files' bytes unchanged, **and every number it
/// measures is the first run's to the bit**; and a phase asked for at another protocol is refused by
/// name.
#[cfg(not(target_family = "wasm"))]
#[test]
fn an_apportioning_resumes_to_the_same_bits() {
    let s = small3(0x5EC5);
    let sur = complex_surroundings(&s);
    let labels = complex_labels(&s, the_system(), &sur);
    let decoupling = s.decoupling();
    let template = s.dynamics(LEG_FLEXIBILITY);
    let start = || s.positions().to_vec();
    let complex = Phase {
        name: "complex",
        sur: &sur,
        labels: &labels,
        decoupling: &decoupling,
        template: &template,
        lambda: Lambda::new(0.0, 1.0, 1.0),
        identity: format!("built {:016x}", positions_digest(s.positions())),
        start: &start,
    };
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC6, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    let free = water_surroundings(&b, &w);
    let free_labels = water_labels(&free);
    let (wd, wt) = (w.decoupling(), w.dynamics(&b));
    let water_start = || w.at.clone();
    let water = Phase {
        name: "water",
        sur: &free,
        labels: &free_labels,
        decoupling: &wd,
        template: &wt,
        lambda: Lambda::COUPLED,
        identity: format!("built {:016x}", positions_digest(&w.at)),
        start: &water_start,
    };
    let p = Protocol {
        time_step: 0.5 * FS,
        equilibration: 2,
        stride: 2,
        samples: 3,
        ..FRAMES
    };
    let dir = fresh_dir("apportion");
    std::fs::write(dir.join("complex.partial"), b"a killed run's frames").unwrap();
    let first = apportion(&dir, p, [&complex, &water], 1);
    assert!(!dir.join("complex.partial").exists());
    let files: Vec<Vec<u8>> = ["complex.frames", "water.frames"]
        .iter()
        .map(|f| std::fs::read(dir.join(f)).expect("saved"))
        .collect();
    let again = apportion(&dir, p, [&complex, &water], 1);
    for (f, bytes) in ["complex.frames", "water.frames"].iter().zip(&files) {
        assert_eq!(&std::fs::read(dir.join(f)).unwrap(), bytes, "{f} rewritten");
    }
    let log = std::fs::read_to_string(dir.join("variants.log")).unwrap();
    assert_eq!(log.matches("frames read back").count(), 2, "{log}");
    let bits = |r: &Reparametrised| -> Vec<u64> {
        r.du.iter()
            .chain(&r.u_variants)
            .chain(&r.classes)
            .chain(&r.residues)
            .flatten()
            .chain(&r.u)
            .chain([
                &r.against_sums,
                &r.parts,
                &r.against_field,
                &r.against_hybrid,
            ])
            .map(|x| x.to_bits())
            .collect()
    };
    for (a, b) in first.iter().zip(&again) {
        assert_eq!(a.du.len(), 3);
        assert_eq!(bits(a), bits(b));
        assert!(a.against_sums <= 1.0 && a.parts <= 1.0, "{a:?}");
        assert!(a.against_field <= 1.0 && a.against_hybrid <= 1.0, "{a:?}");
    }
    // The tables the measurement prints, on these three frames: what would otherwise first run
    // after hours of dynamics.
    report_apportioned(p.kt(), &first[0], &first[1], &labels);
    let other = Protocol {
        seed: p.seed + 1,
        ..p
    };
    let refused = refusal(|| {
        phase_frames(&dir, &complex, other);
    });
    assert!(refused.contains("was run differently"), "{refused}");
    let _ = std::fs::remove_dir_all(&dir);
}

// ---------------------------------------------------------------------------------------------
// The restraint interval refined: a window at λ_r 0.25

/// The candidate the refined window runs at: λ_r 0.25 with everything coupled, `grid()[4]`. **What
/// the insertion would choose first**: `refine_schedule` puts `⌊(0 + 8)/2⌋` = 4 between 0 and 8.
/// Chosen over λ_r 0.5 (`grid()[6]`) because the unrestrained end is the hard one: read from windows
/// 0 and 8 alone, one-sided estimates toward 4 kept 378 and 442 effective samples of 2000, toward 6
/// 247 and 1 263. A window at 4 is the one the first round would run, to the header and to the bit
/// (`the_refined_window_is_the_one_the_insertion_runs`).
const REFINED_CANDIDATE: usize = 4;

/// The schedule's windows read back, each refused unless its file is there — none is run — and the
/// window at [`REFINED_CANDIDATE`] read back or run, through [`window`]: resumable, written as
/// `window_04.txt`, its header the leg's.
fn refined_records(
    leg: &Leg,
    dir: &std::path::Path,
    schedule: &[usize],
) -> (Vec<Recorded>, Recorded) {
    for &g in schedule {
        let path = dir.join(format!("window_{g:02}.txt"));
        assert!(
            path.exists(),
            "{} is not there: the refined window reads the leg's windows and runs none of them",
            path.display()
        );
    }
    let records = schedule.iter().map(|&g| window(leg, dir, g)).collect();
    (records, window(leg, dir, REFINED_CANDIDATE))
}

/// The leg's records with the refined window in its place, path order.
fn with_refined(records: &[Recorded], refined: &Recorded) -> Vec<Recorded> {
    let mut all = records.to_vec();
    let at = all
        .iter()
        .position(|r| r.grid > refined.grid)
        .expect("a window after the refined one");
    assert!(at > 0 && all[at - 1].grid < refined.grid);
    all.insert(at, refined.clone());
    all
}

/// **The refined window's header is window 0's but for its candidate**: on the 3 Å pocket, the
/// leg's [`header`] at [`REFINED_CANDIDATE`] is its header at 0 with `candidate 0 lambda 0e0 1e0
/// 1e0` replaced by `candidate 4 lambda 2.5e-1 1e0 1e0`, and nothing else; and the protocol and the
/// state count in it are the ones `window_00.txt` of the production run carries, literally, so the
/// refined window can only be read back beside windows that ran under them.
#[test]
fn the_refined_windows_header_is_window_zeros_but_its_candidate() {
    assert_eq!(grid()[REFINED_CANDIDATE], Lambda::new(0.25, 1.0, 1.0));
    let s = small3(0x2E5);
    let leg = production_leg(&s, s.positions().to_vec());
    let zero = header(&leg, 0);
    let refined = header(&leg, REFINED_CANDIDATE);
    assert!(zero.contains(" protocol 2e-15 2.9815e2 1e12 10000 50 2000 281214 atoms "));
    assert!(zero.contains(" flexibility Mobile states 33 start "));
    assert_eq!(
        refined,
        zero.replacen(
            "# candidate 0 lambda 0e0 1e0 1e0 ",
            "# candidate 4 lambda 2.5e-1 1e0 1e0 ",
            1
        )
    );
    assert_ne!(refined, zero);
}

/// **The refined window is the one the insertion runs**, to the header and to the bit: on the 3 Å
/// pocket with [`a_killed_leg_resumes_to_the_same_bits`]'s short stages and windows, the leg's
/// [`refine_schedule`] from windows 0 and 8, told every interval overlaps too little, inserts 4,
/// then 2 and 6, through [`window`] — the call [`the_complex_leg_measured`] makes. In another
/// directory, windows 0 and 8 and then [`refined_records`] run the window at 4 alone. Its file's
/// header line is the inserted one's, and every gradient and energy the same bits. And
/// [`refined_records`] refuses a directory without the schedule's windows rather than run them.
#[cfg(not(target_family = "wasm"))]
#[test]
fn the_refined_window_is_the_one_the_insertion_runs() {
    let restraint_grid = || grid()[..9].to_vec();
    let s = small3(0x2E5);
    let coupled = s.decoupling();
    let potential = AtLambda {
        hamiltonian: &coupled,
        lambda: Lambda::COUPLED,
    };
    let stages = [
        Stage {
            flexibility: Flexibility::Frozen,
            time_step: 0.5e-15,
            steps: 3,
            friction: 50e12,
        },
        Stage {
            flexibility: LEG_FLEXIBILITY,
            time_step: 0.5e-15,
            steps: 2,
            friction: 50e12,
        },
    ];
    let protocol = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: 0x2E6,
    };
    let leg_in = |dir: &std::path::Path| {
        let start = prepared(
            &s,
            &potential,
            &stages,
            PREPARATION_SEED,
            KELVIN,
            dir,
            "prepared.txt",
            "complex",
        );
        leg_from(&s, start, protocol, restraint_grid())
    };
    let inserted = fresh_dir("inserted");
    let leg = leg_in(&inserted);
    let mut schedule = vec![0, 8];
    let mut records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| window(&leg, &inserted, g))
        .collect();
    refine_schedule(
        &mut schedule,
        &mut records,
        OVERLAP_THRESHOLD,
        INSERTION_ROUNDS,
        |_, _| 0.0,
        |g| window(&leg, &inserted, g),
    );
    assert_eq!(schedule, [0, 2, 4, 6, 8]);

    let alone = fresh_dir("refined");
    let leg = leg_in(&alone);
    let why = refusal(|| drop(refined_records(&leg, &alone, &[0, 8])));
    assert!(why.contains("runs none of them"), "{why}");
    for g in [0, 8] {
        window(&leg, &alone, g);
    }
    let (read, refined) = refined_records(&leg, &alone, &[0, 8]);
    assert_eq!(read.len(), 2);
    let first = |d: &std::path::Path| {
        std::fs::read_to_string(d.join("window_04.txt"))
            .unwrap()
            .lines()
            .next()
            .unwrap()
            .to_string()
    };
    assert_eq!(first(&alone), first(&inserted));
    assert_eq!(first(&alone), header(&leg, REFINED_CANDIDATE));
    let other = &records[2];
    assert_eq!(other.grid, REFINED_CANDIDATE);
    let bits = |r: &Recorded| -> Vec<u64> {
        let mut v: Vec<u64> = r.gradients.iter().flatten().map(|x| x.to_bits()).collect();
        for e in &r.energies {
            v.extend(e.iter().map(|x| x.to_bits()));
        }
        v
    };
    assert_eq!(bits(&refined), bits(other));
    assert_eq!(
        with_refined(&read, &refined)
            .iter()
            .map(|r| r.grid)
            .collect::<Vec<_>>(),
        [0, 4, 8]
    );
    for d in [inserted, alone] {
        let _ = std::fs::remove_dir_all(d);
    }
}

/// **The restraint interval refined: a window at λ_r 0.25** (release, one core, about nine hours,
/// one window). **Run once**, at `2b80d79`, exit 0 on 2026-10-16 after 36 185 s (a first attempt
/// was cut 25 minutes short by a Windows Update restart): 0 → 0.25 +1.314 ± 0.162 at overlap 0.198,
/// 0.25 → 1 +0.823 ± 0.013 at 0.412, the restraint on +2.137 ± 0.163 against the single interval's
/// +2.205 ± 0.187, and `ΔG°_bind` −8.975 ± 0.279. Nothing here asserts it. The leg's box and start
/// as [`the_complex_leg_measured`] has them —
/// [`whole_complex`], [`production_start`] reading `prepared.txt` back through its header check,
/// [`production_leg`] — in the same results directory; the fourteen windows of [`SCHEDULE`] read
/// back through [`window`]'s header check, none run; and the window at [`REFINED_CANDIDATE`] run
/// through [`window`], written to `window_04.partial` and renamed `window_04.txt` when complete,
/// so a stopped run resumes, its progress in `refined.log`. Before it runs, its header is checked
/// against `window_00.txt`'s first line with the candidate and λ replaced. Printed: BAR over
/// 0 → 0.25 → 1 with each interval's overlap, against the single interval's +2.205; the restraint-on
/// segment both ways; the complex leg and `ΔG°_bind` recomputed with the refined window.
/// `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p the_restraint_interval_refined -- --ignored --nocapture`.
#[test]
#[ignore = "one window of the complex leg at λ_r 0.25, about ten hours with --release, run once (2026-10-16)"]
fn the_restraint_interval_refined() {
    let dir = results_dir();
    let lock = Lock::take(&dir, "complex");
    log(
        &dir,
        "refined",
        &format!(
            "the refined window begins; its lock is {} — remove it by hand only once no run is \
             writing here",
            lock.0.display()
        ),
    );
    assert!(
        dir.join("prepared.txt").exists(),
        "no prepared.txt in {}: the refined window starts from the leg's",
        dir.display()
    );
    let t = std::time::Instant::now();
    let s = whole_complex();
    let mut leg = production_leg(&s, production_start(&s, &dir));
    leg.log = "refined";
    let zero = std::fs::read_to_string(dir.join("window_00.txt")).expect("window_00.txt");
    let zero = zero.lines().next().unwrap_or("");
    assert_eq!(
        header(&leg, REFINED_CANDIDATE),
        zero.replacen(
            "# candidate 0 lambda 0e0 1e0 1e0 ",
            "# candidate 4 lambda 2.5e-1 1e0 1e0 ",
            1
        ),
        "the refined window's header is window 0's but for its candidate"
    );
    let kt = BOLTZMANN.to_si() * KELVIN;
    let to_kcal = |x: f64| kcal(kt * x);
    let release = kcal(leg.restraint.release_free_energy(KELVIN));
    let (records, refined) = refined_records(&leg, &dir, &SCHEDULE);
    let all = with_refined(&records, &refined);
    for (name, path) in [
        ("one interval", vec![records[0].clone(), records[1].clone()]),
        (
            "refined",
            vec![records[0].clone(), refined.clone(), records[1].clone()],
        ),
    ] {
        for k in 0..path.len() - 1 {
            let bar = interval_bar(&path[k], &path[k + 1], kt);
            log(
                &dir,
                "refined",
                &format!(
                    "{name}: λ_r {} → {}: BAR {:+.3} ± {:.3} kcal/mol, overlap {:.3}",
                    path[k].lambda.restraint,
                    path[k + 1].lambda.restraint,
                    to_kcal(bar.delta),
                    to_kcal(bar.error()),
                    bar.overlap
                ),
            );
        }
        let c = chain(&path, kt);
        let (ti, te) = trapezoid(&path);
        log(
            &dir,
            "refined",
            &format!(
                "{name}: the restraint on, BAR {:+.3} ± {:.3}, TI {ti:+.3} ± {te:.3} kcal/mol",
                to_kcal(c.value),
                to_kcal(c.error)
            ),
        );
    }
    for (name, path) in [("one interval", &records), ("refined", &all)] {
        let total = chain(path, kt);
        let (bar, bar_error) = (to_kcal(total.value), to_kcal(total.error));
        let bind = SOLVENT_LEG.0 - bar - release;
        let error = (SOLVENT_LEG.1 * SOLVENT_LEG.1 + bar_error * bar_error).sqrt();
        log(
            &dir,
            "refined",
            &format!(
                "{name}: {} windows, complex leg BAR {bar:+.3} ± {bar_error:.3}; ΔG°_bind = \
                 {:+.3} − ({bar:+.3}) − ({release:+.3}) = {bind:+.3} ± {error:.3} kcal/mol \
                 (the single interval's run: +2.205 ± 0.187 for the restraint, −9.043 ± 0.294)",
                path.len(),
                SOLVENT_LEG.0
            ),
        );
    }
    log(
        &dir,
        "refined",
        &format!(
            "the refined window took {:.0} s, {:.0} s in all",
            refined.seconds,
            t.elapsed().as_secs_f64()
        ),
    );
}

// ---------------------------------------------------------------------------------------------
// QEq's benzene charges against FreeSolv's AM1-BCC, perturbed on the saved frames

/// The ligand of `b` as a graph for [`freesolv::charges_onto`]: its elements, and its bonds among
/// themselves in ligand indices. No bond leaves it.
fn ligand_graph(b: &Binding) -> (Vec<Element>, Vec<[usize; 2]>) {
    let ligand = b.ligand_range();
    let bonds = b.component().bonds();
    assert!(
        bonds
            .iter()
            .all(|bond| ligand.contains(&bond.atoms[0]) == ligand.contains(&bond.atoms[1])),
        "a bond leaves the ligand"
    );
    (
        b.elements()[ligand.clone()].to_vec(),
        bonds
            .iter()
            .filter(|bond| ligand.contains(&bond.atoms[0]))
            .map(|bond| bond.atoms.map(|i| i - ligand.start))
            .collect(),
    )
}

/// **FreeSolv's AM1-BCC charges on a benzene whose atoms are `graph`'s**, by element and bonding
/// ([`freesolv::charges_onto`]).
fn am1bcc(graph: &(Vec<Element>, Vec<[usize; 2]>)) -> freesolv::Mapped {
    freesolv::charges_onto(&freesolv::benzene(), &graph.0, &graph.1)
        .expect("FreeSolv's benzene is this one")
}

/// **The same charges neutral and equal by symmetry**: every carbon −0.1301 and every hydrogen
/// +0.1301, FreeSolv's C1's extra −0.0001 left out. A sensitivity, not FreeSolv's.
fn am1bcc_neutral(elements: &[Element]) -> Vec<f64> {
    elements
        .iter()
        .map(|e| match e {
            Element::C => -0.1301,
            Element::H => 0.1301,
            other => panic!("benzene has no {other:?}"),
        })
        .collect()
}

/// Benzene's cross Coulomb energy with everything else at `at`, joules: the field's own, the
/// couplings' `∂U/∂λ_e` at full coupling — its real-space pairs across, the reciprocal sum's cross
/// terms and the background's.
fn cross_coulomb(d: &PeriodicDecoupling, at: &[[f64; 3]]) -> f64 {
    d.couplings(at, &[Lambda::COUPLED])[0].gradient[1]
}

/// **Benzene's charges swapped in a box**: the decoupling as given, and the same field with
/// benzene's charges — the decoupled group's, in its order — replaced, every other atom's charge,
/// every van der Waals term, the mesh and the exclusions untouched.
struct ChargeSwap {
    from: PeriodicDecoupling,
    to: PeriodicDecoupling,
}

impl ChargeSwap {
    fn new(d: &PeriodicDecoupling, charges: &[f64]) -> ChargeSwap {
        let group = d.group();
        let members: Vec<usize> = (0..group.len()).filter(|&i| group[i]).collect();
        assert_eq!(members.len(), charges.len(), "one charge per benzene atom");
        let mut q = d.field().charges().to_vec();
        for (k, &i) in members.iter().enumerate() {
            q[i] = charges[k];
        }
        ChargeSwap {
            from: d.clone(),
            to: PeriodicDecoupling::new(&d.field().clone().with_charges(q), group)
                .expect("benzene decouples"),
        }
    }

    /// **`ΔU_elec` at `at`**, joules: benzene's cross Coulomb energy under the new charges less
    /// under the old ([`cross_coulomb`]). Benzene's own Coulomb pairs are not in it: the field
    /// computes them in vacuum at every λ, so they stay QEq's in the perturbed Hamiltonian, whose
    /// decoupled end state is then the same as the run's in both legs.
    fn delta(&self, at: &[[f64; 3]]) -> f64 {
        cross_coulomb(&self.to, at) - cross_coulomb(&self.from, at)
    }
}

/// **Benzene's own Coulomb energy, by hand**: every pair of its atoms the field does not exclude
/// (not 1-2, not 1-3), in vacuum, `k_e q_i q_j / r` on `at` as given — benzene is whole in every box
/// here; with the number of pairs and `Σ |term|`.
fn own_coulomb(d: &PeriodicDecoupling, at: &[[f64; 3]]) -> (f64, usize, f64) {
    let g = d.group();
    let q = d.field().charges();
    let members: Vec<usize> = (0..g.len()).filter(|&i| g[i]).collect();
    let excluded: std::collections::BTreeSet<[usize; 2]> = d
        .field()
        .excluded_pairs()
        .iter()
        .filter(|p| g[p[0]] && g[p[1]])
        .map(|&[i, j]| [i.min(j), i.max(j)])
        .collect();
    let (mut e, mut n, mut m) = (0.0, 0, 0.0);
    for (p, &i) in members.iter().enumerate() {
        for &j in &members[p + 1..] {
            if excluded.contains(&[i, j]) {
                continue;
            }
            let u = COULOMB * q[i] * q[j] / len(sub(at[i], at[j]));
            e += u;
            n += 1;
            m += u.abs();
        }
    }
    (e, n, m)
}

/// One frame's swap held to the field.
#[derive(Clone, Copy, Debug)]
struct Held {
    /// `ΔU_elec`, the couplings'.
    du: f64,
    /// Benzene's own Coulomb energy's change, by hand ([`own_coulomb`]).
    own: f64,
    own_pairs: usize,
    /// `|ΔU_elec − [(E_to − E_from) − ΔU_own]|` over `64 ε` of both evaluations' parts and the own
    /// pairs' — W1's and W3's bound for the cross terms through the mesh.
    ratio: f64,
    /// That bound's scale: both evaluations' parts and the own pairs', joules.
    scale: f64,
}

/// **The swap held to two whole evaluations** at `at`: each decoupling's electrostatic energy at
/// full coupling — the rest's Ewald sum, benzene's cross terms through the evaluation's own pair
/// search and mesh, and benzene's own pairs in vacuum — differenced, benzene's own pairs' change
/// by hand taken off, against [`ChargeSwap::delta`].
fn held(swap: &ChargeSwap, at: &[[f64; 3]]) -> Held {
    let du = swap.delta(at);
    let to = swap.to.energy(at, Lambda::COUPLED);
    let from = swap.from.energy(at, Lambda::COUPLED);
    let (own_to, own_pairs, m_to) = own_coulomb(&swap.to, at);
    let (own_from, _, m_from) = own_coulomb(&swap.from, at);
    let whole = (to.electrostatic - from.electrostatic) - (own_to - own_from);
    let scale = magnitude(&to) + magnitude(&from) + m_to + m_from;
    Held {
        du,
        own: own_to - own_from,
        own_pairs,
        ratio: (du - whole).abs() / (64.0 * EPS * scale),
        scale,
    }
}

/// A box the swap is checked on: its name, the box taken apart, benzene's graph ([`ligand_graph`])
/// and the positions.
type SwapBox = (
    &'static str,
    Surroundings,
    (Vec<Element>, Vec<[usize; 2]>),
    Vec<[f64; 3]>,
);

/// The boxes the swap is checked on, taken apart ([`Surroundings`]) with benzene's graph and the
/// positions as built: the 3 Å pocket in water, whose reciprocal sum is the mesh's, and W3's small
/// box, classical Ewald.
fn swap_boxes() -> Vec<SwapBox> {
    let s = small3(0x5EC5);
    let pocket = complex_surroundings(&s);
    assert!(
        pocket.whole.field().ewald().mesh().is_some(),
        "the pocket's mesh"
    );
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC6, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    let free = water_surroundings(&b, &w);
    assert!(
        free.whole.field().ewald().mesh().is_none(),
        "W3's classical sum"
    );
    vec![
        (
            "3 Å pocket, mesh",
            pocket,
            ligand_graph(s.binding()),
            s.positions().to_vec(),
        ),
        (
            "W3's small box, Ewald",
            free,
            freesolv::component_graph(&b.component),
            w.at.clone(),
        ),
    ]
}

/// **A `.mol2` is read with its charge column**: a literal water, no disk — its name, every atom's
/// name, element, type, position and charge to the bit, and its bonds as indices — and each thing
/// the reader refuses, refused by name; then FreeSolv's benzene as kept ([`freesolv`]), twelve
/// atoms of GAFF's `ca` and `ha` and twelve bonds.
#[test]
fn a_mol2_is_read_with_its_charges() {
    let water =
        "@<TRIPOS>MOLECULE\nwater\n    3     2     1     0     0\nSMALL\nUSER_CHARGES\n\n\n\
                 @<TRIPOS>ATOM\n\
                 \x20     1 O1          0.0000    0.0000    0.0000 ow        1 WAT     -0.834000\n\
                 \x20     2 H1          0.9572    0.0000    0.0000 hw        1 WAT      0.417000\n\
                 \x20     3 H2         -0.2400    0.9266    0.0000 hw        1 WAT      0.417000\n\
                 @<TRIPOS>BOND\n\
                 \x20    1    1    2 1\n\
                 \x20    2    3    1 1\n\
                 @<TRIPOS>SUBSTRUCTURE\n\
                 \x20    1 WAT         1 TEMP              0 ****  ****    0 ROOT\n";
    let m = freesolv::parse(water).expect("a water");
    assert_eq!(m.name, "water");
    let names: Vec<&str> = m.atoms.iter().map(|a| a.name.as_str()).collect();
    assert_eq!(names, ["O1", "H1", "H2"]);
    let elements: Vec<Element> = m.atoms.iter().map(|a| a.element).collect();
    assert_eq!(elements, [Element::O, Element::H, Element::H]);
    let kinds: Vec<&str> = m.atoms.iter().map(|a| a.kind.as_str()).collect();
    assert_eq!(kinds, ["ow", "hw", "hw"]);
    let charges: Vec<u64> = m.atoms.iter().map(|a| a.charge.to_bits()).collect();
    assert_eq!(
        charges,
        [-0.834f64, 0.417, 0.417].map(f64::to_bits),
        "the charge column"
    );
    assert_eq!(
        m.atoms[2].at.map(f64::to_bits),
        [-0.24f64, 0.9266, 0.0].map(f64::to_bits)
    );
    assert_eq!(m.bonds, [[0, 1], [2, 0]]);

    let refused = |text: String, why: &str| {
        let e = freesolv::parse(&text).expect_err(why);
        assert!(e.contains(why), "{e:?} does not say {why:?}");
    };
    refused(
        water.replace("WAT      0.417000\n      3", "\n      3"),
        "columns",
    );
    refused(water.replace("-0.834000", "-0.83x"), "O1's charge");
    refused(
        water.replace("    3     2  ", "    4     2  "),
        "3 atoms where",
    );
    refused(
        water.replace("    3     2  ", "    3     3  "),
        "2 bonds where",
    );
    refused(
        water.replace("    2    3    1 1", "    2    4    1 1"),
        "not there",
    );
    refused(
        water.replace("    2    3    1 1", "    2    2    1 1"),
        "repeats",
    );
    refused(
        water.replace("    2    3    1 1", "    2    1    1 1"),
        "itself",
    );
    refused(water.replace("      2 H1", "      7 H1"), "numbered 7");
    refused(water.replace("O1  ", "Xx1 "), "not an element");
    refused(water.replace("0.9572", "inf"), "finite");
    refused(format!("{water}{water}"), "a second molecule");

    let b = freesolv::benzene();
    assert_eq!(b.name, "benzene");
    assert_eq!((b.atoms.len(), b.bonds.len()), (12, 12));
    for a in &b.atoms {
        let want = if a.element == Element::C { "ca" } else { "ha" };
        assert_eq!(a.kind, want, "{}", a.name);
    }
}

/// A benzene written out of order: carbons and hydrogens interleaved, the bonds shuffled, each
/// carbon's charge its own (−0.1301 to −0.1306 by its name) and each hydrogen's minus its
/// carbon's. Mapped by file order, the crate's C1 would take a hydrogen's charge.
const SHUFFLED_BENZENE: &str =
    "@<TRIPOS>MOLECULE\nshuffled\n   12    12     1     0     0\nSMALL\n\
     USER_CHARGES\n\n\n@<TRIPOS>ATOM\n\
     \x20 1 H4  0 0 0 ha 1 MOL  0.1304\n\
     \x20 2 C2  0 0 0 ca 1 MOL -0.1302\n\
     \x20 3 H1  0 0 0 ha 1 MOL  0.1301\n\
     \x20 4 C5  0 0 0 ca 1 MOL -0.1305\n\
     \x20 5 C1  0 0 0 ca 1 MOL -0.1301\n\
     \x20 6 H6  0 0 0 ha 1 MOL  0.1306\n\
     \x20 7 C4  0 0 0 ca 1 MOL -0.1304\n\
     \x20 8 H3  0 0 0 ha 1 MOL  0.1303\n\
     \x20 9 C6  0 0 0 ca 1 MOL -0.1306\n\
     \x2010 H2  0 0 0 ha 1 MOL  0.1302\n\
     \x2011 C3  0 0 0 ca 1 MOL -0.1303\n\
     \x2012 H5  0 0 0 ha 1 MOL  0.1305\n\
     @<TRIPOS>BOND\n\
     \x20 1  7  1 1\n\x20 2  2 11 ar\n\x20 3  9  5 ar\n\x20 4  5  3 1\n\
     \x20 5 11  7 ar\n\x20 6  4 12 1\n\x20 7  5  2 ar\n\x20 8  9  6 1\n\
     \x20 9  2 10 1\n\x2010  4  9 ar\n\x2011 11  8 1\n\x2012  7  4 ar\n";

/// **FreeSolv's charges map by bonding, not by file order.** [`SHUFFLED_BENZENE`] onto the
/// dictionary's `BNZ`: twelve automorphisms; every carbon takes a carbon's charge and every
/// hydrogen minus its own carbon's; every bond of `BNZ` is a bond of the file. Then FreeSolv's own
/// file, onto `BNZ` and onto the 3 Å pocket's benzene: every carbon −0.1302 or −0.1301 and one
/// of them −0.1302, every hydrogen +0.1301, a spread of 0.0001 over the automorphisms, and a sum
/// of −0.0001 e to its rounding.
#[test]
fn freesolvs_charges_map_by_bonding_not_by_file_order() {
    let bnz = Component::from_ccd(hydration::BNZ).expect("BNZ");
    let graph = freesolv::component_graph(&bnz);
    let shuffled = freesolv::parse(SHUFFLED_BENZENE).expect("the shuffled benzene");
    let m = freesolv::charges_onto(&shuffled, &graph.0, &graph.1).expect("benzene onto benzene");
    assert_eq!(m.automorphisms, 12);
    let carbons = [-0.1301, -0.1302, -0.1303, -0.1304, -0.1305, -0.1306];
    for (i, a) in bnz.atoms().iter().enumerate() {
        let q = m.charges[i];
        match a.element {
            Element::C => assert!(carbons.contains(&q), "{}: {q}", a.name),
            Element::H => {
                let (c, _) = bnz.neighbours(i).next().expect("a hydrogen's carbon");
                assert_eq!(q, -m.charges[c], "{}", a.name);
            }
            other => panic!("{other:?}"),
        }
    }
    for b in bnz.bonds() {
        let [x, y] = b.atoms.map(|i| m.map[i]);
        assert!(
            shuffled.bonds.contains(&[x, y]) || shuffled.bonds.contains(&[y, x]),
            "{:?}",
            b.atoms
        );
    }

    for (name, graph) in [
        ("BNZ", graph.clone()),
        ("the 3 Å pocket's", ligand_graph(pocket3())),
    ] {
        let real = am1bcc(&graph);
        let (mut c1, mut sum, mut size) = (0, 0.0, 0.0);
        for (e, &q) in graph.0.iter().zip(&real.charges) {
            match e {
                Element::C => {
                    assert!(q == -0.1301 || q == -0.1302, "{name}: {q}");
                    c1 += usize::from(q == -0.1302);
                }
                Element::H => assert_eq!(q, 0.1301, "{name}"),
                other => panic!("{other:?}"),
            }
            sum += q;
            size += q.abs();
        }
        assert_eq!(c1, 1, "{name}: one carbon at −0.1302");
        assert_eq!(real.automorphisms, 12, "{name}");
        assert!(
            (real.spread - 1e-4).abs() <= 4.0 * EPS * 0.1302,
            "{name}: {}",
            real.spread
        );
        assert!(
            (sum + 1e-4).abs() <= 16.0 * EPS * size,
            "{name}: Σq = {sum:e}"
        );
        println!(
            "{name}: {:?}, Σq {sum:+.3e}, spread {:.1e} over {} automorphisms",
            real.charges, real.spread, real.automorphisms
        );
    }
}

/// **Benzene's charges swapped for its own change nothing, to the bit**: on the 3 Å pocket in
/// water and on W3's small box, the swap ([`ChargeSwap`]) to benzene's own QEq charges gives
/// `ΔU_elec` = 0 exactly, whole and in every class ([`CLASSES`]) with an atom; and the cross Coulomb
/// energy it differences is not zero.
#[test]
fn benzenes_charges_swapped_for_its_own_change_nothing() {
    for (name, sur, _, at) in swap_boxes() {
        let own: Vec<f64> = sur
            .ligand
            .clone()
            .map(|i| sur.whole.field().charges()[i])
            .collect();
        let mut swaps = vec![ChargeSwap::new(&sur.whole, &own)];
        swaps.extend(
            sur.classes
                .iter()
                .flatten()
                .map(|d| ChargeSwap::new(d, &own)),
        );
        assert!(swaps.len() >= 2, "{name}: a class");
        let du: Vec<f64> = swaps.iter().map(|s| s.delta(&at)).collect();
        println!("{name}: ΔU_elec {du:?} J, whole and by class");
        for (k, u) in du.iter().enumerate() {
            assert_eq!(u.to_bits(), 0, "{name}: swap {k}");
        }
        assert!(
            cross_coulomb(&sur.whole, &at).abs() > 1e-3 * KCAL_PER_MOL,
            "{name}: a cross Coulomb energy"
        );
    }
}

/// **Benzene's charges swapped are two whole evaluations apart**: on the 3 Å pocket in water through
/// the mesh and on W3's small box through the classical sum, FreeSolv's charges for QEq's,
/// `ΔU_elec` against the two decouplings' whole electrostatic energies at full coupling less benzene's
/// own 36 pairs' change by hand, within 64 ε of their parts ([`held`]); the classes' swaps summing to
/// the whole's within the same; and neither `ΔU_elec` nor the own pairs' change zero. Both boxes
/// are printed before either is asserted.
#[test]
fn benzenes_charges_swapped_are_two_whole_evaluations_apart() {
    let mut results = Vec::new();
    for (name, sur, graph, at) in swap_boxes() {
        let q = am1bcc(&graph).charges;
        let swap = ChargeSwap::new(&sur.whole, &q);
        let h = held(&swap, &at);
        let classes: f64 = sur
            .classes
            .iter()
            .flatten()
            .map(|d| ChargeSwap::new(d, &q).delta(&at))
            .sum();
        let by_class = (classes - h.du).abs() / (64.0 * EPS * h.scale);
        println!(
            "{name}: ΔU_elec {:+.5} kcal/mol, {:.2e} of the bound from the evaluations; the \
             classes' {:.2e}; benzene's own {} pairs {:+.5}",
            kcal(h.du),
            h.ratio,
            by_class,
            h.own_pairs,
            kcal(h.own)
        );
        results.push((name, h, by_class));
    }
    for (name, h, by_class) in results {
        assert_eq!(
            h.own_pairs, 36,
            "{name}: 66 pairs less 12 bonds and 18 angles"
        );
        assert!(h.ratio <= 1.0, "{name}: {h:?}");
        assert!(by_class <= 1.0, "{name}: {by_class}");
        assert!(h.du.abs() > 1e-3 * KCAL_PER_MOL, "{name}: {h:?}");
        assert!(h.own.abs() > 1e-3 * KCAL_PER_MOL, "{name}: {h:?}");
    }
}

/// [`Phase`]'s frames **read back, never run**: `NAME.frames` in `dir` under the header the phase
/// would have written ([`phase_header`]), refused by name when it is absent, not whole, or was run
/// differently. Nothing is written.
fn saved_frames(dir: &std::path::Path, phase: &Phase, p: Protocol) -> Frames {
    let path = dir.join(format!("{}.frames", phase.name));
    let want = frames_line(&phase_header(phase, p), p.samples, phase.sur.kinds.len());
    let f = read_frames(&path)
        .unwrap_or_else(|| {
            panic!(
                "{} is not there: this measurement runs no dynamics",
                path.display()
            )
        })
        .unwrap_or_else(|why| panic!("{} is not whole: {why}", path.display()));
    assert_eq!(f.line, want, "{} was run differently", path.display());
    assert!(f.steps.windows(2).all(|s| s[0] < s[1]), "steps that rise");
    f
}

/// The leg's prepared start **read back, never prepared**: `prepared.txt` in `dir` under
/// [`prepared_header`]'s header, refused by name when it is absent, not whole, or prepared
/// differently.
fn saved_start(s: &SolvatedComplex, dir: &std::path::Path) -> Vec<[f64; 3]> {
    let path = dir.join("prepared.txt");
    let want = prepared_header(
        s,
        &SolvatedComplex::preparation(LEG_FLEXIBILITY),
        PREPARATION_SEED,
        KELVIN,
    );
    let (head, at) = read_positions(&path)
        .unwrap_or_else(|| panic!("{} is not there", path.display()))
        .unwrap_or_else(|why| panic!("{} is not whole: {why}", path.display()));
    assert_eq!(head, want, "{} was prepared differently", path.display());
    assert_eq!(at.len(), s.positions().len());
    at
}

/// A phase's frames with benzene's charges swapped: per frame `ΔU_elec` for FreeSolv's charges and
/// for the neutral set ([`am1bcc_neutral`]), and FreeSolv's by class; every `every`th frame benzene's
/// own Coulomb change and the swap held to two whole evaluations ([`held`]); and the worst of each
/// identity, over its bound.
struct SwapFrames {
    elec: Vec<f64>,
    neutral: Vec<f64>,
    classes: Vec<[f64; 3]>,
    own: Vec<f64>,
    against_classes: f64,
    against_evaluations: f64,
}

// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn swap_frames(
    sur: &Surroundings,
    frames: &Frames,
    charges: [&[f64]; 2],
    every: usize,
    name: &str,
    dir: &std::path::Path,
) -> SwapFrames {
    let swap = ChargeSwap::new(&sur.whole, charges[0]);
    let neutral = ChargeSwap::new(&sur.whole, charges[1]);
    let classes: Vec<Option<ChargeSwap>> = sur
        .classes
        .iter()
        .map(|d| d.as_ref().map(|d| ChargeSwap::new(d, charges[0])))
        .collect();
    let mut out = SwapFrames {
        elec: Vec::new(),
        neutral: Vec::new(),
        classes: Vec::new(),
        own: Vec::new(),
        against_classes: 0.0,
        against_evaluations: 0.0,
    };
    let t = std::time::Instant::now();
    let mut scale = 0.0;
    for (k, at) in frames.at.iter().enumerate() {
        if k % every == 0 {
            let h = held(&swap, at);
            assert_eq!(h.own_pairs, 36, "{name}: benzene's own pairs");
            scale = h.scale;
            out.own.push(h.own);
            out.against_evaluations = out.against_evaluations.max(h.ratio);
        }
        let du = swap.delta(at);
        let mut c = [0.0; 3];
        for (j, s) in classes.iter().enumerate() {
            if let Some(s) = s {
                c[j] = s.delta(at);
            }
        }
        let sum: f64 = c.iter().sum();
        out.against_classes = out
            .against_classes
            .max((sum - du).abs() / (64.0 * EPS * scale));
        out.elec.push(du);
        out.neutral.push(neutral.delta(at));
        out.classes.push(c);
        if (k + 1) % 80 == 0 {
            log(
                dir,
                "charges",
                &format!(
                    "{name}: {} frames swapped, {:.0} s",
                    k + 1,
                    t.elapsed().as_secs_f64()
                ),
            );
        }
    }
    out
}

/// The free energy of `ΔU` along `U₀ + s ΔU` at `stages + 1` evenly spaced `s` **from the same
/// frames**, a one-sided chain: each stage's `−kT ln(⟨e^(−β s_(k+1) ΔU)⟩ / ⟨e^(−β s_k ΔU)⟩)`, kcal/mol,
/// with the effective samples of the frames reweighted to `s_(k+1)`. **It telescopes**: the stages sum
/// to the one-step estimate, and nothing in it reaches configurations the frames do not have. For
/// the charges, `ΔU` is linear in `s` exactly — benzene's cross Coulomb energy is linear in its
/// charges — so these are the states a staged run would sample.
fn staged(du: &[f64], kt: f64, stages: usize) -> Vec<(f64, f64, f64)> {
    let at = |s: f64| -> (f64, f64) {
        let b: Vec<f64> = du.iter().map(|u| s * u / kt).collect();
        let m = b.iter().copied().fold(f64::INFINITY, f64::min);
        let w: Vec<f64> = b.iter().map(|y| (-(y - m)).exp()).collect();
        let sum: f64 = w.iter().sum();
        let ess = sum * sum / w.iter().map(|v| v * v).sum::<f64>();
        (kcal(kt * (m - (sum / du.len() as f64).ln())), ess)
    };
    let mut out = Vec::new();
    let mut last = 0.0;
    for k in 1..=stages {
        let s = k as f64 / stages as f64;
        let (g, ess) = at(s);
        out.push((s, g - last, ess));
        last = g;
    }
    out
}

/// **The staged chain telescopes to the one step**: on 320 made-up `ΔU` spread over ±2.4 kT,
/// [`staged`]'s five stages sum to [`zwanzig`]'s one-step estimate within 64 ε of their parts, the
/// last state is `s` = 1 with Zwanzig's effective samples, and the first stage is not the whole.
#[test]
fn the_staged_chain_telescopes_to_one_step() {
    let kt = BOLTZMANN.to_si() * KELVIN;
    let du: Vec<f64> = (0..320usize)
        .map(|k| kt * 0.4 * ((k * 37 % 13) as f64 - 6.0))
        .collect();
    let z = zwanzig(&du, kt);
    let st = staged(&du, kt, 5);
    let sum: f64 = st.iter().map(|x| x.1).sum();
    let scale = st.iter().map(|x| x.1.abs()).sum::<f64>() + z.exp.abs();
    println!(
        "stages {st:?}, sum {sum:+.6} against one step {:+.6}",
        z.exp
    );
    assert_eq!(st.len(), 5);
    assert_eq!(st[4].0, 1.0);
    assert!(
        (sum - z.exp).abs() <= 64.0 * EPS * scale,
        "{sum} against {}",
        z.exp
    );
    assert!((st[4].2 - z.ess).abs() <= 64.0 * EPS * z.ess, "{st:?}");
    assert!((st[0].1 - z.exp).abs() > 1e-3, "{st:?}");
}

/// W3's hydration free energy, kcal/mol: minus its decoupling ([`SOLVENT_LEG`]).
fn w3_hydration() -> (f64, f64) {
    (-SOLVENT_LEG.0, SOLVENT_LEG.1)
}

/// `ΔG°_bind` of benzene in L99A by Mobley et al. 2007 (*J. Mol. Biol.* **371**, 1118, Table 3,
/// read in PMC2104542) with GAFF and AM1-BCC charges on the ligand and AMBER96 on the protein,
/// kcal/mol, and its error: −3.95 ± 0.20, 1.24 short of experiment's −5.19 (their AM1-CM2 Table 1
/// gives −4.56 ± 0.20).
const MOBLEY_2007: (f64, f64) = (-3.95, 0.20);

/// One perturbation's summary: `ΔΔG_bind` to first order, by the second cumulant and by Zwanzig
/// with its error and halves; σ(βΔU) and the effective samples bound and in water; the water side's
/// Zwanzig with its error.
struct SwapSummary {
    name: &'static str,
    first: f64,
    cumulant: f64,
    exp: f64,
    error: f64,
    halves: [f64; 2],
    spread: [f64; 2],
    ess: [f64; 2],
    water: (f64, f64),
}

/// The table for one perturbation, bound and in water — `⟨ΔU⟩`, σ(βΔU), the second cumulant,
/// Zwanzig ± its error, its halves and the effective samples — then `ΔΔG_bind`, the implied
/// `ΔG_hyd` and `ΔG°_bind`; and, when σ(βΔU) is above [`ONE_STEP_SPREAD`] on either side, the
/// staged chain ([`staged`]) with four intermediate states.
fn report_swap(name: &'static str, kt: f64, bound: &[f64], water: &[f64]) -> SwapSummary {
    let n = bound.len();
    let (zb, zw) = (zwanzig(bound, kt), zwanzig(water, kt));
    println!("\n**{name}**\n\n| | bound | in water |\n| --- | --- | --- |");
    println!("| ⟨ΔU⟩ | {} | {} |", pm(&zb.mean), pm(&zw.mean));
    println!("| σ(βΔU) | {:.3} | {:.3} |", zb.spread, zw.spread);
    println!(
        "| ⟨ΔU⟩ − βσ²/2 | {:+.3} | {:+.3} |",
        zb.cumulant, zw.cumulant
    );
    println!(
        "| Zwanzig −kT ln⟨e^(−βΔU)⟩ | {:+.3} ± {:.3} | {:+.3} ± {:.3} |",
        zb.exp, zb.exp_error, zw.exp, zw.exp_error
    );
    println!(
        "| Zwanzig, halves | {:+.3}, {:+.3} | {:+.3}, {:+.3} |",
        zb.halves[0], zb.halves[1], zw.halves[0], zw.halves[1]
    );
    println!(
        "| effective samples | {:.1} of {n} | {:.1} of {} |",
        zb.ess,
        zw.ess,
        water.len()
    );
    let exp = zb.exp - zw.exp;
    let error = (zb.exp_error * zb.exp_error + zw.exp_error * zw.exp_error).sqrt();
    let (hyd, hyd_error) = w3_hydration();
    let hydration = hyd + zw.exp;
    let bind = W4_BINDING.0 + exp;
    println!(
        "\nΔΔG_bind: first order {:+.3}, second cumulant {:+.3}, Zwanzig {exp:+.3} ± {error:.3} \
         (halves {:+.3}, {:+.3}). Implied ΔG_hyd {hyd:+.2} + ({:+.3}) = {hydration:+.3} ± {:.3}; \
         ΔG°_bind {:+.3} + ({exp:+.3}) = {bind:+.3} ± {:.3} kcal/mol.",
        zb.mean.mean - zw.mean.mean,
        zb.cumulant - zw.cumulant,
        zb.halves[0] - zw.halves[0],
        zb.halves[1] - zw.halves[1],
        zw.exp,
        (hyd_error * hyd_error + zw.exp_error * zw.exp_error).sqrt(),
        W4_BINDING.0,
        (W4_BINDING.1 * W4_BINDING.1 + error * error).sqrt()
    );
    for (side, z, du) in [("bound", &zb, bound), ("in water", &zw, water)] {
        if z.spread > ONE_STEP_SPREAD {
            println!(
                "\n**σ(βΔU) {:.3} {side}, above {ONE_STEP_SPREAD}: one step is not to be trusted.** \
                 The chain through U₀ + s ΔU from the same frames, four intermediate states — it \
                 telescopes to the one step and adds no coverage:\n\n| s | stage (kcal/mol) | \
                 effective samples at s |\n| --- | --- | --- |",
                z.spread
            );
            let mut sum = 0.0;
            for (s, g, ess) in staged(du, kt, 5) {
                sum += g;
                println!("| {s:.1} | {g:+.3} | {ess:.1} of {} |", du.len());
            }
            println!("| sum | {sum:+.3} (one step {:+.3}) | |", z.exp);
        }
    }
    SwapSummary {
        name,
        first: zb.mean.mean - zw.mean.mean,
        cumulant: zb.cumulant - zw.cumulant,
        exp,
        error,
        halves: [zb.halves[0] - zw.halves[0], zb.halves[1] - zw.halves[1]],
        spread: [zb.spread, zw.spread],
        ess: [zb.ess, zw.ess],
        water: (zw.exp, zw.exp_error),
    }
}

/// **QEq's benzene charges against FreeSolv's AM1-BCC, on the saved frames** (release, one core,
/// minutes: the whole complex is built and the frames read, and no dynamics is run). Zwanzig's
/// one-step perturbation from UFF/QEq's ensemble, bound and in water, of:
///
/// - **e**: benzene's cross Coulomb energy with everything else, QEq's ±0.098 e for FreeSolv's
///   AM1-BCC charges ([`freesolv`], mapped by element and bonding), through the field's own Ewald
///   and mesh: `ΔU_elec` ([`ChargeSwap::delta`]), its reciprocal-space and background terms
///   included. Benzene's own Coulomb pairs are left QEq's ([`ChargeSwap::delta`] says why), and
///   their change is printed;
/// - **e′**: the same with the neutral, symmetric set ([`am1bcc_neutral`]), a sensitivity;
/// - **3**: 0ff5194's variant 3, GAFF's benzene and ff14SB's protein van der Waals by
///   Lorentz–Berthelot, again on the same frames, for reference;
/// - **4**: variant 3's van der Waals and FreeSolv's charges together, `ΔU = ΔU_vdW(3) + ΔU_elec`
///   per frame: GAFF/AM1-BCC benzene against an ff14SB-van-der-Waals protein with UFF's bonded
///   terms and QEq's charges.
///
/// Each with `⟨ΔU⟩`, σ(βΔU), the second cumulant, Zwanzig with its error and halves, and the
/// effective samples ([`report_swap`]); `ΔΔG_bind`, the implied `ΔG_hyd` (W3's −1.57 plus the
/// water side) and `ΔG°_bind` (W4's −8.975 plus `ΔΔG_bind`); and **e** by class bound. Every
/// frame's class sum and every twentieth frame's two whole evaluations held to their bounds, and
/// [`reparametrise`]'s identities, asserted after everything is printed.
///
/// **Nothing is written beside the leg**: `prepared.txt` and the frames are read back through their
/// headers and refused if absent ([`saved_start`], [`saved_frames`]); the progress goes to
/// `frames/charges.log`.
///
/// **Run once**, on `3b5aeb7`'s tree with this change, exit 0 on 2026-10-19 after 346 s, the
/// frames' SHA-256 the same after as before. Zwanzig's `ΔΔG_bind`: 3 +3.365 ± 0.073 (0ff5194's,
/// reproduced), e +0.569 ± 0.104, e′ +0.570 ± 0.103, and 4 +3.953 ± 0.113, which puts `ΔG°_bind`
/// at −5.02 against experiment's −5.19 and the implied `ΔG_hyd` at −0.62 against −0.90. σ(βΔU) is
/// at most 1.12, so no variant is flagged. Nothing here asserts it; the CHANGELOG has the tables.
///
/// `PANTOMETRY_W4_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_bound_in_tip3p benzenes_charges_perturbed_measured -- --ignored --nocapture`.
#[test]
#[ignore = "W4's saved frames, read only; about six minutes with --release, run once (2026-10-19)"]
fn benzenes_charges_perturbed_measured() {
    let dir = results_dir();
    let frames_dir = dir.join("frames");
    for f in ["complex.frames", "water.frames"] {
        assert!(
            frames_dir.join(f).exists(),
            "{} is not there: this measurement runs no dynamics",
            frames_dir.join(f).display()
        );
    }
    let t = std::time::Instant::now();
    let kt = BOLTZMANN.to_si() * KELVIN;
    let s = whole_complex();
    let leg = production_leg(&s, saved_start(&s, &dir));
    log(
        &frames_dir,
        "charges",
        &format!(
            "the leg's box and start read back in {:.0} s",
            t.elapsed().as_secs_f64()
        ),
    );
    let g = grid();
    assert_eq!(g[0], Lambda::new(0.0, 1.0, 1.0));
    let no_dynamics = || -> Vec<[f64; 3]> { panic!("this measurement runs no dynamics") };
    let bound = complex_surroundings(&s);
    let bound_labels = complex_labels(&s, the_system(), &bound);
    let complex = Phase {
        name: "complex",
        sur: &bound,
        labels: &bound_labels,
        decoupling: &leg.decoupling,
        template: &leg.template,
        lambda: g[0],
        identity: format!(
            "start {:016x} restraint {:016x} flexibility {:?}",
            positions_digest(&leg.start),
            restraint_digest(&leg.restraint),
            leg.flexibility
        ),
        start: &no_dynamics,
    };
    let benzene = hydration::Benzene::new(Some(pocket3().ligand_force_field().charges().to_vec()));
    let w3 = hydration::w3_solvated(&benzene);
    let w3_decoupling = w3.decoupling();
    let w3_template = w3.dynamics(&benzene);
    let free = water_surroundings(&benzene, &w3);
    let free_labels = water_labels(&free);
    let water = Phase {
        name: "water",
        sur: &free,
        labels: &free_labels,
        decoupling: &w3_decoupling,
        template: &w3_template,
        lambda: Lambda::COUPLED,
        identity: format!(
            "built {:016x} melted w3_melt {}",
            positions_digest(&w3.at),
            exact(KELVIN)
        ),
        start: &no_dynamics,
    };
    let complex_frames = saved_frames(&frames_dir, &complex, FRAMES);
    let water_frames = saved_frames(&frames_dir, &water, FRAMES);
    log(
        &frames_dir,
        "charges",
        &format!(
            "{} and {} frames read back, {:.0} s",
            complex_frames.at.len(),
            water_frames.at.len(),
            t.elapsed().as_secs_f64()
        ),
    );

    // The charges, each benzene's in its own order.
    let graph_b = ligand_graph(s.binding());
    let graph_w = freesolv::component_graph(&benzene.component);
    let (qb, qw) = (am1bcc(&graph_b), am1bcc(&graph_w));
    let (nb, nw) = (am1bcc_neutral(&graph_b.0), am1bcc_neutral(&graph_w.0));
    let qeq_b: Vec<f64> = bound
        .ligand
        .clone()
        .map(|i| bound.whole.field().charges()[i])
        .collect();
    let qeq_w: Vec<f64> = free
        .ligand
        .clone()
        .map(|i| free.whole.field().charges()[i])
        .collect();
    println!(
        "\nFreeSolv's benzene: {}, SHA-256 {}\n\n| | bound: element | QEq | AM1-BCC | in water: \
         element | QEq | AM1-BCC |\n| --- | --- | --- | --- | --- | --- | --- |",
        freesolv::SOURCE,
        freesolv::SHA256
    );
    for k in 0..qb.charges.len() {
        println!(
            "| {k} | {:?} | {:+.4} | {:+.4} | {:?} | {:+.4} | {:+.4} |",
            graph_b.0[k], qeq_b[k], qb.charges[k], graph_w.0[k], qeq_w[k], qw.charges[k]
        );
    }
    let sum = |q: &[f64]| q.iter().sum::<f64>();
    println!(
        "\nΣq: QEq {:+.2e}, {:+.2e}; AM1-BCC {:+.2e}, {:+.2e}; spread over {} automorphisms {:.1e}",
        sum(&qeq_b),
        sum(&qeq_w),
        sum(&qb.charges),
        sum(&qw.charges),
        qb.automorphisms,
        qb.spread
    );

    let rb = reparametrise(&complex, &complex_frames, 20);
    let rw = reparametrise(&water, &water_frames, 20);
    let sb = swap_frames(
        &bound,
        &complex_frames,
        [&qb.charges, &nb],
        20,
        "complex",
        &frames_dir,
    );
    let sw = swap_frames(
        &free,
        &water_frames,
        [&qw.charges, &nw],
        20,
        "water",
        &frames_dir,
    );
    log(
        &frames_dir,
        "charges",
        &format!("swapped, {:.0} s", t.elapsed().as_secs_f64()),
    );

    let v3 = |r: &Reparametrised| -> Vec<f64> { r.du.iter().map(|d| d[3]).collect() };
    let joint =
        |a: &[f64], b: &[f64]| -> Vec<f64> { a.iter().zip(b).map(|(x, y)| x + y).collect() };
    let (v3b, v3w) = (v3(&rb), v3(&rw));
    let summaries = [
        report_swap(
            "3: GAFF benzene, ff14SB protein, Lorentz–Berthelot (van der Waals only, 0ff5194's)",
            kt,
            &v3b,
            &v3w,
        ),
        report_swap(
            "e: QEq → FreeSolv's AM1-BCC charges, UFF's van der Waals",
            kt,
            &sb.elec,
            &sw.elec,
        ),
        report_swap(
            "e′: QEq → AM1-BCC neutral and symmetric (±0.1301)",
            kt,
            &sb.neutral,
            &sw.neutral,
        ),
        report_swap(
            "4: variant 3's van der Waals and FreeSolv's AM1-BCC charges, jointly",
            kt,
            &joint(&v3b, &sb.elec),
            &joint(&v3w, &sw.elec),
        ),
    ];

    println!(
        "\n**e by class, bound**, first order:\n\n| class | ⟨ΔU_elec,c⟩ (kcal/mol) |\n| --- | --- |"
    );
    for (c, class) in CLASSES.iter().enumerate() {
        if bound.classes[c].is_some() {
            println!(
                "| {class} | {} |",
                pm(&estimate_kcal(sb.classes.iter().map(|x| x[c])))
            );
        }
    }
    println!(
        "| in water, all water | {} |",
        pm(&estimate_kcal(sw.elec.iter().copied()))
    );
    println!(
        "\nbenzene's own Coulomb energy's change, QEq → AM1-BCC, left out of every ΔU: bound {}, \
         in water {} kcal/mol, every twentieth frame",
        pm(&estimate_kcal(sb.own.iter().copied())),
        pm(&estimate_kcal(sw.own.iter().copied()))
    );

    let (hyd, _) = w3_hydration();
    println!(
        "\n| perturbation | ΔΔG_bind: first order | second cumulant | Zwanzig | halves | σ(βΔU) \
         bound, water | effective samples bound, water | implied ΔG_hyd | ΔG°_bind | σ ≤ \
         {ONE_STEP_SPREAD} |\n| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
    );
    for x in &summaries {
        println!(
            "| {} | {:+.3} | {:+.3} | {:+.3} ± {:.3} | {:+.3}, {:+.3} | {:.2}, {:.2} | {:.0}, {:.0} \
             | {:+.2} | {:+.2} | {} |",
            x.name,
            x.first,
            x.cumulant,
            x.exp,
            x.error,
            x.halves[0],
            x.halves[1],
            x.spread[0],
            x.spread[1],
            x.ess[0],
            x.ess[1],
            hyd + x.water.0,
            W4_BINDING.0 + x.exp,
            if x.spread[0] <= ONE_STEP_SPREAD && x.spread[1] <= ONE_STEP_SPREAD {
                "yes"
            } else {
                "**no**"
            }
        );
    }
    println!(
        "\nagainst: ΔG_hyd experiment {:+.2} ± {:.2}, FreeSolv's GAFF/AM1-BCC {:+.2} ± {:.2}, W3 \
         {hyd:+.2}; ΔG°_bind experiment {:+.2} ± {:.2}, Mobley et al. 2007 {:+.2} ± {:.2}, W4 \
         {:+.3}",
        freesolv::EXPERIMENT.0,
        freesolv::EXPERIMENT.1,
        freesolv::CALCULATED.0,
        freesolv::CALCULATED.1,
        EXPERIMENT.0,
        EXPERIMENT.1,
        MOBLEY_2007.0,
        MOBLEY_2007.1,
        W4_BINDING.0
    );
    println!(
        "\nthe identities, worst over the frames (1 is the bound): the classes' sum {:.3} bound, \
         {:.3} in water; two whole evaluations {:.3}, {:.3}; variant 3's ΔU against its sums {:.3}, \
         {:.3}; the field's coupling {:.3}, {:.3}. {:.0} s in all",
        sb.against_classes,
        sw.against_classes,
        sb.against_evaluations,
        sw.against_evaluations,
        rb.against_sums,
        rw.against_sums,
        rb.against_field,
        rw.against_field,
        t.elapsed().as_secs_f64()
    );
    log(
        &frames_dir,
        "charges",
        &format!("done in {:.0} s", t.elapsed().as_secs_f64()),
    );
    for x in [&sb, &sw] {
        assert!(x.against_classes <= 1.0, "{}", x.against_classes);
        assert!(x.against_evaluations <= 1.0, "{}", x.against_evaluations);
    }
    for r in [&rb, &rw] {
        assert!(r.against_sums <= 1.0, "{}", r.against_sums);
        assert!(r.parts <= 1.0, "{}", r.parts);
        assert!(r.against_field <= 1.0, "{}", r.against_field);
        assert!(r.against_hybrid <= 1.0, "{}", r.against_hybrid);
    }
}

// ---------------------------------------------------------------------------------------------
// Variant 4's one step, checked by staged sampling: H(s) = H_UFF + s ΔU₄

/// The states of `H(s) = H_UFF + s ΔU₄` the staged check samples, in path order.
const STAGES: [f64; 5] = [0.0, 0.25, 0.5, 0.75, 1.0];

/// Each state's index for [`window_seed`]. **`s` = 0 is index 0**: under [`PROTOCOL`] that is the
/// leg's window 0 — the same box, start, state, dynamics and seed — and under [`FRAMES`] the saved
/// frames', so that state's trajectory is theirs to the bit if nothing it runs on has changed
/// since. The others are past [`grid`]'s 33 candidates, seeded apart from every window the leg has
/// run or could insert.
const STAGE_INDEX: [usize; 5] = [0, 33, 34, 35, 36];

/// How often a staged state holds its dynamics' `H(s) − H(0)` to the frames' `s ΔU₄`, in samples:
/// the first sample and every 200th after it.
const STAGE_CHECK_EVERY: usize = 200;

/// Variant 4's one-step Zwanzig `ΔΔG_bind` on the saved frames and its error, kcal/mol
/// (`benzenes_charges_perturbed_measured`, run once on 2026-10-19).
const VARIANT_4_ONE_STEP: (f64, f64) = (3.953, 0.113);

/// What that one step implied, kcal/mol: `ΔG°_bind` and `ΔG_hyd`.
const VARIANT_4_IMPLIED: [f64; 2] = [-5.02, -0.62];

/// A Lennard-Jones pair's `dE/dr` at `r` for `[x, D]`, as the crate's `Pair` writes it:
/// `12 D (s⁶ − s¹²)/r` with `s = x/r`.
fn lj_slope(p: [f64; 2], r: f64) -> f64 {
    let s = p[0] / r;
    let s6 = (s * s) * (s * s) * (s * s);
    let s12 = s6 * s6;
    12.0 * p[1] * (s6 - s12) / r
}

/// **The van der Waals half of `ΔU₄`, with its forces**: benzene's cross van der Waals from the
/// current model to variant 3 ([`variants`]' fourth: GAFF's benzene, ff14SB's protein,
/// Lorentz–Berthelot). Every pair of a benzene atom with an atom of the rest inside the cutoff,
/// found by the same wrapped minimum image as [`cross_deltas`], gives `u₃(r) − u(r)` and its
/// derivative; an atom with no well in either model gives exactly zero and is left out. **The
/// long-range corrections' difference** is a constant: it depends on the box's volume and the
/// parameters and on no position, so at a fixed volume it has no force, and none is applied.
struct VdwSwap<'a> {
    sur: &'a Surroundings,
    from: CrossModel,
    to: CrossModel,
    /// The atoms of the rest with a well in either model.
    partners: Vec<usize>,
    /// The long-range corrections' difference, class by class as [`cross_deltas`] takes it, with
    /// its `Σ|term|` and terms.
    tail: CrossSum,
}

impl<'a> VdwSwap<'a> {
    fn new(sur: &'a Surroundings, labels: &Labels) -> VdwSwap<'a> {
        let from = current_model(sur);
        let to = variants(sur, labels).swap_remove(3);
        let partners = (0..sur.kinds.len())
            .filter(|&j| sur.kinds[j].is_some() && (from.vdw[j][1] != 0.0 || to.vdw[j][1] != 0.0))
            .collect();
        let (a, b) = (
            class_tails(sur, labels, &from),
            class_tails(sur, labels, &to),
        );
        let mut tail = CrossSum::default();
        for c in 0..a.len() {
            tail.tail += b[c].tail - a[c].tail;
            tail.magnitude += a[c].magnitude + b[c].magnitude;
            tail.terms += b[c].terms;
        }
        VdwSwap {
            sur,
            from,
            to,
            partners,
            tail,
        }
    }

    /// `ΔU_vdW` at `at`, its pairs and its long-range term apart; with `forces`, `scale` times
    /// its force added into them.
    fn evaluate(
        &self,
        at: &[[f64; 3]],
        scale: f64,
        mut forces: Option<&mut [[f64; 3]]>,
    ) -> CrossSum {
        let sur = self.sur;
        let w: Vec<[f64; 3]> = at.iter().map(|p| sur.cell.wrap(*p)).collect();
        let rc2 = sur.cutoff * sur.cutoff;
        let mut sum = CrossSum::default();
        for g in sur.ligand.clone() {
            for &j in &self.partners {
                let (a, b) = (g.min(j), g.max(j));
                let (d, r2) = sur.separation(&w, a, b);
                if r2 >= rc2 {
                    continue;
                }
                let r = r2.sqrt();
                let p0 = self.from.mixing.pair(self.from.vdw[a], self.from.vdw[b]);
                let p1 = self.to.mixing.pair(self.to.vdw[a], self.to.vdw[b]);
                let energy = |p: [f64; 2]| {
                    Pair {
                        atoms: [a, b],
                        distance: p[0],
                        well: p[1],
                    }
                    .energy(r)
                };
                let (u0, u1) = (energy(p0), energy(p1));
                sum.pairs += u1 - u0;
                sum.magnitude += u0.abs() + u1.abs();
                sum.terms += 1;
                if let Some(f) = forces.as_deref_mut() {
                    // −(dΔu/dr)/r along d = r_a − r_b, on a; the opposite on b.
                    let pull = -scale * (lj_slope(p1, r) - lj_slope(p0, r)) / r;
                    for x in 0..3 {
                        f[a][x] += pull * d[x];
                        f[b][x] -= pull * d[x];
                    }
                }
            }
        }
        sum.tail = self.tail.tail;
        sum.magnitude += self.tail.magnitude;
        sum.terms += self.tail.terms;
        sum
    }
}

/// **Benzene's own Coulomb pairs taken back to `old`'s charges.** The field computes them in
/// vacuum, on whole molecules, with its own charges at every λ; `H(s)` keeps them QEq's, as the
/// one step did. So each own pair's `k_e q_i q_j / r` under `charged`'s charges less under `old`'s
/// is returned, to be subtracted, and with `forces` its force is subtracted from them. With the
/// pairs' `Σ|term|`.
fn own_back(
    charged: &PeriodicDecoupling,
    old: &[f64],
    at: &[[f64; 3]],
    mut forces: Option<&mut [[f64; 3]]>,
) -> (f64, f64) {
    let q = charged.field().charges();
    let whole = charged.field().whole(at);
    let (mut back, mut size) = (0.0, 0.0);
    for p in charged.own_pairs() {
        let [i, j] = p.atoms;
        let d = sub(whole[i], whole[j]);
        let r = len(d);
        let (now, then) = (
            pantometry_forcefield::energy::coulomb(q[i], q[j], r),
            pantometry_forcefield::energy::coulomb(old[i], old[j], r),
        );
        let delta = now - then;
        back += delta;
        size += now.abs() + then.abs();
        if let Some(f) = forces.as_deref_mut() {
            // The field put k_e q_i q_j d/r³ on i and its opposite on j.
            let push = delta / (r * r);
            for x in 0..3 {
                f[i][x] -= push * d[x];
                f[j][x] += push * d[x];
            }
        }
    }
    (back, size)
}

/// **`H(s) = H_UFF + s ΔU₄`**, the Hamiltonian a staged state samples, with `ΔU₄` variant 4's:
/// benzene's cross van der Waals made variant 3's and its cross charges FreeSolv's AM1-BCC, its
/// own pairs left as they are.
///
/// **One whole evaluation a step, the base's cost.** Benzene's cross Coulomb energy `C(q)` — its
/// real-space pairs across, the reciprocal sum's cross terms and the background's — is linear in
/// benzene's charges, so `C((1 − s) q_old + s q_new) = C(q_old) + s [C(q_new) − C(q_old)]`
/// exactly, real space, mesh or classical sum and background alike, and the rest's own terms do
/// not see benzene's charges. So the field is evaluated once with benzene's charges
/// `(1 − s) q_old + s q_new` — every force of `s ΔU_elec`, the reciprocal ones on atoms far from
/// benzene included, comes out of the field's own evaluation — and two things that are not
/// linear, or not the field's, are put right by hand: **benzene's own Coulomb pairs**, quadratic
/// in its charges, taken back to QEq's ([`own_back`], 36 pairs), and **`s ΔU_vdW`**, added
/// ([`VdwSwap`], 12 × N pairs). At `s` = 0 the charges are QEq's to the bit, so `H(0)` is the
/// field's energy and forces to the bit; at `s` = 1 they are FreeSolv's to the bit.
struct Staged<'a> {
    s: f64,
    lambda: Lambda,
    /// The field `H(s)` is a change of, at `s` = 0.
    base: &'a PeriodicDecoupling,
    /// The same with benzene's charges `(1 − s) q_old + s q_new`.
    charged: PeriodicDecoupling,
    vdw: &'a VdwSwap<'a>,
}

impl<'a> Staged<'a> {
    /// `H(s)` on `base`, whose restraint is `restraint` (`base` gives no way to read it back, so
    /// it is given again; [`the_staged_hamiltonian_at_zero_is_the_field`] holds the two equal),
    /// with benzene's new charges `to` in the decoupled group's order, at `lambda`.
    fn new(
        base: &'a PeriodicDecoupling,
        restraint: Option<Boresch>,
        to: &[f64],
        vdw: &'a VdwSwap<'a>,
        lambda: Lambda,
        s: f64,
    ) -> Staged<'a> {
        let group = base.group();
        let members: Vec<usize> = (0..group.len()).filter(|&i| group[i]).collect();
        assert_eq!(members.len(), to.len(), "one charge per benzene atom");
        let mut q = base.field().charges().to_vec();
        for (k, &i) in members.iter().enumerate() {
            q[i] = (1.0 - s) * q[i] + s * to[k];
        }
        let mut charged = PeriodicDecoupling::new(&base.field().clone().with_charges(q), group)
            .expect("benzene decouples");
        if let Some(r) = restraint {
            charged = charged.with_restraint(r);
        }
        Staged {
            s,
            lambda,
            base,
            charged,
            vdw,
        }
    }

    /// `H(s)` at `at`, joules, its forces into `forces`.
    fn evaluate(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let e = Alchemical::energy_and_forces(&self.charged, at, self.lambda, forces);
        let (own, _) = own_back(
            &self.charged,
            self.base.field().charges(),
            at,
            Some(&mut *forces),
        );
        let v = self.vdw.evaluate(at, self.s, Some(&mut *forces));
        e - own + self.s * v.total()
    }
}

impl pantometry_forcefield::Potential for Staged<'_> {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        self.evaluate(at, forces)
    }
}

/// One configuration's `ΔU₄` **as the frames compute it**, joules: variant 3's `ΔU_vdW` pair by
/// pair and `ΔU_elec` through the couplings, with the van der Waals sum's `Σ|term|` and terms.
#[derive(Clone, Copy, Debug)]
struct FrameChange {
    vdw: f64,
    elec: f64,
    magnitude: f64,
    terms: usize,
}

impl FrameChange {
    /// `ΔU₄`, as `benzenes_charges_perturbed_measured` adds it: the van der Waals first.
    fn total(&self) -> f64 {
        self.vdw + self.elec
    }
}

/// **`ΔU₄` by the frames' own computation**, not the dynamics': [`cross_deltas`] from the current
/// model to variant 3 with [`class_tails`]' corrections, which is [`reparametrise`]'s fourth `ΔU`
/// to the bit, and [`ChargeSwap::delta`] for FreeSolv's charges, `C(q_new) − C(q_old)` from two
/// passes of the couplings. What each staged state records and what BAR is run on.
struct FramesDelta<'a> {
    sur: &'a Surroundings,
    labels: &'a Labels,
    from: CrossModel,
    to: CrossModel,
    from_tails: Vec<CrossSum>,
    to_tails: Vec<CrossSum>,
    swap: ChargeSwap,
}

impl<'a> FramesDelta<'a> {
    fn new(sur: &'a Surroundings, labels: &'a Labels, charges: &[f64]) -> FramesDelta<'a> {
        let from = current_model(sur);
        let to = variants(sur, labels).swap_remove(3);
        FramesDelta {
            from_tails: class_tails(sur, labels, &from),
            to_tails: class_tails(sur, labels, &to),
            from,
            to,
            swap: ChargeSwap::new(&sur.whole, charges),
            sur,
            labels,
        }
    }

    fn at(&self, at: &[[f64; 3]]) -> FrameChange {
        let d = cross_deltas(
            self.sur,
            self.labels,
            at,
            (&self.from, &self.from_tails),
            &[(&self.to, &self.to_tails)],
        )
        .swap_remove(0);
        FrameChange {
            vdw: d.total(),
            elec: self.swap.delta(at),
            magnitude: d.magnitude,
            terms: d.terms,
        }
    }
}

/// **The dynamics' `ΔU₄` held to the frames'**: `|H(s) − H(0) − s ΔU₄|` at `at`, `H(0)` the base
/// field's own evaluation and `ΔU₄` `change`, over `64 ε` of both evaluations' parts and the own
/// pairs' — W1's and W3's bound for the cross terms through the mesh — plus the summation bound of
/// the two van der Waals sums. A ratio above 1 is a failure; `NaN` is one too.
fn staged_held(st: &Staged, at: &[[f64; 3]], change: &FrameChange) -> f64 {
    let mut f = vec![[0.0; 3]; at.len()];
    let hs = st.evaluate(at, &mut f);
    let h0 = Alchemical::energy_and_forces(st.base, at, st.lambda, &mut f);
    let (_, own) = own_back(&st.charged, st.base.field().charges(), at, None);
    let v = st.vdw.evaluate(at, 0.0, None);
    let parts = magnitude(&st.charged.energy(at, st.lambda))
        + magnitude(&st.base.energy(at, st.lambda))
        + own;
    let bound = 64.0 * EPS * parts
        + summation_bound(v.terms + change.terms, v.magnitude + change.magnitude);
    ((hs - h0) - st.s * change.total()).abs() / bound
}

/// One leg of the staged check: the field `H(s)` changes and its restraint, FreeSolv's charges in
/// benzene's order, `ΔU_vdW` with forces and `ΔU₄` as the frames compute it, the dynamics and the
/// state, what says what `start` is (in the header), the start — called only when a state has to
/// run — and the log's name.
struct StagedPhase<'a> {
    name: &'static str,
    base: &'a PeriodicDecoupling,
    restraint: Option<Boresch>,
    charges: &'a [f64],
    vdw: &'a VdwSwap<'a>,
    frames: &'a FramesDelta<'a>,
    template: &'a MolecularDynamics,
    lambda: Lambda,
    identity: String,
    start: &'a dyn Fn() -> Vec<[f64; 3]>,
    log: &'static str,
}

/// **The potential state `k` runs on**: `H(STAGES[k])` on the phase. The one place a state's `s`
/// is turned into a Hamiltonian; [`the_staged_states_are_what_their_headers_say`] holds it to the
/// header's `s`.
fn stage_potential<'a>(phase: &StagedPhase<'a>, k: usize) -> Staged<'a> {
    Staged::new(
        phase.base,
        phase.restraint,
        phase.charges,
        phase.vdw,
        phase.lambda,
        STAGES[k],
    )
}

/// State `k`'s header: the phase, the state and its `s` among all of them, its seed's index and
/// the seed, the protocol literally, the phase's λ, the atom count, the phase's identity, and
/// digests of benzene's two sets of charges and of both van der Waals models, so that a state run
/// differently is refused rather than mixed with these.
fn staged_header(phase: &StagedPhase, p: Protocol, k: usize) -> String {
    let l = phase.lambda;
    let group = phase.base.group();
    let qeq = (0..group.len())
        .filter(|&i| group[i])
        .map(|i| phase.base.field().charges()[i]);
    let rule = |m: Mixing| f64::from(u8::from(m == Mixing::LorentzBerthelot));
    let models = digest(
        phase
            .vdw
            .from
            .vdw
            .iter()
            .chain(&phase.vdw.to.vdw)
            .flatten()
            .copied()
            .chain([rule(phase.vdw.from.mixing), rule(phase.vdw.to.mixing)]),
    );
    format!(
        "# staged {} state {k} s {} of {} index {} seed {:016x} protocol {} {} {} {} {} {} {:#x} \
         lambda {} {} {} atoms {} {} charges {:016x} {:016x} vdw {models:016x}",
        phase.name,
        exact(STAGES[k]),
        STAGES.map(exact).join(","),
        STAGE_INDEX[k],
        window_seed(p.seed, STAGE_INDEX[k]),
        exact(p.time_step),
        exact(p.temperature),
        exact(p.friction),
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        exact(l.restraint),
        exact(l.electrostatics),
        exact(l.van_der_waals),
        group.len(),
        phase.identity,
        digest(qeq),
        digest(phase.charges.iter().copied())
    )
}

/// The state, `s`, seed index and seed a staged header says; `None` for a line that does not say
/// them.
fn staged_state_of(header: &str) -> Option<(usize, f64, usize, u64)> {
    let w: Vec<&str> = header.split_whitespace().collect();
    let field = |name: &str| w.get(w.iter().position(|&x| x == name)? + 1).copied();
    let state: usize = field("state")?.parse().ok()?;
    (state < STAGES.len()).then_some(())?;
    Some((
        state,
        field("s")?.parse().ok()?,
        field("index")?.parse().ok()?,
        u64::from_str_radix(field("seed")?, 16).ok()?,
    ))
}

/// One staged state, read back from its file: per sample its step, `ΔU₄` and its two halves
/// ([`FramesDelta`]), the coupling at the phase's state — the number the leg's window 0 recorded
/// at its own state — and a digest of the positions; the worst [`staged_held`] over the checked
/// samples, and the seconds it took.
#[derive(Clone, Debug)]
struct StagedRecord {
    state: usize,
    s: f64,
    steps: Vec<u64>,
    du: Vec<f64>,
    vdw: Vec<f64>,
    elec: Vec<f64>,
    coupling: Vec<f64>,
    digests: Vec<u64>,
    held: f64,
    seconds: f64,
}

impl StagedRecord {
    fn empty(state: usize) -> StagedRecord {
        StagedRecord {
            state,
            s: STAGES[state],
            steps: Vec::new(),
            du: Vec::new(),
            vdw: Vec::new(),
            elec: Vec::new(),
            coupling: Vec::new(),
            digests: Vec::new(),
            held: 0.0,
            seconds: f64::NAN,
        }
    }

    /// Every number it holds, as floats: what its trailer's digest is of.
    fn numbers(&self) -> Vec<f64> {
        let mut v: Vec<f64> = self.steps.iter().map(|&s| s as f64).collect();
        for xs in [&self.du, &self.vdw, &self.elec, &self.coupling] {
            v.extend_from_slice(xs);
        }
        v.extend(self.digests.iter().map(|&d| f64::from_bits(d)));
        v.push(self.held);
        v.push(self.seconds);
        v
    }

    /// The same record with only the samples in `range`.
    fn slice(&self, range: std::ops::Range<usize>) -> StagedRecord {
        StagedRecord {
            steps: self.steps[range.clone()].to_vec(),
            du: self.du[range.clone()].to_vec(),
            vdw: self.vdw[range.clone()].to_vec(),
            elec: self.elec[range.clone()].to_vec(),
            coupling: self.coupling[range.clone()].to_vec(),
            digests: self.digests[range].to_vec(),
            ..self.clone()
        }
    }
}

/// A staged file's last line: the sample count, the seconds, the worst identity and a digest of
/// every number in the file, the digest last and sixteen digits wide, so that a file cut anywhere
/// is refused.
fn staged_trailer(r: &StagedRecord) -> String {
    format!(
        "# done {} seconds {} held {} digest {:016x}",
        r.du.len(),
        exact(r.seconds),
        exact(r.held),
        digest(r.numbers())
    )
}

/// A staged file's reading: its header line's state, then each sample line — the step, `ΔU₄`,
/// `ΔU_vdW`, `ΔU_elec`, the coupling, the positions' digest — then the trailer, which must hold
/// the count and the digest of what was read. Refused, and why, when it is not whole.
fn parse_staged(text: &str) -> Result<StagedRecord, String> {
    let mut lines = text.lines();
    let head = lines.next().ok_or_else(|| "an empty file".to_string())?;
    let (state, _, _, _) =
        staged_state_of(head).ok_or_else(|| format!("a header without its state: {head:?}"))?;
    let mut r = StagedRecord::empty(state);
    for line in lines {
        if let Some(rest) = line.strip_prefix("# done ") {
            let w: Vec<&str> = rest.split_whitespace().collect();
            let [count, "seconds", secs, "held", held, "digest", d] = w[..] else {
                return Err(format!("a trailer cut short: {line:?}"));
            };
            r.seconds = secs.parse().map_err(|_| "the seconds".to_string())?;
            r.held = held.parse().map_err(|_| "the identity".to_string())?;
            if count.parse::<usize>().ok() != Some(r.du.len()) {
                return Err(format!(
                    "{} samples where the trailer says {count}",
                    r.du.len()
                ));
            }
            if d.len() != 16 || u64::from_str_radix(d, 16).ok() != Some(digest(r.numbers())) {
                return Err("the samples' digest is not the trailer's".into());
            }
            return Ok(r);
        }
        let w: Vec<&str> = line.split_whitespace().collect();
        let sample = || -> Option<(u64, [f64; 4], u64)> {
            let [step, du, vdw, elec, coupling, d] = w[..] else {
                return None;
            };
            Some((
                step.parse().ok()?,
                [
                    du.parse().ok()?,
                    vdw.parse().ok()?,
                    elec.parse().ok()?,
                    coupling.parse().ok()?,
                ],
                u64::from_str_radix(d, 16).ok()?,
            ))
        };
        let (step, x, d) = sample().ok_or_else(|| format!("not a sample: {line:?}"))?;
        r.steps.push(step);
        r.du.push(x[0]);
        r.vdw.push(x[1]);
        r.elec.push(x[2]);
        r.coupling.push(x[3]);
        r.digests.push(d);
    }
    Err("no trailer: the file is not whole".into())
}

/// **State `k` of a phase: read back from `NAME_TAG_K.txt` in `dir` when a run completed it under
/// the same header** ([`staged_header`]), refused by name when it was run differently or is not
/// whole; or run now — [`sample_on`] with [`stage_potential`]'s `H(s)` from the phase's start at
/// [`STAGE_INDEX`]`[k]`'s seed, every sample's `ΔU₄` by [`FramesDelta`], streamed to
/// `NAME_TAG_K.partial`, which is begun afresh, and renamed once whole, then read back from the
/// file. Every [`STAGE_CHECK_EVERY`]th sample is held to the frames' `s ΔU₄` ([`staged_held`]) on
/// a second `H(s)`, so that the check's evaluations leave the dynamics' own field untouched.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn staged_state(
    dir: &std::path::Path,
    tag: &str,
    phase: &StagedPhase,
    p: Protocol,
    k: usize,
) -> StagedRecord {
    use std::io::Write;
    let path = dir.join(format!("{}_{tag}_{k}.txt", phase.name));
    let want = staged_header(phase, p, k);
    if let Ok(text) = std::fs::read_to_string(&path) {
        assert_eq!(
            text.lines().next().unwrap_or(""),
            want,
            "{} was run differently: remove it to run the state again",
            path.display()
        );
        let r = parse_staged(&text).unwrap_or_else(|why| {
            panic!(
                "{} is not whole ({why}): remove it to run the state again",
                path.display()
            )
        });
        assert_eq!(r.du.len(), p.samples, "{}", path.display());
        log(
            dir,
            phase.log,
            &format!("{} state {k}: {} samples read back", phase.name, p.samples),
        );
        return r;
    }
    let start = (phase.start)();
    assert_eq!(start.len(), phase.base.group().len());
    let potential = stage_potential(phase, k);
    let checker = stage_potential(phase, k);
    let t = std::time::Instant::now();
    let partial = path.with_extension("partial");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
    writeln!(file, "{want}").unwrap();
    log(
        dir,
        phase.log,
        &format!(
            "{} state {k} at s = {} begins: {} steps",
            phase.name,
            STAGES[k],
            p.steps_per_window()
        ),
    );
    let mut r = StagedRecord::empty(k);
    let mut last = (0u64, 0.0f64);
    sample_on(
        &potential,
        &start,
        phase.template,
        p,
        STAGE_INDEX[k],
        |step, at| {
            let change = phase.frames.at(at);
            let coupling = phase.base.couplings(at, &[phase.lambda])[0].energy;
            let d = positions_digest(at);
            if r.du.len() % STAGE_CHECK_EVERY == 0 {
                let held = staged_held(&checker, at, &change);
                if held.is_nan() || held > r.held {
                    r.held = held;
                }
            }
            writeln!(
                file,
                "{step} {} {} {} {} {d:016x}",
                exact(change.total()),
                exact(change.vdw),
                exact(change.elec),
                exact(coupling)
            )
            .unwrap();
            r.steps.push(step);
            r.du.push(change.total());
            r.vdw.push(change.vdw);
            r.elec.push(change.elec);
            r.coupling.push(coupling);
            r.digests.push(d);
            if r.du.len() % 100 == 0 {
                file.flush().unwrap();
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    phase.log,
                    &format!(
                        "{} state {k}: {} of {} samples, step {step}, {secs:.0} s, {:.1} ms a step \
                         since the last report",
                        phase.name,
                        r.du.len(),
                        p.samples,
                        (secs - last.1) * 1e3 / (step - last.0) as f64
                    ),
                );
                last = (step, secs);
            }
        },
    );
    r.seconds = t.elapsed().as_secs_f64();
    writeln!(file, "{}", staged_trailer(&r)).unwrap();
    file.flush().unwrap();
    drop(file);
    std::fs::rename(&partial, &path).expect("rename");
    log(
        dir,
        phase.log,
        &format!(
            "{} state {k} done in {:.0} s, the identity's worst {:.3} of its bound",
            phase.name, r.seconds, r.held
        ),
    );
    parse_staged(&std::fs::read_to_string(&path).unwrap()).expect("just written")
}

/// Every state of a phase, in path order ([`staged_state`]).
fn staged_states(
    dir: &std::path::Path,
    tag: &str,
    phase: &StagedPhase,
    p: Protocol,
) -> Vec<StagedRecord> {
    (0..STAGES.len())
        .map(|k| staged_state(dir, tag, phase, p, k))
        .collect()
}

/// BAR between two neighbouring staged states, reduced: `u_b − u_a = (s_b − s_a) βΔU₄` on both
/// states' samples, the delta method's variance from each side's autocorrelation.
fn staged_interval(a: &StagedRecord, b: &StagedRecord, kt: f64) -> pantometry_forcefield::Bar {
    let ds = b.s - a.s;
    let f: Vec<f64> = a.du.iter().map(|u| ds * u / kt).collect();
    let r: Vec<f64> = b.du.iter().map(|u| ds * u / kt).collect();
    pantometry_forcefield::free_energy::bar_correlated(&f, &r)
}

/// BAR along the states, kcal/mol, with its error: [`bennett_chain`], window by window, as the leg
/// is chained.
fn staged_chain(records: &[StagedRecord], kt: f64) -> (f64, f64) {
    let n = records.len();
    let windows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|k| {
            let towards = |j: usize| -> Vec<f64> {
                let ds = records[j].s - records[k].s;
                records[k].du.iter().map(|u| ds * u / kt).collect()
            };
            (
                if k > 0 { towards(k - 1) } else { Vec::new() },
                if k + 1 < n {
                    towards(k + 1)
                } else {
                    Vec::new()
                },
            )
        })
        .collect();
    let fe = bennett_chain(&windows);
    (kcal(fe.value * kt), kcal(fe.error * kt))
}

/// TI on the trapezoid, kcal/mol: `∫₀¹ ⟨ΔU₄⟩_s ds`, `dH/ds` being `ΔU₄` exactly, with
/// [`Estimate`]'s errors.
fn staged_ti(records: &[StagedRecord]) -> (f64, f64) {
    let n = records.len();
    let (mut value, mut variance) = (0.0, 0.0);
    for k in 0..n {
        let lo = records[k.saturating_sub(1)].s;
        let hi = records[(k + 1).min(n - 1)].s;
        let w = 0.5 * (hi - lo);
        let e = estimate_kcal(records[k].du.iter().copied());
        value += w * e.mean;
        variance += w * w * e.error * e.error;
    }
    (value, variance.sqrt())
}

/// What one phase's staged states give, kcal/mol.
struct StagedSummary {
    /// BAR along the chain and its error.
    bar: (f64, f64),
    /// The same on each half of every state's samples.
    halves: [f64; 2],
    ti: (f64, f64),
    /// Zwanzig's one step from the fresh `s` = 0 state forward, and from `s` = 1 backward.
    one_step: [f64; 2],
    /// The lowest overlap of the four intervals.
    overlap: f64,
}

/// The tables for one phase: each state's `⟨ΔU₄⟩` and its halves, σ(βΔU₄) and the seconds; each
/// interval's BAR and overlap; and the chain, its halves, TI, and the one step from either end.
fn report_staged(name: &str, records: &[StagedRecord], kt: f64) -> StagedSummary {
    println!(
        "\n**{name}: the states**\n\n| s | samples | ⟨ΔU₄⟩ | ⟨ΔU_vdW⟩ | ⟨ΔU_elec⟩ | σ(βΔU₄) | \
         identity (of its bound) | seconds |\n| --- | --- | --- | --- | --- | --- | --- | --- |"
    );
    for r in records {
        println!(
            "| {:.2} | {} | {} | {} | {} | {:.3} | {:.3} | {:.0} |",
            r.s,
            r.du.len(),
            pm(&estimate_kcal(r.du.iter().copied())),
            pm(&estimate_kcal(r.vdw.iter().copied())),
            pm(&estimate_kcal(r.elec.iter().copied())),
            zwanzig(&r.du, kt).spread,
            r.held,
            r.seconds
        );
    }
    println!("\n| interval | BAR (kcal/mol) | overlap |\n| --- | --- | --- |");
    let mut overlap = f64::INFINITY;
    for k in 0..records.len() - 1 {
        let b = staged_interval(&records[k], &records[k + 1], kt);
        overlap = overlap.min(b.overlap);
        println!(
            "| {:.2} → {:.2} | {:+.3} ± {:.3} | {:.3} |",
            records[k].s,
            records[k + 1].s,
            kcal(b.delta * kt),
            kcal(b.error() * kt),
            b.overlap
        );
    }
    let bar = staged_chain(records, kt);
    let half = records[0].du.len() / 2;
    let first: Vec<StagedRecord> = records.iter().map(|r| r.slice(0..half)).collect();
    let second: Vec<StagedRecord> = records.iter().map(|r| r.slice(half..r.du.len())).collect();
    let halves = [staged_chain(&first, kt).0, staged_chain(&second, kt).0];
    let ti = staged_ti(records);
    let last = records.last().expect("a state");
    let back: Vec<f64> = last.du.iter().map(|u| -u).collect();
    let one_step = [zwanzig(&records[0].du, kt).exp, -zwanzig(&back, kt).exp];
    println!(
        "\n{name}: ΔG_swap by BAR along the chain {:+.3} ± {:.3} kcal/mol (halves {:+.3}, {:+.3}); \
         TI {:+.3} ± {:.3}; Zwanzig's one step from s = 0 {:+.3}, from s = 1 back {:+.3}; the lowest \
         overlap {overlap:.3}",
        bar.0, bar.1, halves[0], halves[1], ti.0, ti.1, one_step[0], one_step[1]
    );
    StagedSummary {
        bar,
        halves,
        ti,
        one_step,
        overlap,
    }
}

/// **How many of a state's numbers are another run's to the bit**: equal, and compared.
fn to_the_bit(a: &[u64], b: &[u64]) -> (usize, usize) {
    let n = a.len().min(b.len());
    ((0..n).filter(|&i| a[i] == b[i]).count(), n)
}

/// **The staged check of variant 4, at protocol `p`, in `staged/` under the results directory.**
/// The water first, W3's box melted as W3 melts it, and then the complex, the leg's box and start;
/// each through [`staged_states`], every state resumable on its own, under the lock and the log
/// `TAG`. Then BAR along each chain, `ΔΔG_bind = ΔG_swap(bound) − ΔG_swap(water)`, and what it
/// implies, against the one step on the saved frames and against experiment. And `s` = 0 against
/// what it reproduces if nothing it runs on has changed: under [`PROTOCOL`] the leg's window 0
/// (the coupling at its state, sample by sample), under [`FRAMES`] the saved frames (the
/// positions' digests); reported, not asserted. Every state's identity is asserted after
/// everything is printed. **Nothing outside `staged/` is written**: `prepared.txt` is read back
/// through its header ([`saved_start`]).
fn staged_measured(tag: &'static str, p: Protocol) {
    let dir = results_dir();
    assert!(
        dir.join("prepared.txt").exists(),
        "{} is not there: this measurement reads the leg's preparation and writes nothing beside \
         it",
        dir.join("prepared.txt").display()
    );
    let staged = dir.join("staged");
    std::fs::create_dir_all(&staged).expect("the staged directory");
    let _lock = Lock::take(&staged, tag);
    let t = std::time::Instant::now();
    let kt = p.kt();

    let benzene = hydration::Benzene::new(Some(pocket3().ligand_force_field().charges().to_vec()));
    assert_eq!(
        &benzene.types[..],
        &pocket3().types()[pocket3().ligand_range()]
    );
    let w3 = hydration::w3_solvated(&benzene);
    let w3_decoupling = w3.decoupling();
    let w3_template = w3.dynamics(&benzene);
    let free = water_surroundings(&benzene, &w3);
    let free_labels = water_labels(&free);
    let qw = am1bcc(&freesolv::component_graph(&benzene.component)).charges;
    let free_vdw = VdwSwap::new(&free, &free_labels);
    let free_frames = FramesDelta::new(&free, &free_labels, &qw);
    let melted: std::cell::OnceCell<Vec<[f64; 3]>> = std::cell::OnceCell::new();
    let water_start = || {
        melted
            .get_or_init(|| {
                let mut at = w3.at.clone();
                let tm = std::time::Instant::now();
                hydration::w3_melt(&w3_decoupling, &w3_template, KELVIN, &mut at);
                log(
                    &staged,
                    tag,
                    &format!(
                        "W3's box, {} waters, melted in {:.0} s",
                        w3.waters,
                        tm.elapsed().as_secs_f64()
                    ),
                );
                at
            })
            .clone()
    };
    let water = StagedPhase {
        name: "water",
        base: &w3_decoupling,
        restraint: None,
        charges: &qw,
        vdw: &free_vdw,
        frames: &free_frames,
        template: &w3_template,
        lambda: Lambda::COUPLED,
        identity: format!(
            "built {:016x} melted w3_melt {}",
            positions_digest(&w3.at),
            exact(KELVIN)
        ),
        start: &water_start,
        log: tag,
    };
    let records_w = staged_states(&staged, tag, &water, p);

    let s = whole_complex();
    let leg = production_leg(&s, saved_start(&s, &dir));
    log(
        &staged,
        tag,
        &format!(
            "the leg's box and start read back in {:.0} s",
            t.elapsed().as_secs_f64()
        ),
    );
    let g = grid();
    assert_eq!(g[0], Lambda::new(0.0, 1.0, 1.0));
    let bound = complex_surroundings(&s);
    let bound_labels = complex_labels(&s, the_system(), &bound);
    let qb = am1bcc(&ligand_graph(s.binding())).charges;
    let bound_vdw = VdwSwap::new(&bound, &bound_labels);
    let bound_frames = FramesDelta::new(&bound, &bound_labels, &qb);
    let complex_start = || leg.start.clone();
    let complex = StagedPhase {
        name: "complex",
        base: &leg.decoupling,
        restraint: Some(leg.restraint),
        charges: &qb,
        vdw: &bound_vdw,
        frames: &bound_frames,
        template: &leg.template,
        lambda: g[0],
        identity: format!(
            "start {:016x} restraint {:016x} flexibility {:?}",
            positions_digest(&leg.start),
            restraint_digest(&leg.restraint),
            leg.flexibility
        ),
        start: &complex_start,
        log: tag,
    };
    let records_b = staged_states(&staged, tag, &complex, p);

    println!(
        "\n**Variant 4 staged**, H(s) = H_UFF + s ΔU₄ at s = {STAGES:?}, protocol {p:?}: {:.0} s \
         in all",
        t.elapsed().as_secs_f64()
    );
    let b = report_staged("bound", &records_b, kt);
    let w = report_staged("in water", &records_w, kt);
    let ddg = b.bar.0 - w.bar.0;
    let ddg_error = (b.bar.1 * b.bar.1 + w.bar.1 * w.bar.1).sqrt();
    let (hyd, hyd_error) = w3_hydration();
    let hydration = hyd + w.bar.0;
    let bind = W4_BINDING.0 + ddg;
    println!(
        "\n| | staged (BAR) | its halves | TI | Zwanzig from s = 0, from s = 1 | the one step on the \
         saved frames | experiment |\n| --- | --- | --- | --- | --- | --- | --- |"
    );
    println!(
        "| ΔΔG_bind | {ddg:+.3} ± {ddg_error:.3} | {:+.3}, {:+.3} | {:+.3} | {:+.3}, {:+.3} | {:+.3} \
         ± {:.3} | |",
        b.halves[0] - w.halves[0],
        b.halves[1] - w.halves[1],
        b.ti.0 - w.ti.0,
        b.one_step[0] - w.one_step[0],
        b.one_step[1] - w.one_step[1],
        VARIANT_4_ONE_STEP.0,
        VARIANT_4_ONE_STEP.1
    );
    println!(
        "| ΔG°_bind = W4's {:+.3} + ΔΔG_bind | {bind:+.3} ± {:.3} | | {:+.3} | | {:+.2} | {:+.2} ± \
         {:.2} |",
        W4_BINDING.0,
        (W4_BINDING.1 * W4_BINDING.1 + ddg_error * ddg_error).sqrt(),
        W4_BINDING.0 + b.ti.0 - w.ti.0,
        VARIANT_4_IMPLIED[0],
        EXPERIMENT.0,
        EXPERIMENT.1
    );
    println!(
        "| ΔG_hyd = W3's {hyd:+.3} + ΔG_swap(water) | {hydration:+.3} ± {:.3} | | {:+.3} | | {:+.2} \
         | {:+.2} ± {:.2} |",
        (hyd_error * hyd_error + w.bar.1 * w.bar.1).sqrt(),
        hyd + w.ti.0,
        VARIANT_4_IMPLIED[1],
        freesolv::EXPERIMENT.0,
        freesolv::EXPERIMENT.1
    );
    println!(
        "\nstaged − one step: ΔΔG_bind {:+.3} ({:.1} of their errors combined); the lowest overlap \
         {:.3} bound, {:.3} in water",
        ddg - VARIANT_4_ONE_STEP.0,
        (ddg - VARIANT_4_ONE_STEP.0).abs()
            / (ddg_error * ddg_error + VARIANT_4_ONE_STEP.1 * VARIANT_4_ONE_STEP.1).sqrt(),
        b.overlap,
        w.overlap
    );

    // s = 0 against what it should reproduce, if nothing it runs on has changed since.
    let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
    if p == PROTOCOL {
        if let Ok(text) = std::fs::read_to_string(dir.join("window_00.txt")) {
            let w0 = parse(&text);
            let theirs: Vec<f64> = w0.energies.iter().map(|e| e[0]).collect();
            let (same, of) = to_the_bit(&bits(&records_b[0].coupling), &bits(&theirs));
            println!(
                "\ns = 0, bound, against the leg's window 0: the coupling at its state the same to \
                 the bit at {same} of {of} samples"
            );
        }
    }
    if p == FRAMES {
        for (name, records) in [("complex", &records_b), ("water", &records_w)] {
            let path = dir.join("frames").join(format!("{name}.frames"));
            if let Some(Ok(f)) = read_frames(&path) {
                let theirs: Vec<u64> = f.at.iter().map(|at| positions_digest(at)).collect();
                let (same, of) = to_the_bit(&records[0].digests, &theirs);
                println!(
                    "\ns = 0, {name}, against the saved frames: the positions the same to the bit \
                     at {same} of {of} samples; the steps {}",
                    if records[0].steps == f.steps {
                        "the same"
                    } else {
                        "**not the same**"
                    }
                );
            }
        }
    }
    log(
        &staged,
        tag,
        &format!("done in {:.0} s", t.elapsed().as_secs_f64()),
    );
    for r in records_b.iter().chain(&records_w) {
        assert_eq!(r.du.len(), p.samples);
        assert!(r.held <= 1.0, "state {}: {}", r.state, r.held);
    }
}

/// **Variant 4's one step checked by staged sampling, at W3's and the leg's protocol** (release,
/// one core; not run). `H(s) = H_UFF + s ΔU₄` ([`Staged`]) sampled at `s` = 0, 0.25, 0.5, 0.75
/// and 1 in W3's water and in the complex, each state 20 ps discarded and 200 ps sampled every
/// 100 fs ([`PROTOCOL`]), and BAR along each chain ([`staged_measured`]). Written to
/// `staged/water_long_K.txt` and `staged/complex_long_K.txt` under `PANTOMETRY_W4_DIR`, each state
/// read back if a run completed it; progress in `staged/long.log`. `PANTOMETRY_W4_DIR=… cargo test
/// --release -p pantometry-forcefield --test benzene_bound_in_tip3p
/// variant_4_staged_measured -- --ignored --nocapture`.
#[test]
#[ignore = "five states of the complex and of W3's box, 110 000 steps each: about two days with --release, not run"]
fn variant_4_staged_measured() {
    staged_measured("long", PROTOCOL);
}

/// **The same at the frames' protocol**, 4 ps discarded and 16 ps sampled every 50 fs per state
/// ([`FRAMES`]), its `s` = 0 the saved frames' own run: a fifth of a day rather than two
/// (`staged/*_short_K.txt`, `staged/short.log`). `PANTOMETRY_W4_DIR=… cargo test --release -p
/// pantometry-forcefield --test benzene_bound_in_tip3p variant_4_staged_short_measured --
/// --ignored --nocapture`.
///
/// **Run once**, at `a464806`, exit 0 on 2026-10-20 after 17 565 s. ΔΔG_bind is +4.370 ± 0.063 by
/// BAR along the chain, against the one step's +3.953. ΔG°_bind is −4.61 ± 0.29 against −5.19,
/// and ΔG_hyd −0.82 ± 0.15 against −0.90 (FreeSolv's GAFF/AM1-BCC calculation: −0.81). Every
/// interval overlapped at 0.485 or more, and s = 0 was the saved frames to the bit. Nothing here
/// asserts it.
#[test]
#[ignore = "five states of the complex and of W3's box, 10 000 steps each: about five hours with --release, run once (2026-10-20)"]
fn variant_4_staged_short_measured() {
    staged_measured("short", FRAMES);
}

/// Milliseconds a step of `potential` from `start` under `template` at 2 fs in a 1 ps⁻¹ bath, over
/// `steps` steps after ten that are not timed.
fn ms_per_step(
    potential: &dyn pantometry_forcefield::Potential,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    steps: u64,
) -> f64 {
    let seed = 0x5_7A6E;
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 1e12,
            seed,
        })
        .thermalised(start, KELVIN, seed);
    let mut at = start.to_vec();
    md.prepare(potential, &at);
    for _ in 0..10 {
        md.step(potential, &mut at, 2.0 * FS);
    }
    let t = std::time::Instant::now();
    for _ in 0..steps {
        md.step(potential, &mut at, 2.0 * FS);
    }
    t.elapsed().as_secs_f64() * 1e3 / steps as f64
}

/// One phase's costs, milliseconds: a step of the field, a step of `H(0.5)`, a sample's records
/// and a check's.
struct StagedCost {
    field: f64,
    staged: f64,
    sample: f64,
    check: f64,
}

/// [`StagedCost`] for a phase from `start`: `steps` of the field and of `H(0.5)`, then five
/// samples' records and two checks, timed.
fn staged_cost(phase: &StagedPhase, start: &[[f64; 3]], steps: u64) -> StagedCost {
    let field = ms_per_step(
        &AtLambda {
            hamiltonian: phase.base,
            lambda: phase.lambda,
        },
        start,
        phase.template,
        steps,
    );
    let half = stage_potential(phase, 2);
    assert_eq!(half.s, 0.5);
    let staged = ms_per_step(&half, start, phase.template, steps);
    let t = std::time::Instant::now();
    for _ in 0..5 {
        let change = phase.frames.at(start);
        let coupling = phase.base.couplings(start, &[phase.lambda])[0].energy;
        assert!(change.total().is_finite() && coupling.is_finite());
        let _ = positions_digest(start);
    }
    let sample = t.elapsed().as_secs_f64() * 1e3 / 5.0;
    let change = phase.frames.at(start);
    let t = std::time::Instant::now();
    for _ in 0..2 {
        let held = staged_held(&half, start, &change);
        assert!(held <= 1.0, "{}: {held}", phase.name);
    }
    StagedCost {
        field,
        staged,
        sample,
        check: t.elapsed().as_secs_f64() * 1e3 / 2.0,
    }
}

/// Hours for one phase's five states at `p` at `cost`.
fn staged_hours(cost: &StagedCost, p: Protocol) -> f64 {
    let checks = (p.samples - 1) / STAGE_CHECK_EVERY + 1;
    let ms = p.steps_per_window() as f64 * cost.staged
        + p.samples as f64 * cost.sample
        + checks as f64 * cost.check;
    STAGES.len() as f64 * ms / 3.6e6
}

/// **What the staged check costs** (release, one core, about fifteen minutes; writes nothing):
/// the leg's box built and its start read back, W3's box built and melted, then on each a few
/// hundred steps of the field and of `H(0.5)` ([`Staged`]), a sample's records and a check's, and
/// from them the whole run at [`PROTOCOL`] and at [`FRAMES`]. `PANTOMETRY_W4_DIR=… cargo test
/// --release -p pantometry-forcefield --test benzene_bound_in_tip3p the_staged_step_costed --
/// --ignored --nocapture`.
#[test]
#[ignore = "the whole complex built and W3's box melted, a few hundred steps of each: about fifteen minutes with --release"]
fn the_staged_step_costed() {
    let dir = results_dir();
    let t = std::time::Instant::now();
    let s = whole_complex();
    let leg = production_leg(&s, saved_start(&s, &dir));
    let built = t.elapsed().as_secs_f64();
    let bound = complex_surroundings(&s);
    let bound_labels = complex_labels(&s, the_system(), &bound);
    let qb = am1bcc(&ligand_graph(s.binding())).charges;
    let bound_vdw = VdwSwap::new(&bound, &bound_labels);
    let bound_frames = FramesDelta::new(&bound, &bound_labels, &qb);
    let no_start = || -> Vec<[f64; 3]> { panic!("costing runs no state") };
    let complex = StagedPhase {
        name: "complex",
        base: &leg.decoupling,
        restraint: Some(leg.restraint),
        charges: &qb,
        vdw: &bound_vdw,
        frames: &bound_frames,
        template: &leg.template,
        lambda: grid()[0],
        identity: String::new(),
        start: &no_start,
        log: "cost",
    };
    let set_up = t.elapsed().as_secs_f64();
    let c = staged_cost(&complex, &leg.start, 200);
    println!(
        "complex: {} atoms, built and read back in {built:.0} s, taken apart in {:.0} s; {} \
         partners of benzene's van der Waals",
        leg.start.len(),
        set_up - built,
        bound_vdw.partners.len()
    );

    let benzene = hydration::Benzene::new(Some(pocket3().ligand_force_field().charges().to_vec()));
    let w3 = hydration::w3_solvated(&benzene);
    let w3_decoupling = w3.decoupling();
    let w3_template = w3.dynamics(&benzene);
    let free = water_surroundings(&benzene, &w3);
    let free_labels = water_labels(&free);
    let qw = am1bcc(&freesolv::component_graph(&benzene.component)).charges;
    let free_vdw = VdwSwap::new(&free, &free_labels);
    let free_frames = FramesDelta::new(&free, &free_labels, &qw);
    let water = StagedPhase {
        name: "water",
        base: &w3_decoupling,
        restraint: None,
        charges: &qw,
        vdw: &free_vdw,
        frames: &free_frames,
        template: &w3_template,
        lambda: Lambda::COUPLED,
        identity: String::new(),
        start: &no_start,
        log: "cost",
    };
    let tm = std::time::Instant::now();
    let mut at = w3.at.clone();
    hydration::w3_melt(&w3_decoupling, &w3_template, KELVIN, &mut at);
    let melt = tm.elapsed().as_secs_f64();
    let w = staged_cost(&water, &at, 500);
    println!(
        "water: {} atoms, melted in {melt:.0} s\n\n| | the field (ms a step) | H(0.5) (ms a step) | \
         over the field | a sample's records (ms) | a check (ms) |\n| --- | --- | --- | --- | --- | \
         --- |",
        at.len()
    );
    for (name, x) in [("complex", &c), ("water", &w)] {
        println!(
            "| {name} | {:.1} | {:.1} | {:+.1}% | {:.1} | {:.1} |",
            x.field,
            x.staged,
            100.0 * (x.staged / x.field - 1.0),
            x.sample,
            x.check
        );
    }
    println!(
        "\n| protocol | steps a state | complex (h) | water (h) | both, with the builds and the \
         melt (h) |\n| --- | --- | --- | --- | --- |"
    );
    for (name, p) in [
        ("PROTOCOL, 20 + 200 ps", PROTOCOL),
        ("FRAMES, 4 + 16 ps", FRAMES),
    ] {
        let (hc, hw) = (staged_hours(&c, p), staged_hours(&w, p));
        println!(
            "| {name} | {} | {hc:.1} | {hw:.1} | {:.1} |",
            p.steps_per_window(),
            hc + hw + (set_up + melt) / 3600.0
        );
    }
}

/// A small box for the staged checks: the box taken apart, its labels, the decoupling `H(s)` is a
/// change of and its restraint, FreeSolv's charges in benzene's order, the dynamics, the state the
/// production samples it at, and the positions as built.
struct SmallStaged {
    name: &'static str,
    sur: Surroundings,
    labels: Labels,
    base: PeriodicDecoupling,
    restraint: Option<Boresch>,
    charges: Vec<f64>,
    template: MolecularDynamics,
    lambda: Lambda,
    at: Vec<[f64; 3]>,
}

impl SmallStaged {
    /// The box as a [`StagedPhase`], started from `start`.
    fn phase<'a>(
        &'a self,
        vdw: &'a VdwSwap<'a>,
        frames: &'a FramesDelta<'a>,
        start: &'a dyn Fn() -> Vec<[f64; 3]>,
    ) -> StagedPhase<'a> {
        StagedPhase {
            name: "small",
            base: &self.base,
            restraint: self.restraint,
            charges: &self.charges,
            vdw,
            frames,
            template: &self.template,
            lambda: self.lambda,
            identity: format!("built {:016x}", positions_digest(&self.at)),
            start,
            log: "staged",
        }
    }
}

/// The 3 Å pocket in water through the mesh, the leg's flexibility, the restraint on and the
/// leg's state, `λ_r` = 0; and W3's small box through the classical sum, coupled.
fn small_staged() -> Vec<SmallStaged> {
    let s = small3(0x5EC7);
    let sur = complex_surroundings(&s);
    let labels = complex_labels(&s, the_system(), &sur);
    let restraint = anchors(s.binding(), &s.positions()[..s.complex_atoms()]);
    let pocket = SmallStaged {
        name: "3 Å pocket, mesh",
        base: s.decoupling().with_restraint(restraint),
        restraint: Some(restraint),
        charges: am1bcc(&ligand_graph(s.binding())).charges,
        template: s.dynamics(LEG_FLEXIBILITY),
        lambda: grid()[0],
        at: s.positions().to_vec(),
        sur,
        labels,
    };
    assert!(pocket.sur.whole.field().ewald().mesh().is_some());
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC8, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    let free = water_surroundings(&b, &w);
    let free_labels = water_labels(&free);
    let water = SmallStaged {
        name: "W3's small box, Ewald",
        base: w.decoupling(),
        restraint: None,
        charges: am1bcc(&freesolv::component_graph(&b.component)).charges,
        template: w.dynamics(&b),
        lambda: Lambda::COUPLED,
        at: w.at.clone(),
        sur: free,
        labels: free_labels,
    };
    assert!(water.sur.whole.field().ewald().mesh().is_none());
    vec![pocket, water]
}

/// The positions' bits, atom by atom.
fn position_bits(at: &[[f64; 3]]) -> Vec<[u64; 3]> {
    at.iter().map(|p| p.map(f64::to_bits)).collect()
}

/// **The forces of `ΔU₄` are its gradient**: on the 3 Å pocket through the mesh and W3's small box
/// through the classical sum, `F_H(1) − F_H(0)` from the dynamics' own [`Staged`] against central
/// differences of `ΔU₄` **as the frames compute it** ([`FramesDelta`], the couplings' pass and the
/// van der Waals sums: another path, with no force in it), within [`finite_difference`]'s
/// tolerance at h = 10⁻⁵ Å plus `64 ε` of the larger forces for the two evaluations' rounding.
/// Every axis of every benzene atom — its own Coulomb pairs taken back, the van der Waals swapped,
/// the cross Coulomb — and of the protein or water atom nearest it; and on the pocket **an atom
/// beyond the cutoff of every benzene atom**, whose force from `ΔU₄` is the reciprocal sum's alone
/// and must be ten times the tolerance or more, so that a run without it would fail. No pair a moved atom is in lies within 8 h of the cutoff.
/// Both boxes are printed before either is asserted.
#[test]
fn the_staged_forces_are_the_gradient_of_its_change() {
    let h = 1e-5 * ANGSTROM;
    let mut results = Vec::new();
    for b in small_staged() {
        let vdw = VdwSwap::new(&b.sur, &b.labels);
        let frames = FramesDelta::new(&b.sur, &b.labels, &b.charges);
        let one = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, 1.0);
        let zero = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, 0.0);
        let n = b.at.len();
        let (mut f1, mut f0) = (vec![[0.0; 3]; n], vec![[0.0; 3]; n]);
        one.evaluate(&b.at, &mut f1);
        zero.evaluate(&b.at, &mut f0);
        let rounding = 64.0 * EPS * (largest_force(&f1) + largest_force(&f0));
        let change = frames.at(&b.at);
        let scale = b.sur.scale(&b.at) + change.magnitude;
        let w: Vec<[f64; 3]> = b.at.iter().map(|p| b.sur.cell.wrap(*p)).collect();
        let ligand = b.sur.ligand.clone();
        let distance = |i: usize, j: usize| b.sur.separation(&w, i.min(j), i.max(j)).1.sqrt();
        let to_benzene = |j: usize| {
            ligand
                .clone()
                .map(|g| distance(g, j))
                .fold(f64::INFINITY, f64::min)
        };
        let rest: Vec<usize> = (0..n).filter(|&j| !ligand.contains(&j)).collect();
        let near = *rest
            .iter()
            .filter(|&&j| b.sur.vdw[j][1] != 0.0)
            .min_by(|&&x, &&y| to_benzene(x).total_cmp(&to_benzene(y)))
            .expect("an atom near benzene");
        let far = rest.iter().copied().find(|&j| {
            b.sur.whole.field().charges()[j] != 0.0 && to_benzene(j) > b.sur.cutoff + ANGSTROM
        });
        let mut atoms: Vec<usize> = ligand.clone().collect();
        atoms.push(near);
        atoms.extend(far);
        for &i in &atoms {
            let partners: Vec<usize> = if ligand.contains(&i) {
                rest.clone()
            } else {
                ligand.clone().collect()
            };
            for j in partners {
                let r = distance(i, j);
                assert!(
                    (r - b.sur.cutoff).abs() > 8.0 * h,
                    "{}: atoms {i} and {j} at the cutoff",
                    b.name
                );
            }
        }
        let (mut worst, mut far_force, mut largest) = (0.0f64, 0.0f64, 0.0f64);
        for &i in &atoms {
            for x in 0..3 {
                let energy = |dx: f64| {
                    let mut p = b.at.clone();
                    p[i][x] += dx;
                    (frames.at(&p).total(), scale)
                };
                let (d, tol) = finite_difference(energy, h);
                let f = f1[i][x] - f0[i][x];
                let ratio = (f + d).abs() / (tol + rounding);
                if ratio.is_nan() || ratio > worst {
                    worst = ratio;
                }
                largest = largest.max(f.abs() / (tol + rounding));
                if Some(i) == far {
                    far_force = far_force.max(f.abs() / (tol + rounding));
                }
            }
        }
        println!(
            "{}: {} atoms by three axes, the worst {worst:.2e} of the tolerance, the largest force \
             {largest:.1e} times it; atom {near} nearest benzene at {:.2} Å; the far atom {far:?} \
             at {:.2} Å, its force {far_force:.1e} times the tolerance",
            b.name,
            atoms.len(),
            to_benzene(near) / ANGSTROM,
            far.map_or(f64::NAN, |j| to_benzene(j) / ANGSTROM)
        );
        results.push((b.name, worst, far, far_force));
    }
    for (name, worst, _, _) in &results {
        assert!(*worst <= 1.0, "{name}: {worst}");
    }
    let (name, _, far, far_force) = &results[0];
    assert!(
        far.is_some(),
        "{name}: an atom beyond the cutoff of every benzene atom"
    );
    assert!(
        *far_force > 10.0,
        "{name}: the far atom's force is the reciprocal sum's, far above the tolerance: \
         {far_force}"
    );
}

/// **At `s` = 0, `H(s)` is the field, to the bit**: on both small boxes, [`Staged`] at `s` = 0
/// holds a decoupling equal to the field's — charges, restraint, soft core, group, pairs — and
/// gives its energy and every force with the same bits, at the production's state and at
/// `λ = (1, 0.5, 0.7)`, where the pocket's restraint, the charges and the soft core all count; and
/// three samples of the dynamics under it, through [`sample_on`], are the field's own
/// ([`sample_frames`]) to the bit, while under `H(0.5)` they are not.
#[test]
fn the_staged_hamiltonian_at_zero_is_the_field() {
    for b in small_staged() {
        let vdw = VdwSwap::new(&b.sur, &b.labels);
        let n = b.at.len();
        for lambda in [b.lambda, Lambda::new(1.0, 0.5, 0.7)] {
            let zero = Staged::new(&b.base, b.restraint, &b.charges, &vdw, lambda, 0.0);
            assert!(
                zero.charged == b.base,
                "{}: the decoupling at s = 0 is the field's",
                b.name
            );
            let (mut f, mut g) = (vec![[0.0; 3]; n], vec![[0.0; 3]; n]);
            let e = zero.evaluate(&b.at, &mut f);
            let e0 = Alchemical::energy_and_forces(&b.base, &b.at, lambda, &mut g);
            assert_eq!(e.to_bits(), e0.to_bits(), "{}: {lambda:?}", b.name);
            assert!(
                position_bits(&f) == position_bits(&g),
                "{}: every force, {lambda:?}",
                b.name
            );
        }
        let p = Protocol {
            time_step: 0.5e-15,
            temperature: KELVIN,
            friction: 50e12,
            equilibration: 2,
            stride: 2,
            samples: 3,
            seed: 0x5A60,
        };
        let mut field = Vec::new();
        sample_frames(&b.base, &b.at, &b.template, p, 0, b.lambda, |s, at| {
            field.push((s, position_bits(at)))
        });
        let run = |s: f64| {
            let st = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, s);
            let mut v = Vec::new();
            sample_on(&st, &b.at, &b.template, p, 0, |s, at| {
                v.push((s, position_bits(at)))
            });
            v
        };
        assert_eq!(field.len(), 3);
        assert!(run(0.0) == field, "{}: the dynamics at s = 0", b.name);
        assert!(run(0.5) != field, "{}: and at s = 0.5 not", b.name);
        println!(
            "{}: H(0) is the field to the bit, energy, forces and three samples",
            b.name
        );
    }
}

/// **`ΔU₄` in the dynamics is the frames' `ΔU₄`**: on both small boxes, as built and after four
/// steps at `s` = 1, `H(s) − H(0)` from [`Staged`] is `s ΔU₄` from [`FramesDelta`] at `s` = 0.25,
/// 0.5, 0.75 and 1, within [`staged_held`]'s bound — so `H` is linear in `s` and its states are
/// the ones BAR is told they are — and the van der Waals half of [`Staged`]'s own sum
/// ([`VdwSwap`]) is the frames' within the summation bound. Neither half of `ΔU₄` is zero.
#[test]
fn the_staged_change_is_the_frames_change() {
    let mut results = Vec::new();
    for b in small_staged() {
        let vdw = VdwSwap::new(&b.sur, &b.labels);
        let frames = FramesDelta::new(&b.sur, &b.labels, &b.charges);
        let one = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, 1.0);
        let p = Protocol {
            time_step: 0.5e-15,
            temperature: KELVIN,
            friction: 50e12,
            equilibration: 0,
            stride: 4,
            samples: 1,
            seed: 0x5A61,
        };
        let mut moved = Vec::new();
        sample_on(&one, &b.at, &b.template, p, 0, |_, at| moved = at.to_vec());
        for (where_, at) in [("as built", &b.at), ("after four steps", &moved)] {
            let change = frames.at(at);
            let v = vdw.evaluate(at, 0.0, None);
            let vdw_ratio = (v.total() - change.vdw).abs()
                / summation_bound(v.terms + change.terms, v.magnitude + change.magnitude);
            let mut worst = 0.0f64;
            for &s in &STAGES[1..] {
                let st = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, s);
                let held = staged_held(&st, at, &change);
                if held.is_nan() || held > worst {
                    worst = held;
                }
            }
            println!(
                "{} {where_}: ΔU₄ {:+.5} kcal/mol (van der Waals {:+.5}, electrostatics {:+.5}); \
                 H(s) − H(0) − s ΔU₄ at worst {worst:.1e} of its bound, the van der Waals sums \
                 {vdw_ratio:.1e}",
                b.name,
                kcal(change.total()),
                kcal(change.vdw),
                kcal(change.elec)
            );
            results.push((b.name, where_, change, worst, vdw_ratio));
        }
    }
    for (name, where_, change, worst, vdw_ratio) in results {
        assert!(worst <= 1.0, "{name} {where_}: {worst}");
        assert!(vdw_ratio <= 1.0, "{name} {where_}: {vdw_ratio}");
        assert!(change.vdw != 0.0 && change.elec != 0.0, "{name} {where_}");
    }
}

/// **The `s`-states are what their headers say**: [`STAGES`] runs from 0 to 1 upwards and
/// [`STAGE_INDEX`] gives `s` = 0 window 0's seed and every other state one past [`grid`]'s
/// candidates, all distinct; on the 3 Å pocket each state's header ([`staged_header`]) reads back
/// as its state, its `s` to the bit and its seed; the potential the run builds for it
/// ([`stage_potential`]) is `H` at that `s`, and at that `s` it is `H(0) + s ΔU₄` within
/// [`staged_held`]'s bound; and each header is state 0's but for its state, `s`, index and seed.
#[test]
fn the_staged_states_are_what_their_headers_say() {
    assert_eq!(STAGES[0], 0.0);
    assert_eq!(STAGES[STAGES.len() - 1], 1.0);
    assert!(STAGES.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(STAGE_INDEX[0], 0);
    assert!(STAGE_INDEX[1..].iter().all(|&i| i >= grid().len()));
    let distinct: std::collections::BTreeSet<usize> = STAGE_INDEX.iter().copied().collect();
    assert_eq!(distinct.len(), STAGE_INDEX.len());
    let b = small_staged().swap_remove(0);
    let vdw = VdwSwap::new(&b.sur, &b.labels);
    let frames = FramesDelta::new(&b.sur, &b.labels, &b.charges);
    let start = || b.at.clone();
    let phase = b.phase(&vdw, &frames, &start);
    let change = frames.at(&b.at);
    let headers: Vec<String> = (0..STAGES.len())
        .map(|k| staged_header(&phase, PROTOCOL, k))
        .collect();
    for (k, head) in headers.iter().enumerate() {
        let (state, s, index, seed) = staged_state_of(head).expect("a state");
        assert_eq!(state, k, "{head}");
        assert_eq!(s.to_bits(), STAGES[k].to_bits(), "{head}");
        assert_eq!(index, STAGE_INDEX[k], "{head}");
        assert_eq!(seed, window_seed(PROTOCOL.seed, index), "{head}");
        let st = stage_potential(&phase, k);
        assert_eq!(st.s.to_bits(), s.to_bits(), "state {k}");
        let held = staged_held(&st, &b.at, &change);
        assert!(held <= 1.0, "state {k}: {held}");
        let own = |j: usize| {
            format!(
                "state {j} s {} of {} index {} seed {:016x} ",
                exact(STAGES[j]),
                STAGES.map(exact).join(","),
                STAGE_INDEX[j],
                window_seed(PROTOCOL.seed, STAGE_INDEX[j])
            )
        };
        assert!(head.contains(&own(k)), "{head}");
        assert_eq!(head.replace(&own(k), &own(0)), headers[0], "state {k}");
    }
    println!("{}", headers[2]);
}

/// **A killed staged run resumes to the same bits**, through the production's own
/// [`staged_state`] on the 3 Å pocket at a protocol of eight 0.5 fs steps and three samples a
/// state, in temporary directories: every state run uninterrupted; and states 0 and 1 run, then
/// what a run killed during state 2 leaves — its `.partial`, a header and one sample — then the
/// whole phase again. States 0 and 1 are read, not run (their files' bytes unchanged), the start is
/// asked for only by the states that run, the partial is begun afresh and renamed, and **every
/// number of every state is the uninterrupted run's to the bit**; and state 2's positions are
/// [`sample_on`] of `H` at the `s` and seed index its header says, so a state runs what its header
/// says it ran. And each guard refuses by name:
/// a state asked for at another seed, and a state file cut short at every line.
#[cfg(not(target_family = "wasm"))]
#[test]
fn a_staged_run_resumes_to_the_same_bits() {
    let b = small_staged().swap_remove(0);
    let vdw = VdwSwap::new(&b.sur, &b.labels);
    let frames = FramesDelta::new(&b.sur, &b.labels, &b.charges);
    let started = std::cell::Cell::new(0usize);
    let start = || {
        started.set(started.get() + 1);
        b.at.clone()
    };
    let phase = b.phase(&vdw, &frames, &start);
    let p = Protocol {
        time_step: 0.5e-15,
        temperature: KELVIN,
        friction: 50e12,
        equilibration: 2,
        stride: 2,
        samples: 3,
        seed: 0x5A62,
    };
    let whole = fresh_dir("staged-whole");
    let full = staged_states(&whole, "t", &phase, p);
    assert_eq!(started.get(), STAGES.len());
    let killed = fresh_dir("staged-killed");
    for k in 0..2 {
        staged_state(&killed, "t", &phase, p, k);
    }
    let bytes = |d: &std::path::Path, k: usize| std::fs::read(d.join(format!("small_t_{k}.txt")));
    let before = [bytes(&killed, 0).unwrap(), bytes(&killed, 1).unwrap()];
    let text = std::fs::read_to_string(whole.join("small_t_2.txt")).unwrap();
    let partial: String = text.lines().take(2).map(|l| format!("{l}\n")).collect();
    std::fs::write(killed.join("small_t_2.partial"), partial).unwrap();
    started.set(0);
    let resumed = staged_states(&killed, "t", &phase, p);
    assert_eq!(started.get(), STAGES.len() - 2, "states 0 and 1 read back");
    assert!(bytes(&killed, 0).unwrap() == before[0] && bytes(&killed, 1).unwrap() == before[1]);
    assert!(!killed.join("small_t_2.partial").exists());
    let all = |r: &StagedRecord| -> Vec<u64> {
        let mut v: Vec<u64> = r.numbers().iter().map(|x| x.to_bits()).collect();
        v.pop(); // the seconds
        v
    };
    for (a, c) in full.iter().zip(&resumed) {
        assert_eq!((a.state, a.s.to_bits()), (c.state, c.s.to_bits()));
        assert_eq!(a.du.len(), p.samples);
        assert!(a.held <= 1.0, "state {}: {}", a.state, a.held);
        assert!(all(a) == all(c), "state {}", a.state);
    }
    // Each state's samples are its own: the states differ. And a state ran at its header's `s`
    // and seed: state 2's positions are `sample_on` of `H` at the `s` and index its header says.
    assert!(full[0].digests != full[4].digests);
    let (_, s, index, _) = staged_state_of(&staged_header(&phase, p, 2)).expect("a state");
    let st = Staged::new(&b.base, b.restraint, &b.charges, &vdw, b.lambda, s);
    let mut digests = Vec::new();
    sample_on(&st, &b.at, &b.template, p, index, |_, at| {
        digests.push(positions_digest(at))
    });
    assert!(
        digests == full[2].digests,
        "state 2 ran at its header's s and seed"
    );
    println!(
        "resumed: {} states, {} numbers each the uninterrupted run's",
        resumed.len(),
        all(&resumed[0]).len()
    );

    let other = Protocol {
        seed: p.seed + 1,
        ..p
    };
    let why = refusal(|| drop(staged_state(&killed, "t", &phase, other, 0)));
    assert!(why.contains("was run differently"), "{why}");
    let file = std::fs::read_to_string(whole.join("small_t_3.txt")).unwrap();
    let cut = fresh_dir("staged-cut");
    for k in 1..file.len() - 1 {
        if !file.is_char_boundary(k) {
            continue;
        }
        assert!(
            parse_staged(&file[..k]).is_err(),
            "cut at byte {k} of {} read as whole",
            file.len()
        );
    }
    let short = &file[..file.len() - 3];
    std::fs::write(
        cut.join("small_t_0.txt"),
        short.replacen(
            &staged_header(&phase, p, 3),
            &staged_header(&phase, p, 0),
            1,
        ),
    )
    .unwrap();
    let why = refusal(|| drop(staged_state(&cut, "t", &phase, p, 0)));
    assert!(why.contains("not whole"), "{why}");
    for d in [whole, killed, cut] {
        let _ = std::fs::remove_dir_all(d);
    }
}

// ---------------------------------------------------------------------------------------------
// Step 4c-2's direct water route: benzene's model the staged run's

/// **The direct route's benzene is variant 4's, to the bit**: `benzene_partitioned_into_octanol.rs`
/// decouples GAFF/AM1-BCC benzene from water and from octanol directly, with its model built by
/// [`gaff::benzene_in_water`] from `gaff.dat`'s `ca` and `ha` read again and FreeSolv's types. On
/// W3's small box, every atom's `[x, D]` is the one variant 4's staged run swapped in
/// ([`VdwSwap`]'s variant 3, [`gaff_benzene`] written from the same file separately), the waters
/// TIP3P's, the mixing Lorentz–Berthelot in both, and benzene's charges [`am1bcc`]'s, the waters'
/// TIP3P's, all to the bit. And **the direct field's cross van der Waals is the frames' variant 3**:
/// `U(coupled) − U(decoupled)`'s van der Waals is [`cross_sum`]'s pairs, by the frames' own loop and
/// mixing, within both sums' rounding. What is not the same, and is not meant to be: benzene's own
/// pairs, GAFF's, AM1-BCC's and 1-4 scaled there, UFF's and QEq's in the staged run.
#[test]
fn the_direct_routes_benzene_is_variant_4s() {
    let b = hydration::Benzene::new(None);
    let w = hydration::Solvated::new(&b, 3, 0x5EC6, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    let free = water_surroundings(&b, &w);
    let labels = water_labels(&free);
    let swap = VdwSwap::new(&free, &labels);
    let q4 = am1bcc(&freesolv::component_graph(&b.component)).charges;
    let direct = gaff::benzene_in_water(w.field.clone(), &b.component);
    let n = w.at.len();
    let bits2 =
        |v: &[[f64; 2]]| -> Vec<[u64; 2]> { v.iter().map(|x| x.map(f64::to_bits)).collect() };
    let bits = |v: &[f64]| -> Vec<u64> { v.iter().map(|x| x.to_bits()).collect() };
    assert_eq!(bits2(direct.van_der_waals()), bits2(&swap.to.vdw[..n]));
    assert_eq!(swap.to.mixing, Mixing::LorentzBerthelot);
    assert_eq!(
        direct.mixing(),
        pantometry_forcefield::Mixing::LorentzBerthelot
    );
    assert_eq!(bits(&direct.charges()[..12]), bits(&q4));
    assert_eq!(
        bits(&direct.charges()[12..]),
        bits(&w.field.charges()[12..])
    );
    // The cross van der Waals, by the field and by the frames.
    let d = PeriodicDecoupling::new(&direct, &w.group).expect("benzene decouples");
    let on = d.energy(&w.at, Lambda::COUPLED);
    let off = d.energy(&w.at, Lambda::new(0.0, 0.0, 0.0));
    let frames = cross_sum(&free, &w.at, &swap.to, CrossSum::default());
    let field = on.van_der_waals - off.van_der_waals;
    let bound = 64.0 * EPS * (magnitude(&on) + magnitude(&off))
        + summation_bound(frames.terms, frames.magnitude);
    println!(
        "the cross van der Waals: the field {:+.12e} kcal/mol, the frames {:+.12e}, {:.3} of the \
         bound; {} pairs",
        kcal(field),
        kcal(frames.pairs),
        (field - frames.pairs).abs() / bound,
        frames.terms
    );
    assert!((field - frames.pairs).abs() <= bound);
    assert!(
        frames.terms > 0 && frames.pairs.abs() > 1e3 * bound,
        "the sum is seen"
    );
}
