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
//!   melt and then the release, run as the stages say; and its schedule.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod protein;

use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::free_energy::{bennett_chain, refine_schedule, window_seed};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::water::{self, WaterBox};
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Binding, Boresch, Complex, Component, Coupling, Element, Estimate,
    EwaldParameters, Flexibility, Lambda, MolecularDynamics, PeriodicBox, PeriodicDecoupling,
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

/// Writes positions one atom a line, exactly (`{:e}` round-trips).
fn write_positions(path: &std::path::Path, at: &[[f64; 3]]) {
    let text: String = at
        .iter()
        .map(|p| format!("{:e} {:e} {:e}\n", p[0], p[1], p[2]))
        .collect();
    std::fs::write(path, text).expect("positions");
}

fn read_positions(path: &std::path::Path) -> Option<Vec<[f64; 3]>> {
    let text = std::fs::read_to_string(path).ok()?;
    Some(
        text.lines()
            .map(|l| {
                let v: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
                [v[0], v[1], v[2]]
            })
            .collect(),
    )
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

/// `stages` of [`SolvatedComplex::prepare`] on `potential` from the box as built, read from
/// `file` in `dir` if a run has written it, and each stage logged as it ends.
fn prepared(
    s: &SolvatedComplex,
    potential: &dyn pantometry_forcefield::Potential,
    stages: &[Stage],
    dir: &std::path::Path,
    file: &str,
    name: &str,
) -> Vec<[f64; 3]> {
    let path = dir.join(file);
    if let Some(at) = read_positions(&path) {
        assert_eq!(at.len(), s.positions().len(), "{file} is another system's");
        return at;
    }
    let t = std::time::Instant::now();
    let mut at = s.positions().to_vec();
    s.prepare(
        potential,
        stages,
        KELVIN,
        PREPARATION_SEED,
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
    write_positions(&path, &at);
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
/// and the couplings at `lambda` followed by every one of `states`. Returns where it ended. W3's.
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
    let seed = window_seed(p.seed, index);
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: p.temperature,
            friction: p.friction,
            seed,
        })
        .thermalised(start, p.temperature, seed);
    let potential = AtLambda {
        hamiltonian: h,
        lambda,
    };
    let mut at = start.to_vec();
    md.prepare(&potential, &at);
    let mut lambdas = vec![lambda];
    lambdas.extend_from_slice(states);
    let mut taken = 0;
    // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
    // builds on 1.78.
    #[allow(clippy::manual_is_multiple_of)]
    while taken < p.samples {
        md.step(&potential, &mut at, p.time_step);
        let s = md.steps();
        if s > p.equilibration && (s - p.equilibration) % p.stride == 0 {
            on(s, &h.couplings(&at, &lambdas));
            taken += 1;
        }
    }
    at
}

fn header(g: usize, l: Lambda, flexibility: Flexibility, atoms: usize) -> String {
    let p = PROTOCOL;
    format!(
        "# candidate {g} lambda {:e} {:e} {:e} protocol {:e} {:e} {:e} {} {} {} {} atoms {atoms} \
         flexibility {flexibility:?} states {}",
        l.restraint,
        l.electrostatics,
        l.van_der_waals,
        p.time_step,
        p.temperature,
        p.friction,
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        grid().len()
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

/// The leg's inputs.
struct Leg {
    decoupling: PeriodicDecoupling,
    template: MolecularDynamics,
    start: Vec<[f64; 3]>,
    flexibility: Flexibility,
}

/// Window `g`'s record: from its file if a run completed it with this protocol, or run now —
/// written sample by sample to a `.partial` file, renamed when complete. W3's.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn window(leg: &Leg, dir: &std::path::Path, g: usize) -> Recorded {
    use std::io::Write;
    let path = dir.join(format!("window_{g:02}.txt"));
    let grid = grid();
    let want = header(g, grid[g], leg.flexibility, leg.start.len());
    if let Ok(text) = std::fs::read_to_string(&path) {
        assert_eq!(
            text.lines().next().unwrap_or(""),
            want,
            "{} was run differently",
            path.display()
        );
        let r = parse(&text);
        assert_eq!(r.energies.len(), PROTOCOL.samples);
        return r;
    }
    let t = std::time::Instant::now();
    let partial = path.with_extension("partial");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
    writeln!(file, "{want}").unwrap();
    let total = PROTOCOL.steps_per_window();
    log(
        dir,
        "complex",
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
        PROTOCOL,
        g,
        grid[g],
        &grid,
        |s, c| {
            let mut line = format!(
                "{s} {:e} {:e} {:e}",
                c[0].gradient[0], c[0].gradient[1], c[0].gradient[2]
            );
            for x in &c[1..] {
                line.push_str(&format!(" {:e}", x.energy));
            }
            writeln!(file, "{line}").unwrap();
            taken += 1;
            if taken % 100 == 0 {
                file.flush().unwrap();
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    "complex",
                    &format!(
                        "candidate {g}: {taken} of {} samples, step {s} of {total}, {secs:.0} s, \
                         {:.1} ms per step",
                        PROTOCOL.samples,
                        secs * 1e3 / s as f64
                    ),
                );
            }
        },
    );
    let seconds = t.elapsed().as_secs_f64();
    writeln!(file, "# done seconds {seconds:e}").unwrap();
    drop(file);
    std::fs::rename(&partial, &path).expect("rename");
    log(
        dir,
        "complex",
        &format!("window at candidate {g} done in {seconds:.0} s"),
    );
    parse(&std::fs::read_to_string(&path).unwrap())
}

/// **Benzene decoupled from T4 lysozyme L99A in TIP3P, the complex leg, and the binding free
/// energy it closes with W3's solvent leg.** **Not run yet**: written, its pieces tested, its cost
/// measured (`the_cost_of_the_mobile_leg_measured`). The whole of 181L as [`whole_complex`]
/// builds it — the hydrogens relaxed in vacuum first, the neighbour list on — **the protein free**
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
#[ignore = "the complex leg: days with --release, and not to be run until production is started"]
fn the_complex_leg_measured() {
    let flexibility = LEG_FLEXIBILITY;
    let dir = results_dir();
    let _lock = Lock::take(&dir, "complex");
    let t = std::time::Instant::now();
    let s = whole_complex();
    let coupled = s.decoupling();
    let start = prepared(
        &s,
        &AtLambda {
            hamiltonian: &coupled,
            lambda: Lambda::COUPLED,
        },
        &SolvatedComplex::preparation(flexibility),
        &dir,
        "prepared.txt",
        "complex",
    );
    let restraint = anchors(s.binding(), &start[..s.complex_atoms()]);
    let kt = BOLTZMANN.to_si() * KELVIN;
    let release = kcal(restraint.release_free_energy(KELVIN));
    log(
        &dir,
        "complex",
        &format!(
            "{} atoms, flexibility {flexibility:?}; Boresch: receptor {:?}, ligand {:?}; r₀ {:.3} \
             Å, θ_A {:.1}°, θ_B {:.1}°; release to 1 M {release:+.3} kcal/mol (extended {:+.3})",
            start.len(),
            restraint.receptor,
            restraint.ligand,
            restraint.reference[0] / ANGSTROM,
            restraint.reference[1].to_degrees(),
            restraint.reference[2].to_degrees(),
            kcal(restraint.release_free_energy_extended(KELVIN))
        ),
    );
    let leg = Leg {
        decoupling: s.decoupling().with_restraint(restraint),
        template: s.dynamics(flexibility),
        start,
        flexibility,
    };
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
