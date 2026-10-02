//! Aspirin, every atom of it, relaxing out of the clash the chemical dictionary's coordinates hold.
//!
//! ```text
//! cargo run --release --example aspirin_relaxes        # numbers, checked
//! ```
//!
//! The wwPDB Chemical Component Dictionary entry for aspirin (`AIN`) gives two sets of
//! coordinates, and the *ideal* ones — generated, not measured — put the acetyl oxygen `O4`
//! 1.645 Å from the ring hydrogen `H1`. Under the Universal Force Field that one pair is most of
//! the molecule's van der Waals energy. This reads the entry, types every atom for UFF, and
//! minimises: L-BFGS with an Armijo line search, one step at a time, until the largest force on
//! any atom is below the tolerance.
//!
//! # What is checked, because a relaxed energy has no closed form
//!
//! Where this force field's minimum is, and what its energy is there, are the force field's
//! answer: nothing outside it says what they should be, so neither is asserted. What *is* exact:
//!
//! - **the molecule is the one the entry says it is**: the atoms read from the atom loop, counted
//!   by element, against `_chem_comp.formula` — a second statement in the same file, written by
//!   the dictionary and not computed from the atoms;
//! - **the energy never rises**, step by step, and a step that moved the atoms lowered it — the
//!   Armijo condition's guarantee, which is the reason this minimiser was chosen over a damped
//!   dynamics that overshoots;
//! - **electrostatics is exactly zero**: the charge model is a later step, every partial charge is
//!   zero until it arrives, and `q_i q_j / r` with a zero in it is zero in floating point too;
//! - **it ends converged**, with the largest force recomputed here from the forces the force field
//!   hands back rather than taken from the minimiser's own report;
//! - **the clash opens**, against a threshold somebody else published: Bondi's van der Waals radii,
//!   1.52 Å for oxygen and 1.20 for hydrogen (J. Phys. Chem. 68, 441 (1964)), less the 0.4 Å
//!   overlap MolProbity calls a serious clash (Word et al., J. Mol. Biol. 285, 1735 (1999)). The
//!   start has to be inside it, or the check could pass on a molecule that never clashed.
//!
//! It does not open to a contact. The relaxed pair is inside the Bondi sum, 2.72 Å: the ester is
//! held planar by its conjugation, and the minimum is the compromise between that and this pair.
//! That is printed and not asserted either way.

use pantometry::forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry::forcefield::uff::KCAL_PER_MOL;
use pantometry::forcefield::{Component, Element, Energy, Molecule, Status};

mod common;
use common::{check_between, heading};

const AIN: &str = include_str!("../../pantometry-forcefield/components/AIN.cif");

/// Bondi's van der Waals radii, Å.
const BONDI_O: f64 = 1.52;
const BONDI_H: f64 = 1.20;
/// The overlap MolProbity calls a serious clash, Å.
const SERIOUS_OVERLAP: f64 = 0.4;

/// One ångström, in metres.
const ANGSTROM: f64 = 1e-10;

/// More than the minimum needs. A run that reaches this without converging fails below.
const MAX_STEPS: usize = 2000;

fn kcal(e: &Energy) -> [(&'static str, f64); 7] {
    let k = |j: f64| j / KCAL_PER_MOL;
    [
        ("bond stretch", k(e.bond)),
        ("angle bend", k(e.angle)),
        ("torsion", k(e.torsion)),
        ("inversion", k(e.inversion)),
        ("van der Waals", k(e.van_der_waals)),
        ("electrostatic", k(e.electrostatic)),
        ("total", k(e.total)),
    ]
}

/// `_chem_comp.formula`, `C9 H8 O4`, as element counts. An element written without a number is
/// one of it.
fn formula_counts(formula: &str) -> Vec<(Element, usize)> {
    let mut out: Vec<(Element, usize)> = formula
        .split_whitespace()
        .map(|word| {
            let digits = word
                .find(|c: char| c.is_ascii_digit())
                .unwrap_or(word.len());
            let (symbol, n) = word.split_at(digits);
            let element = Element::from_symbol(symbol)
                .unwrap_or_else(|| panic!("{symbol} in the formula is not an element"));
            let n = if n.is_empty() {
                1
            } else {
                n.parse().expect("a count")
            };
            (element, n)
        })
        .collect();
    out.sort();
    out
}

fn distance(m: &Molecule, a: &str, b: &str) -> f64 {
    let index = |name: &str| {
        m.component()
            .atoms()
            .iter()
            .position(|atom| atom.name == name)
            .unwrap_or_else(|| panic!("no atom called {name}"))
    };
    let (p, q) = (m.at(index(a)), m.at(index(b)));
    (0..3).map(|k| (p[k] - q[k]).powi(2)).sum::<f64>().sqrt() / ANGSTROM
}

fn main() {
    let component = Component::from_ccd(AIN).expect("the dictionary entry reads");

    heading(&format!(
        "{}: {}",
        component.id(),
        component.name().unwrap_or("(no name)")
    ));
    let stated = formula_counts(component.formula().expect("the entry states a formula"));
    let counted: Vec<(Element, usize)> = component.composition().into_iter().collect();
    println!(
        "  {} atoms, {} bonds; the entry says {}",
        component.atoms().len(),
        component.bonds().len(),
        component.formula().unwrap_or("")
    );
    assert_eq!(
        counted, stated,
        "the atoms read are not the molecule the entry's formula states"
    );
    println!("  counted by element from the atom loop: the same, exactly");

    let mut molecule = Molecule::new("aspirin", component);
    let start = molecule.evaluate().expect("UFF describes aspirin").energy;
    let clash_start = distance(&molecule, "O4", "H1");

    // One step at a time, so the guarantee is checked at every step rather than at the end.
    let tolerance = Molecule::DEFAULT_TOLERANCE;
    let mut last = start.total;
    let mut moved = 0;
    let mut stopped = None;
    for _ in 0..MAX_STEPS {
        let before: Vec<[f64; 3]> = (0..molecule.component().atoms().len())
            .map(|i| molecule.at(i))
            .collect();
        let p = molecule.minimise(1).expect("UFF describes aspirin");
        assert!(
            p.energy <= last,
            "step {}: the energy rose, {:.9} kcal/mol after {:.9}",
            p.steps,
            p.energy / KCAL_PER_MOL,
            last / KCAL_PER_MOL
        );
        let after: Vec<[f64; 3]> = (0..before.len()).map(|i| molecule.at(i)).collect();
        if after != before {
            assert!(
                p.energy < last,
                "step {}: the atoms moved and the energy did not fall",
                p.steps
            );
            moved += 1;
        }
        last = p.energy;
        if p.status != Status::Running {
            stopped = Some(p);
            break;
        }
    }
    let stopped = stopped.expect("the minimiser stopped within the step budget");
    let end = molecule.evaluate().expect("UFF describes aspirin");

    heading("Energy by term, kcal/mol");
    println!("  {:<16} {:>12} {:>12}", "term", "ideal start", "relaxed");
    for ((name, a), (_, b)) in kcal(&start).iter().zip(kcal(&end.energy).iter()) {
        println!("  {name:<16} {a:>12.4} {b:>12.4}");
    }
    println!(
        "  {} steps, {moved} of which moved the atoms; never once uphill",
        stopped.steps
    );
    assert!(moved > 0, "nothing moved");
    assert!(
        start.electrostatic == 0.0 && end.energy.electrostatic == 0.0,
        "no charges were given, and the electrostatic energy is not zero"
    );

    heading("Converged");
    assert_eq!(
        stopped.status,
        Status::Converged,
        "the minimiser stopped without converging"
    );
    let largest = end
        .forces
        .iter()
        .map(|f| (f[0] * f[0] + f[1] * f[1] + f[2] * f[2]).sqrt())
        .fold(0.0f64, f64::max);
    // Four decimals cannot tell 7.7e-5 from the tolerance it is under, so the number first.
    println!(
        "  {:.3e} kcal/mol/Å against a tolerance of {:.0e}",
        largest / KCAL_PER_MOL_ANGSTROM,
        tolerance / KCAL_PER_MOL_ANGSTROM
    );
    check_between(
        "largest force on any atom, recomputed",
        largest / KCAL_PER_MOL_ANGSTROM,
        0.0,
        tolerance / KCAL_PER_MOL_ANGSTROM,
        "kcal/mol/Å",
    );

    heading("The O4-H1 clash");
    let threshold = BONDI_O + BONDI_H - SERIOUS_OVERLAP;
    let clash_end = distance(&molecule, "O4", "H1");
    println!("  a serious clash is closer than {BONDI_O} + {BONDI_H} - {SERIOUS_OVERLAP} = {threshold:.2} Å");
    check_between(
        "O4-H1 at the dictionary's ideal coordinates",
        clash_start,
        0.0,
        threshold,
        "Å",
    );
    check_between("O4-H1 relaxed", clash_end, threshold, f64::INFINITY, "Å");
    let contact = BONDI_O + BONDI_H;
    if clash_end < contact {
        println!(
            "  still {:.3} Å inside the Bondi contact of {contact:.2} Å: opened, and not to a \
             contact",
            contact - clash_end
        );
    } else {
        println!("  past the Bondi contact of {contact:.2} Å as well");
    }
}
