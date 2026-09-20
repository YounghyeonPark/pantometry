//! A small molecule closes an enzyme, and what that does to the enzyme's motions.
//!
//! ```text
//! cargo run --release --example ligand_binding            # numbers, checked
//! cargo run --release --example ligand_binding lid.svg    # and the figure
//! ```
//!
//! Adenylate kinase folds two lids over its substrates and opens again. Both ends of that motion
//! are in the Protein Data Bank — `4AKE` open and empty, `1AKE` closed on **AP5A**, a bi-substrate
//! inhibitor that spans both sites — and the two are 7.1 Å apart, which is a conformational change
//! and not a vibration.
//!
//! The question this example is about is the one a structure alone cannot answer: *what does the
//! molecule do to the protein's dynamics?* An elastic network of the alpha carbons answers it
//! twice over, and the two answers are different questions rather than two goes at one.
//!
//! # A ligand adds nodes. It does not simply add stiffness
//!
//! Put the inhibitor's atoms into the network and the Hessian grows: 57 atoms is 171 more
//! coordinates, and the complex has 807 modes where the protein alone had 636. Its eigenvalues are
//! not the protein's shifted — they are a **different problem's**, whose modes include the ligand
//! moving in its pocket, and reading entry `k` of one against entry `k` of the other compares two
//! motions that are not the same motion.
//!
//! What the protein's own coordinates do under a bound ligand is the Schur complement
//! `H_pp − H_pl H_ll⁺ H_lp`: the bound potential minimised over where the ligand goes. That is
//! `Modes::of_the_protein`, it has 636 modes like the free protein, and entry `k` of it and of the
//! free protein *are* comparable. Both are printed below, side by side, because the difference
//! between them is the point.
//!
//! # What is checked, because a picture of a protein is easy to believe
//!
//! - **Six rigid-body modes**, free and bound and after integrating the ligand out. A body with
//!   nothing holding it has six ways to move without stretching anything, and a model that has
//!   five or seven has a broken Hessian rather than an interesting protein;
//! - **the amplitudes against the eigenvalues**: summing every residue's thermal displacement has
//!   to give `k_BT Σ 1/λ` over the non-rigid modes, because each eigenvector is a unit vector. One
//!   side of that walks the eigenvectors and the other never touches them;
//! - **the direction**, which is the check nothing scalar can make: the softest mode of the *open*
//!   structure against the motion the enzyme is observed to perform. The model never sees the
//!   closed structure. In `3n = 642` dimensions a random direction scores `1/√642 = 0.039`;
//! - **where binding is felt**, against a null: the residues the ligand stiffens most should be the
//!   ones near it, and "near" is measured rather than asserted.
//!
//! # Release only, and why it is worth the seconds
//!
//! Two `642 × 642` eigenproblems, about twelve seconds in release and minutes in debug. The
//! library gate and CI both run examples in release.

use pantometry_protein::{correlation, Modes, Network, Structure};
use pantometry_units::Qty;

mod common;
use common::{check, check_between, heading};

/// Ångström, in metres.
const A: f64 = 1e-10;
/// Body temperature, which is where an enzyme works.
const KELVIN: f64 = 310.15;
/// The cutoff every network in this workspace uses: two residues are joined if their alpha carbons
/// are within it.
const CUTOFF: f64 = 15.0 * A;
/// One newton per metre. The spring constant cancels out of every overlap and every ratio here,
/// and sets the scale of nothing that is compared.
const SPRING: f64 = 1.0;

fn main() {
    let open_text = include_str!("../../pantometry-protein/structures/4AKE.pdb");
    let closed_text = include_str!("../../pantometry-protein/structures/1AKE.pdb");

    // ================================================================ the two photographs
    heading("An enzyme caught open and caught closed, and the molecule in between");

    // **Chain A of each.** Both entries are dimers in the asymmetric unit, and a network built
    // across both joins two molecules that are not joined — every mode after that is about the
    // crystal rather than about the enzyme.
    let open = chain_a(open_text, "4AKE");
    let closed = chain_a(closed_text, "1AKE");
    let ligand = Structure::ligand_from_pdb(closed_text, "AP5", Some('A'));
    println!(
        "  {:<34} {:>6} residues, {} of them in each copy",
        "adenylate kinase",
        open.len(),
        open.len()
    );
    println!(
        "  {:<34} {:>6} atoms  bi-substrate inhibitor, chain A's copy",
        "AP5A",
        ligand.len()
    );
    // The file has 121 `AP5` records: both copies, and seven atoms modelled in two places. The
    // reader takes the first alternate location of each atom, as the protein side does.
    check_between(
        "the ligand is a molecule and not a solvent",
        ligand.len() as f64,
        40.0,
        80.0,
        "atoms",
    );

    // Superposed **closed onto open**, because the network is built in the open structure's frame
    // and a direction has to live there too. The other way round measures the same motion in the
    // wrong coordinates and reads 0.19 where this reads 0.80.
    let fitted = closed.superposed_onto(&open).expect("the same length");
    let apart = open.rmsd_to(&fitted).expect("the same length");
    println!(
        "  {:<34} {:>6.2} A  after superposing, from {:.1} before",
        "the two conformations are apart",
        apart / A,
        open.rmsd_to(&closed).expect("the same length") / A
    );
    check_between(
        "this is a conformational change, not a vibration",
        apart / A,
        5.0,
        10.0,
        "A",
    );

    // ================================================================ the motion, predicted
    heading("The softest mode of the open enzyme, against the motion it is seen to perform");

    let free = Network::new(&open, Qty::from_si(CUTOFF), Qty::from_si(SPRING));
    check(
        "the open enzyme is one connected piece",
        free.components() as f64,
        1.0,
        1e-12,
        "",
    );
    let modes_open = Modes::of(&free);
    check(
        "six ways to move without stretching anything",
        modes_open.rigid_body_modes() as f64,
        6.0,
        1e-12,
        "",
    );

    let went = open
        .direction_to(&fitted)
        .expect("two different structures");
    let random = 1.0 / (3.0 * open.len() as f64).sqrt();
    let mut cumulative = 0.0;
    for k in 0..5 {
        let o = modes_open.overlap_with(k, &went).abs();
        cumulative += o * o;
        println!(
            "  mode {k}: {o:>6.3}   {:>6.3} cumulative",
            cumulative.sqrt()
        );
    }
    println!(
        "  {:<34} {random:>6.3}  in {} dimensions",
        "a direction that knows nothing",
        3 * open.len()
    );
    let lowest = modes_open.overlap_with(0, &went).abs();
    // **The model is never shown the closed structure.** It contributes the direction being asked
    // about and nothing else — not the cutoff, not the spring constant, not which mode to look at.
    check_between(
        "the softest mode is the closing motion",
        lowest / random,
        10.0,
        100.0,
        "x a random direction",
    );

    // ================================================================ what binding does
    heading("A ligand adds nodes, and the protein's own coordinates are what change");

    // The same protein coordinates throughout: the closed structure, with and without its
    // inhibitor. Comparing the open enzyme against the bound one would be comparing two shapes as
    // well as two potentials, and this is about the potential.
    let bare = Network::new(&closed, Qty::from_si(CUTOFF), Qty::from_si(SPRING));
    let bound = Network::new(&closed, Qty::from_si(CUTOFF), Qty::from_si(SPRING))
        .with_ligand(ligand.clone());
    let m_bare = Modes::of(&bare);
    let m_complex = Modes::of(&bound);
    let m_after = Modes::of_the_protein(&bound);

    println!(
        "  {:<34} {:>6} springs free, {} bound",
        "the network gains",
        bare.springs(),
        bound.springs()
    );
    println!(
        "  {:<34} {:>6} free, {} for the complex, {} after integrating the ligand out",
        "modes",
        m_bare.len(),
        m_complex.len(),
        m_after.len()
    );
    // 57 atoms is 171 coordinates, and every one of them is a mode of the complex. That is what
    // "a ligand adds nodes" means, and it is why the complex's spectrum is not the protein's.
    check(
        "the complex has three modes per ligand atom more than the protein",
        (m_complex.len() - m_bare.len()) as f64,
        3.0 * ligand.len() as f64,
        1e-12,
        "modes",
    );
    check(
        "and still six rigid-body modes",
        m_complex.rigid_body_modes() as f64,
        6.0,
        1e-12,
        "",
    );
    check(
        "as does the protein after the ligand is integrated out",
        m_after.rigid_body_modes() as f64,
        6.0,
        1e-12,
        "",
    );

    println!(
        "\n  {:>4}  {:>10}  {:>10}  {:>10}",
        "k", "free", "complex", "after"
    );
    for k in 0..5 {
        println!(
            "  {k:>4}  {:>10.4}  {:>10.4}  {:>10.4}",
            m_bare.eigenvalue(6 + k).to_si(),
            m_complex.eigenvalue(6 + k).to_si(),
            m_after.eigenvalue(6 + k).to_si()
        );
    }
    println!("  the middle column is a different problem's spectrum; the outer two are comparable");

    // ================================================================ where it is felt
    heading("Where the molecule is felt, against where it is");

    let before = m_bare.fluctuations(Qty::from_si(KELVIN));
    let after = m_after.fluctuations(Qty::from_si(KELVIN));
    // **Every residue's share of the thermal motion, summed, against the modes' own.**
    // `<dr^2>_i` is `k_BT` times the sum over modes of `u_ik^2 / lambda_k`; adding it over
    // residues leaves `k_BT * sum 1/lambda_k`, because each eigenvector is a unit vector. One side
    // of that walks the eigenvectors and the other never touches them, so they agree only if the
    // solver's vectors are orthonormal -- which is the thing to ask a solver, and the thing an
    // amplitude cannot fake.
    //
    // **The first version of this check counted modes.** It summed `k_BT` once per mode and
    // compared against `(3n-6) k_BT`, which is `count == count` dressed as equipartition, and it
    // did not read the fluctuations at all. It would have passed over a Hessian whose every
    // eigenvector was wrong.
    //
    // The floor is the summation's: `3n` terms into each residue and `n` residues into the total
    // is about `4n` roundings at the scale of the answer, and `n` is 214.
    for (name, modes, f) in [("free", &m_bare, &before), ("bound", &m_after, &after)] {
        let summed: f64 = f.iter().sum();
        let from_values: f64 = modes
            .spectrum()
            .values()
            .iter()
            .skip(modes.rigid_body_modes())
            .map(|lambda| 1.380_649e-23 * KELVIN / lambda)
            .sum();
        let floor = 4.0 * (3 * closed.len()) as f64 * f64::EPSILON;
        println!(
            "  {:<34} {:>6.3} of a floor of {:.1e}",
            format!("{name}: the two routes agree to"),
            (summed / from_values - 1.0).abs() / floor,
            floor
        );
        check_between(
            &format!("{name}: the amplitudes carry k_BT per mode and no more"),
            (summed / from_values - 1.0).abs() / floor,
            0.0,
            1.0,
            "x the floor",
        );
    }

    let nearest = |i: usize| {
        let p = closed.residues()[i].at;
        ligand
            .iter()
            .map(|l| {
                let d = [p[0] - l[0], p[1] - l[1], p[2] - l[2]];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            })
            .fold(f64::INFINITY, f64::min)
    };
    let mut ranked: Vec<(f64, usize)> = (0..closed.len())
        .map(|i| (after[i] / before[i], i))
        .collect();
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));

    // **Not one residue moves more.** The bound Hessian is the free one plus the ligand's
    // springs, and integrating the ligand out minimises the potential over where it sits; the
    // pseudo-inverse of a larger matrix is smaller, so every mean-square displacement can only
    // fall. That is a theorem rather than an observation, which makes it the right thing to
    // assert over all of them: a sign error in the Schur complement, or `H_pl` where `H_lp`
    // belongs, breaks it on the first residue.
    //
    // It was worth counting. Reading the figure, the blue curve looked higher than the grey one
    // in two places -- and blue is drawn over grey, so the grey that shows is the grey *above*
    // blue and the picture was right. The count is what settled it.
    let loosened = (0..closed.len()).filter(|&i| after[i] > before[i]).count();
    let total: f64 = before.iter().sum::<f64>();
    let bound_total: f64 = after.iter().sum::<f64>();
    println!(
        "  {:<34} {:>6.2} A^2 free, {:.2} bound  x{:.3} over the whole chain",
        "mean-square motion, summed",
        total / (A * A),
        bound_total / (A * A),
        bound_total / total
    );
    check(
        "adding springs makes nothing move more",
        loosened as f64,
        0.0,
        1e-12,
        "residues",
    );

    println!("  {:>8}  {:>7}  {:>9}", "residue", "moves", "from AP5A");
    for (ratio, i) in ranked.iter().take(6) {
        println!(
            "  {:>8}  {:>6.2}x  {:>7.1} A",
            closed.residues()[*i].number,
            ratio,
            nearest(*i) / A
        );
    }
    let ten: f64 = ranked
        .iter()
        .take(10)
        .map(|(_, i)| nearest(*i))
        .sum::<f64>()
        / 10.0;
    let all: f64 = (0..closed.len()).map(nearest).sum::<f64>() / closed.len() as f64;
    println!(
        "  {:<34} {:>6.1} A  against {:.1} A for all {} residues",
        "the ten most stiffened average",
        ten / A,
        all / A,
        closed.len()
    );
    // **The null is the protein itself.** A ranking that had nothing to do with the ligand would
    // average the whole protein's distance; this one has to beat that, and by enough that the
    // difference is not ten residues' worth of noise.
    check_between(
        "binding is felt where the molecule is",
        ten / all,
        0.0,
        0.8,
        "x the mean distance",
    );

    // And the crystallographer's own numbers, which are the only measurement here nothing in this
    // model produced. The bound model is compared against the bound structure's B-factors.
    if let Some(measured) = closed.experimental_fluctuations() {
        let free_r = correlation(&before, &measured).expect("the same length");
        let bound_r = correlation(&after, &measured).expect("the same length");
        println!(
            "  {:<34} {:>6.3} free, {:.3} bound  against the deposited B-factors",
            "correlation", free_r, bound_r
        );
    }

    if let Some(path) = common::output_path() {
        common::write(&path, &draw(&closed, &before, &after, &nearest));
        println!("\n  the lid, the molecule under it, and what stopped moving");
    } else {
        println!("\n  give a filename ending .svg for the figure");
    }
}

/// Chain A of an entry, checked to be the entry it says it is.
fn chain_a(text: &str, expect: &str) -> Structure {
    let header = text.lines().next().expect("a HEADER line");
    assert!(
        header.starts_with("HEADER") && header.len() >= 66 && &header[62..66] == expect,
        "this file's HEADER says {:?}, not {expect}",
        header.get(62..66)
    );
    Structure::from_pdb(text)
        .expect("alpha carbons")
        .chain('A')
        .expect("a chain A")
}

/// The enzyme along its sequence: how much each residue moves, free and bound, with the ligand's
/// reach marked.
fn draw(
    protein: &Structure,
    before: &[f64],
    after: &[f64],
    nearest: &dyn Fn(usize) -> f64,
) -> String {
    use common::svg::{document, rgb, ticks, Plot};
    let (w, h) = (880.0, 420.0);
    let n = protein.len();
    let hi = before
        .iter()
        .chain(after.iter())
        .fold(0.0f64, |m, v| m.max(v.sqrt()));
    let first = protein.residues()[0].number as f64;
    let last = protein.residues()[n - 1].number as f64;
    let mut plot =
        Plot::new(w, h, (first, last), (0.0, hi / A * 1.05)).viewport(64.0, 54.0, 800.0, 300.0);

    // Where the ligand reaches, as a band along the sequence. Drawn first, so the curves are over
    // it: it is the question and they are the answer.
    for i in 0..n {
        if nearest(i) < 8.0 * A {
            let x = protein.residues()[i].number as f64;
            plot.cell(
                (x - 0.5, x + 0.5),
                (0.0, hi / A * 1.05),
                &rgb(233, 237, 243),
            );
        }
    }
    for (values, colour, width) in [
        (before, rgb(120, 122, 132), 1.4),
        (after, rgb(74, 132, 238), 1.8),
    ] {
        plot.polyline(
            (0..n).map(|i| (protein.residues()[i].number as f64, values[i].sqrt() / A)),
            &colour,
            width,
        );
    }
    plot.axes(
        &ticks(first, last, 8),
        &ticks(0.0, hi / A * 1.05, 5),
        |v| format!("{v:.0}"),
        |v| format!("{v:.1}"),
    );
    plot.title("what a bound molecule stops from moving");
    plot.caption("residue number against root-mean-square displacement (A) at 310 K");
    plot.footnote(&format!(
        "grey: the enzyme alone. blue: the same coordinates with AP5A bound, the ligand integrated out. the band is the {} residues within 8 A of it",
        (0..n).filter(|&i| nearest(i) < 8.0 * A).count()
    ));
    document(w, h, [plot.into_body()])
}
