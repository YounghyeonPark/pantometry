//! A small molecule closes an enzyme, and what that does to the enzyme's motions.
//!
//! ```text
//! cargo run --release --example ligand_binding             # numbers, checked
//! cargo run --release --example ligand_binding lid.svg     # and the figure
//! cargo run --release --example ligand_binding close.json  # the closing, as a solid
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
//!   ones near it, and "near" is measured rather than asserted;
//! - **the solid itself**, because the animation is a mesh and a mesh is the one thing a picture
//!   cannot check. A hole reads as shadow, an inside-out surface lights like any other, and a tube
//!   that folds through itself renders as a bead. So the sweep is measured against the exact area
//!   and volume of the shape it is — a regular prism, an icosahedron — the winding against the
//!   *sign* of the enclosed volume, the rings against the radius they are supposed to keep, and
//!   the path against the curve it follows.
//!
//! # The picture is a Ca trace, and says so
//!
//! `Structure` holds one atom per residue, so what the tube follows is the alpha-carbon chain and
//! not a molecular surface: there are no side chains in this model to draw. That is the honest
//! limit of a Ca network and the picture does not pretend otherwise.
//!
//! The **ligand** is not like that. It has real atom positions and, now that `ligand_from_pdb`
//! carries the element column, real radii: 20 carbons, 10 nitrogens, 22 oxygens and 5 phosphorus,
//! drawn at Bondi's 1.52 to 1.80 Å. It was 57 identical balls, which drew five phosphorus atoms
//! the size of a carbon — a picture asserting something about the molecule that the data does not.
//!
//! # Release only, and why it is worth the seconds
//!
//! Two `642 × 642` eigenproblems, about twelve seconds in release and minutes in debug. The
//! library gate and CI both run examples in release.

use pantometry::scene::Frame;
use pantometry_core::Reading;
use pantometry_protein::{correlation, Atom, Modes, Network, Structure};
use pantometry_units::Qty;

mod common;
use common::{check, check_between, heading, mesh};

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
/// The radius the backbone is drawn at.
///
/// **Set by the chain's own curvature, not by taste.** A swept tube of radius r following a curve
/// of radius R reaches its own axis when r reaches R, and this backbone turns through six places
/// where R is about 0.9 A -- residues 10, 24, 41, 86, 129 and 149, against a median of 3.9 A over
/// the whole chain. At 1.6 A the tube folded through itself at all six and rendered as a bead
/// each time, with every other check still passing. The run measures the ratio and refuses above
/// one, on every frame of the animation and not only on the structure at rest.
const TUBE_RADIUS: f64 = 0.6 * A;
/// How many sides the backbone's cross-section has.
///
/// Six rather than sixteen because the wire format writes every panel out on every one of 48
/// frames. At this radius the cross-section is a few pixels across and the facets do not read.
///
/// # Where the 18.2 MB actually goes, measured
///
/// "The positions are where the file size goes" stood here and is only half true:
///
/// | | of the file | changes between frames? |
/// | --- | --- | --- |
/// | backbone triangles | 38.8% | **no** |
/// | backbone positions | 37.3% | yes |
/// | backbone values | 8.9% | yes |
/// | the whole molecule | 11.5% | **no** |
///
/// **Half the file is bytes repeated unchanged forty-eight times.** A triangle list is topology
/// and never moves; the ligand is what the protein closes *on* and does not move either. A format
/// where a panel could omit an array and mean "as the frame before" would cut 18.2 MB to about
/// 9.3 — and it is deliberately not that, because every reader would become stateful and a frame
/// that simply forgot an array would render as an empty mesh rather than as an error. That is a
/// silent failure bought for a factor of two on a local artefact nothing ships, and halving
/// `FRAMES` buys the same factor for nothing. The number is written down here so the next person
/// weighing it does not have to measure it again.
const TUBE_SIDES: usize = 6;
/// How many samples of the spline each 3.8 A step becomes before the tube is swept along it.
///
/// **This is the number that decides whether the chain reads as a chain.** A straight polyline
/// through alpha carbons turns as much as 60 degrees in one step, and the miter widens the ring
/// by one over cos of the half-angle exactly there -- so the picture came out as a row of flat
/// lozenges. The spline still passes through every measured atom; it only turns less between
/// them, which the run checks to within a rounding.
const TUBE_SMOOTH: usize = 3;
/// The radius each element is drawn at: Bondi's van der Waals radii, in metres.
///
/// **A drawing convention, which is why it is here and not in `pantometry-protein`.** A van der
/// Waals radius is a property of an element and not a physics this workspace models -- nothing in
/// the crate computes with it, and putting a periodic table in a crate that reads files would be
/// the wrong crate growing the wrong thing.
///
/// This used to be one radius for all 57 atoms, because `ligand_from_pdb` handed back bare
/// coordinates: five phosphorus atoms drawn the size of a carbon, which is a picture asserting
/// something about the molecule that the data does not. The reader carries the element column
/// now, and `the_ligand_is_the_molecule_its_formula_says` holds it against AP5A's own formula.
const BONDI: [(&str, f64); 6] = [
    ("C", 1.70 * A),
    ("N", 1.55 * A),
    ("O", 1.52 * A),
    ("P", 1.80 * A),
    ("S", 1.80 * A),
    ("F", 1.47 * A),
];
/// What an element with no entry in [`BONDI`] is drawn at, and the run says when it uses this.
///
/// Carbon's, because a ligand is mostly carbon -- but a picture that silently drew an iron the
/// size of a carbon would be the defect this change removed, coming back through the table
/// instead of through the reader.
const UNLISTED_RADIUS: f64 = 1.70 * A;
/// How many times each ligand sphere is subdivided: 20 * 4^level triangles per atom.
///
/// Zero -- the icosahedron itself. The molecule does not move and its 57 balls are a few pixels
/// across, but the wire format writes them out on all 48 frames regardless, so one subdivision
/// would be 1710 extra points per frame for a roundness nothing can see. This is the same trade
/// as `TUBE_SIDES` and it is decided by the same fact about the format.
const ATOM_LEVEL: u32 = 0;

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
    let ligand_at: Vec<[f64; 3]> = ligand.iter().map(|a| a.at).collect();
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
    // **What it is made of, and what that makes it look like.** The reader carries the element
    // column now, so the molecule is drawn at Bondi radii rather than as 57 identical balls --
    // and the test `the_ligand_is_the_molecule_its_formula_says` holds this tally against AP5A's
    // own formula, `C20H29N10O22P5` minus the hydrogens a 1.9 A structure does not have.
    let mut tally: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
    for a in &ligand {
        *tally.entry(a.element.as_str()).or_default() += 1;
    }
    let formula: String = tally
        .iter()
        .map(|(e, n)| format!("{e}{n}"))
        .collect::<Vec<_>>()
        .join(" ");
    let unlisted = ligand.iter().filter(|a| !bondi(&a.element).1).count();
    println!(
        "  {:<34} {:>6}  drawn at {:.2} to {:.2} A, {} not in the table",
        "by element",
        formula,
        ligand
            .iter()
            .map(|a| bondi(&a.element).0 / A)
            .fold(f64::MAX, f64::min),
        ligand
            .iter()
            .map(|a| bondi(&a.element).0 / A)
            .fold(0.0f64, f64::max),
        unlisted
    );
    check_between(
        "every atom is drawn at its own element's radius",
        unlisted as f64,
        0.0,
        0.0,
        "fell back to carbon's",
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

    // ================================================================ how far one mode gets
    heading("Following that mode, and where it stops being the motion");

    // **A mode is a direction, and a conformational change is a path.** Walking the open structure
    // along mode 0 carries it towards the closed one -- and only so far: a harmonic approximation
    // is a straight line through a curve, and past the tangent point it leaves again. Nothing here
    // is fitted; the mode, the amplitude and the closed structure are each computed once.
    //
    // The amplitude is the thermal one, `sqrt(2 k_BT / lambda)`, so "x16" is a real unit rather
    // than a knob: sixteen times the excursion this mode actually has at body temperature.
    let shape = modes_open.shape(0);
    let amplitude = (2.0 * 1.380_649e-23 * KELVIN / modes_open.eigenvalue(0).to_si()).sqrt();
    println!(
        "  {:<34} {:>6.2} A  thermal excursion of the softest mode at {KELVIN:.0} K",
        "amplitude",
        amplitude / A
    );
    let walked = |steps: f64| -> f64 {
        let moved: Vec<[f64; 3]> = (0..open.len())
            .map(|i| {
                let p = open.residues()[i].at;
                [
                    p[0] + steps * amplitude * shape[3 * i],
                    p[1] + steps * amplitude * shape[3 * i + 1],
                    p[2] + steps * amplitude * shape[3 * i + 2],
                ]
            })
            .collect();
        // Superposed before measuring, so this is shape and not where the mode translated it.
        let s = Structure::from_positions(moved);
        s.superposed_onto(&fitted)
            .expect("the same length")
            .rmsd_to(&fitted)
            .expect("the same length")
    };
    let at_rest = walked(0.0);
    let mut best = (0.0f64, at_rest);
    for k in 0..=40 {
        let steps = k as f64;
        let d = walked(steps);
        if d < best.1 {
            best = (steps, d);
        }
    }
    for steps in [-8.0, 0.0, best.0, 2.0 * best.0] {
        println!(
            "  {:>6.0} x amplitude   {:>6.1} A of excursion   {:>5.2} A from the closed form",
            steps,
            (steps.abs() * amplitude) / A,
            walked(steps) / A
        );
    }
    // **The sign is the model's, not a fitted one.** Going the other way along the same mode takes
    // the enzyme further from the closed structure, which is what makes "towards" mean something.
    check_between(
        "the mode leads towards the closed form and not away",
        walked(-best.0) / at_rest,
        1.0,
        2.0,
        "x the distance at rest",
    );
    check_between(
        "and gets a third of the way there before it overshoots",
        1.0 - best.1 / at_rest,
        0.25,
        0.6,
        "of the way",
    );

    // ================================================================ what binding does
    heading("A ligand adds nodes, and the protein's own coordinates are what change");

    // The same protein coordinates throughout: the closed structure, with and without its
    // inhibitor. Comparing the open enzyme against the bound one would be comparing two shapes as
    // well as two potentials, and this is about the potential.
    let bare = Network::new(&closed, Qty::from_si(CUTOFF), Qty::from_si(SPRING));
    let bound = Network::new(&closed, Qty::from_si(CUTOFF), Qty::from_si(SPRING))
        .with_ligand(ligand_at.clone());
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
        ligand_at
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

    // ======================================================= the shapes the picture is made of
    heading("The solid, against the closed forms for the shapes it is made of");

    // **Checked against the exact area and volume of the shape it *is*,** not against a second
    // sweep. `PanelData::surface` says in its own documentation that the winding is the caller's
    // to get right and is not checked there, and an inside-out solid renders as a surface either
    // way -- so the *signed* volume is the only number here that can tell the two apart.
    //
    // A tube along a straight path is a right prism on a regular n-gon. Perimeter `2 n r
    // sin(pi/n)`, cross-section `(n/2) r^2 sin(2 pi/n)`: elementary, and exact rather than a
    // limit. Read again as n doubles, the same pair says how fast the polygon becomes the circle.
    let straight = [[0.0, 0.0, 0.0], [0.0, 0.0, 10.0 * A]];
    let round_area = 2.0 * std::f64::consts::PI * 2.0 * 10.0 + 2.0 * std::f64::consts::PI * 4.0;
    let mut shortfall = Vec::new();
    for sides in [8usize, 16, 32] {
        let solid = mesh::tube(&straight, &[0.0, 0.0], 2.0 * A, sides);
        let (area, volume) = mesh::straight_tube_exactly(2.0 * A, 10.0 * A, sides);
        // N terms, one rounding of order eps each. The rigorous bound for a recursive sum of N
        // positive terms is (N-1)*eps/2 and the realistic behaviour is a random walk at
        // sqrt(N)*eps; measured here at 0.007 to 0.019 of this floor, so the 4x is slack and not
        // a fit. **What it does not carry is cancellation**, and that is a property of these
        // meshes rather than of the model: every solid checked against a closed form is built at
        // the origin, where the volume's terms are all one sign and the measured condition number
        // is 1.000. A mesh checked far from the origin would need the factor
        // max|coordinate| / shortest edge on top -- 2.6e4 for this same straight tube at 1e4 A,
        // where the volume check fails on a geometrically perfect mesh.
        let floor = 4.0 * solid.faces.len() as f64 * f64::EPSILON;
        check(
            &format!("a {sides}-sided tube's area"),
            solid.area() / (A * A),
            area / (A * A),
            floor,
            "A^2",
        );
        check(
            &format!("a {sides}-sided tube's volume"),
            solid.volume() / (A * A * A),
            volume / (A * A * A),
            floor,
            "A^3",
        );
        check_between(
            &format!("a {sides}-sided tube is closed"),
            solid.unshared_edges() as f64,
            0.0,
            0.0,
            "open edges",
        );
        // **From the mesh, not from the formula.** This subtracted `straight_tube_exactly`'s
        // area from a literal and called the result a check on the sweep: no mesh entered it, and
        // a `tube` inflated 10% printed 3.9466 and 3.9865 bit for bit.
        shortfall.push(round_area - solid.area() / (A * A));
    }
    // The n-gon's shortfall against the circle goes as 1/n^2, so each doubling quarters it. The
    // shortfall is measured off the swept mesh, so a sweep built at the wrong radius has the
    // wrong shortfall and loses the rate -- and `round_area` restating the radius and the length
    // as literals is what makes that true rather than a coincidence of two formulas agreeing.
    for w in shortfall.windows(2) {
        check_between(
            "and approaches the cylinder as 1/n^2",
            w[0] / w[1],
            3.7,
            4.3,
            "x per doubling",
        );
    }

    // The sphere begins as an icosahedron, whose area and volume are elementary too, and every
    // subdivision quarters what it is short of the sphere it is being pushed onto.
    //
    // At phosphorus's radius, because these check the *mesh* and any radius would do -- and if
    // one has to be named, the largest of the ones actually drawn is the one whose triangles are
    // biggest and whose roundings are worst.
    let ball_radius = bondi("P").0;
    let ball = mesh::sphere([0.0, 0.0, 0.0], ball_radius, 0, 0.0);
    let (ico_area, ico_volume) = mesh::icosahedron_exactly(ball_radius);
    let floor = 4.0 * ball.faces.len() as f64 * f64::EPSILON;
    check(
        "the icosahedron's area",
        ball.area() / (A * A),
        ico_area / (A * A),
        floor,
        "A^2",
    );
    check(
        "the icosahedron's volume",
        ball.volume() / (A * A * A),
        ico_volume / (A * A * A),
        floor,
        "A^3",
    );
    check_between(
        "the icosahedron is closed",
        ball.unshared_edges() as f64,
        0.0,
        0.0,
        "open edges",
    );
    // **And `Mesh::append`, here rather than only in the animation.** It renumbers one mesh's
    // triangles onto the end of another, and the molecule is 57 spheres joined that way -- but
    // `closing` is the only caller and CI runs this example without an output path, so a version
    // that forgot to renumber shipped 57 spheres all indexing the first one's twelve vertices
    // and exited 0 twice.
    //
    // **The two checks catch different failures, and it is not the one I first wrote down.** A
    // missing renumber does *not* collapse the volume: both face sets then enclose the first ball
    // twice over, which is exactly the 2x expected, and the volume check passes. What sees it is
    // the closed count -- every edge belongs to four triangles instead of two, measured at 30
    // open edges. The volume is what sees an offset that is wrong rather than absent, where the
    // triangles index real but wrong vertices and the geometry becomes garbage that is still
    // closed. Neither is redundant and neither covers the other.
    let mut pair = ball.clone();
    pair.append(&mesh::sphere(
        [4.0 * ball_radius, 0.0, 0.0],
        ball_radius,
        0,
        0.0,
    ));
    check(
        "two balls appended hold two balls",
        pair.volume() / (A * A * A),
        2.0 * ico_volume / (A * A * A),
        2.0 * floor,
        "A^3",
    );
    check_between(
        "and the pair is still closed",
        pair.unshared_edges() as f64,
        0.0,
        0.0,
        "open edges",
    );

    let sphere_area = 4.0 * std::f64::consts::PI * (ball_radius / A).powi(2);
    // The quarter-per-level is asymptotic in the subdivision's edge, and level 0's edge is about
    // as long as the radius itself -- nowhere near that regime. What pins level 0 is its own exact
    // area and volume, checked directly above; the rate is asked from level 1, where the edge
    // halves each time and nothing else changes.
    //
    // **The excluded ratio is 3.3226, and it is written here so nobody widens the range to admit
    // it.** The retained ones run 3.8070, 3.9501, 3.9874 and march to 4 rather than sitting
    // anywhere; the first is the tightest point in this section, 0.107 clear of the bound.
    let missing: Vec<f64> = (1..5)
        .map(|level| {
            sphere_area - mesh::sphere([0.0, 0.0, 0.0], ball_radius, level, 0.0).area() / (A * A)
        })
        .collect();
    for w in missing.windows(2) {
        check_between(
            "and approaches the sphere as 1/4 a level",
            w[0] / w[1],
            3.7,
            4.3,
            "x per level",
        );
    }

    // And the thing actually drawn. The backbone at rest, as a solid.
    //
    // **On the smoothed trace, because that is the path the sweep is given.** Every one of these
    // checks would pass on the raw alpha-carbon polyline and say nothing about the mesh the
    // animation writes -- the crowding in particular, which is the one that changes most when the
    // steps get shorter.
    let raw: Vec<[f64; 3]> = closed.residues().iter().map(|r| r.at).collect();
    let (trace, _) = mesh::smoothed(&raw, &vec![0.0; raw.len()], TUBE_SMOOTH);
    // A Catmull-Rom passes through its control points -- that is what makes it an interpolating
    // spline rather than an approximating one, and it is the whole of this picture's claim to
    // still be showing where the atoms are. Exactly, not nearly: sample `k * TUBE_SMOOTH` is
    // control point `k` by construction, at `t = 0`, so any deviation is a bug and not a rounding.
    let strayed = (0..raw.len())
        .map(|k| {
            let p = trace[k * TUBE_SMOOTH];
            let d = [p[0] - raw[k][0], p[1] - raw[k][1], p[2] - raw[k][2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .fold(0.0f64, f64::max);
    // Not exactly zero, and the reason is worth the line. The spline is evaluated by Barry and
    // Goldman's pyramid of lerps. At a control point the weights come out exactly 1 and 0 at
    // every step but one: the second-to-last blends a point with *itself* at a weight strictly
    // between, which is that point in arithmetic and that point plus three roundings in floating
    // point. Worked through, that bounds the 3D deviation by 2.6 eps of a coordinate; the floor
    // of 4 is above it with a cushion, and nothing else in the evaluation is lossy.
    // **The quantity every absolute floor in this section is made of.** Both the spline's
    // deviation and the ring's radius error are roundings at the ulp of a coordinate, so they
    // scale with how far the structure sits from the origin -- 5.7e-9 m for this one -- and not
    // with how big it is.
    let magnitude = trace.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
    check_between(
        "the spline still passes through every atom",
        strayed / (4.0 * f64::EPSILON * magnitude),
        0.0,
        1.0,
        "x a floor of 4 roundings",
    );
    println!(
        "  {:<44} {:>12} {:<8} from {} alpha carbons",
        "the trace the sweep is given",
        trace.len(),
        "points",
        raw.len()
    );
    let backbone = mesh::tube(&trace, &vec![0.0; trace.len()], TUBE_RADIUS, TUBE_SIDES);
    let contour: f64 = trace
        .windows(2)
        .map(|w| {
            let d = [w[1][0] - w[0][0], w[1][1] - w[0][1], w[1][2] - w[0][2]];
            (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
        })
        .sum();
    println!(
        "  {:<44} {:>12} {:<8} over {:.0} A of chain",
        "the backbone tube",
        backbone.faces.len(),
        "triangles",
        contour / A
    );
    // **A tube that folded through itself would pass every check above.** It is closed, its
    // rings keep their radius, and it renders as a bead where the chain doubles back. What says
    // it has not is the curvature: a tube of radius r following a curve of radius R closes on the
    // inside of the bend and reaches its own axis when r reaches R.
    // Printed as the curve itself rather than as the ratio on the raw polyline, because that
    // number is not comparable: it is the same estimator on a 3.8 A sampling, and three points
    // that far apart cannot resolve a radius of 1 A. Coarser sampling reporting a gentler curve
    // is the grid talking, not the chain.
    println!(
        "  {:<44} {:>12.3} {:<8} against a median of {:.2} A along the chain",
        "the tightest curve the chain takes",
        TUBE_RADIUS / mesh::crowding(&trace, TUBE_RADIUS) / A,
        "A",
        median_curve(&trace) / A
    );
    check_between(
        "the tube against the tightest curve it follows",
        mesh::crowding(&trace, TUBE_RADIUS),
        0.0,
        1.0,
        "x the radius of curvature",
    );
    check_between(
        "the backbone tube is closed",
        backbone.unshared_edges() as f64,
        0.0,
        0.0,
        "open edges",
    );
    // **Constant radius about its own segment is what makes a tube a tube**, and it is the one
    // property of a mitered sweep no picture can show: a tube pinched at every turn renders as a
    // smooth solid and still counts zero open edges. Every ring vertex has to sit exactly
    // `TUBE_RADIUS` from the line of the segment arriving at its point and from the line of the
    // one leaving. Exact, so the floor is the arithmetic and nothing else -- a few roundings on
    // numbers of order the protein's size.
    let ring_error = |path: &[[f64; 3]]| {
        let solid = mesh::tube(path, &vec![0.0; path.len()], TUBE_RADIUS, TUBE_SIDES);
        let scale = path.iter().flatten().fold(0.0f64, |m, c| m.max(c.abs()));
        mesh::worst_radius_error(path, &solid, TUBE_RADIUS, TUBE_SIDES)
            / (8.0 * f64::EPSILON * scale)
    };
    let here = ring_error(&trace);
    check_between(
        "every ring keeps the tube's radius",
        here,
        0.0,
        1.0,
        "x a floor of 8 roundings",
    );
    // **And the floor is divided by the quantity the error is actually made of.** It was divided
    // by the box extent, which is translation-invariant while the error is not: a rounding at the
    // ulp of a coordinate grows with the distance from the origin, and this structure's happens
    // to sit close enough that nothing failed. Move the same tube 1000 A -- the same shape, bit
    // for bit, plus a shift -- and against the box extent it reads 18x worse and fails the line
    // above, while against the coordinate magnitude it does not move at all.
    let moved: Vec<[f64; 3]> = trace
        .iter()
        .map(|p| [p[0] + 1000.0 * A, p[1], p[2]])
        .collect();
    check_between(
        "and the floor is tied to the right quantity",
        ring_error(&moved) / here,
        0.2,
        5.0,
        "x, a thousand angstroms away",
    );
    // **And the winding of the mesh that actually ships.** This was a `println!`, on the
    // grounds that whether a bend conserves the straight tube's volume exactly is a claim I had
    // not proved -- and it is not, the ratio is 0.9983. But the *sign* needs no theorem, and
    // leaving it unasserted left the whole many-ring path uncovered: reversing every face of any
    // tube past two rings flipped this number to -0.9983 and changed nothing else. Both runs
    // exited 0, 48 inside-out frames were written, and neither the closed-edge count nor the
    // radius invariant nor the crowding could see it -- a reversed triangle still leaves its
    // edges shared by two. The range is wide because it is not a theorem; it is wide in the one
    // direction that matters not at all.
    let (_, straight_volume) = mesh::straight_tube_exactly(TUBE_RADIUS, contour, TUBE_SIDES);
    check_between(
        "the backbone tube is right-side-out",
        backbone.volume() / straight_volume,
        0.9,
        1.1,
        "x the straight tube of that length",
    );

    // **And the conformations the animation actually draws.** Every check above is aimed at the
    // resting `closed` structure, and not one of the 48 frames is that structure: they are the
    // open one displaced along its softest mode, and measured over all 48 they crowd at 0.698
    // against the 0.671 under test here. `closing` asserts each frame as it builds it, but CI
    // runs this example with no output path and therefore never calls it -- so the frames had
    // nothing covering them.
    //
    // They can be checked from here for free. Crowding is a property of *shape*, and `closing`
    // builds its frames in the closed structure's frame, which is a rigid motion of what
    // `modes_open` gives here -- a rotated network has rotated modes. So the same excursion
    // schedule, walked in both directions, is the same set of shapes, without the second
    // eigenproblem `closing` pays for.
    let worst = (0..FRAMES)
        .flat_map(|f| [-1.0, 1.0].map(|sign| sign * excursion(f)))
        .map(|steps| {
            let displaced: Vec<[f64; 3]> = (0..open.len())
                .map(|i| {
                    let p = open.residues()[i].at;
                    [
                        p[0] + steps * amplitude * shape[3 * i],
                        p[1] + steps * amplitude * shape[3 * i + 1],
                        p[2] + steps * amplitude * shape[3 * i + 2],
                    ]
                })
                .collect();
            let (curve, _) = mesh::smoothed(&displaced, &vec![0.0; displaced.len()], TUBE_SMOOTH);
            mesh::crowding(&curve, TUBE_RADIUS)
        })
        .fold(0.0f64, f64::max);
    check_between(
        "nor on any conformation the animation draws",
        worst,
        0.0,
        1.0,
        "x the radius of curvature",
    );

    // **And the one crowding is blind to.** Crowding is local: it asks whether the tube folds
    // where the path bends. Two *straight* stretches far apart along the chain -- a helix packed
    // against a sheet, which is what a protein is made of -- bend nowhere and can still pass
    // within a tube's width of each other, and a tube that spanned that gap would fuse them into
    // one body that is closed, right-side-out, keeps its radius and does not fold.
    //
    // Excluding samples within `APART` of each other excludes 10 A of path, which is eight times
    // the tube's diameter and well inside what crowding already covers.
    const APART: usize = 8;
    let approach = mesh::nearest_approach(&trace, APART);
    println!(
        "  {:<44} {:>12.2} {:<8} between stretches {APART} samples apart or more",
        "the chain's closest approach to itself",
        approach / A,
        "A"
    );
    check_between(
        "and two stretches of chain do not fuse",
        2.0 * TUBE_RADIUS / approach,
        0.0,
        1.0,
        "x the gap",
    );

    // ================================================================ the caption a machine reads
    //
    // **`docs/protein-app.png` is a photograph of a GPU.** CI has no adapter, so it is the third
    // figure in this repository nothing in CI can refresh -- `docs/bench-app.png` and
    // `docs/editor.png` are the other two, and `docs/README.md` spends a paragraph on why that
    // matters: a picture nothing compares ages in silence, and every change to the thing it shows
    // leaves it a little more wrong.
    //
    // So the geometry the picture is of is written down beside it, and compared here, on the
    // argument-less invocation CI makes of every example on every commit. Every knob that decides
    // what the frame looks like is in it -- the radius, the sides, the smoothing, the atom radius
    // and its subdivision, how many frames there are and how far out they go -- so moving one
    // turns this example red instead of leaving the picture quietly wrong.
    //
    // **It sees no pixels**, which is the same trade the other two captions make. A change to the
    // shading, the colour scale or the camera moves nothing here. What it holds is the failure
    // that was actually coming.
    let caption = format!(
        concat!(
            "residues            {:>10} alpha carbons, {} ligand atoms\n",
            "trace               {:>10} points, spline x{} through every one\n",
            "tube                {:>10.2} A radius, {} sides, {} triangles\n",
            "molecule            {:>10} at Bondi radii, level {}, {} triangles\n",
            "tightest curve      {:>10.3} A, tube at {:.4} of folding\n",
            "closest approach    {:>10.2} A, tube {:.4} of the gap\n",
            "animation           {:>10} frames to {:.0} amplitudes; drawn at {} = {:.2}\n",
        ),
        closed.len(),
        ligand.len(),
        trace.len(),
        TUBE_SMOOTH,
        TUBE_RADIUS / A,
        TUBE_SIDES,
        backbone.faces.len(),
        formula,
        ATOM_LEVEL,
        ligand.len()
            * mesh::sphere([0.0; 3], ball_radius, ATOM_LEVEL, 0.0)
                .faces
                .len(),
        TUBE_RADIUS / mesh::crowding(&trace, TUBE_RADIUS) / A,
        mesh::crowding(&trace, TUBE_RADIUS),
        approach / A,
        2.0 * TUBE_RADIUS / approach,
        FRAMES,
        REACH,
        FIGURE_FRAME,
        excursion(FIGURE_FRAME),
    );
    // **Writing the caption comes before checking it, and producing anything skips the check.**
    // The recipe in the message below is retake, then write this file — and with a stale caption
    // the *retake* died right here, so the one command that fixes a stale figure was the command
    // the staleness blocked. Only the `.txt` arm was guarded against that, which is half of it: a
    // run producing any artefact is a run refreshing the figure, and the compare belongs to the
    // argument-less run, which is the one CI makes on every commit.
    if let Some(path) = common::output_path() {
        if path.ends_with(".txt") {
            common::write(&path, &caption);
            return;
        }
    }
    if common::output_path().is_none() {
        compare_caption(&caption);
    }

    match common::output_path() {
        Some(path) if path.ends_with(".svg") => {
            common::write(&path, &draw(&closed, &before, &after, &nearest));
            println!("\n  the lid, the molecule under it, and what stopped moving");
        }
        Some(path) => {
            let frames = closing(&open, &closed, &ligand);
            let asset = if path.ends_with(".json") {
                pantometry::view::to_json("adenylate kinase closing", &frames)
            } else {
                pantometry::view::html("adenylate kinase closing", &frames)
            };
            common::write(&path, &asset);
            println!(
                "\n  {} frames. `pantometry view` plays it with space, or open the .html and drag",
                frames.len()
            );
        }
        None => println!(
            "\n  .svg draws the figure, .json or .html the animation, .txt the figure's caption"
        ),
    }
}

/// The geometry `docs/protein-app.png` is a picture of, against what is stored beside it.
///
/// A *missing* file is a skip — a checkout without `docs/` is not this example's business, and
/// neither is a packaged crate. A file that disagrees is a failure.
fn compare_caption(caption: &str) {
    let beside = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/protein-app.txt")
        .canonicalize()
        .ok();
    match beside.as_deref().map(std::fs::read_to_string) {
        Some(Ok(stored)) => {
            // Carriage returns are the checkout's, not the geometry's: `core.autocrlf` rewrites
            // the stored file and this string has none.
            let flat = |t: &str| t.replace('\r', "");
            assert_eq!(
                flat(&stored),
                flat(caption),
                "the solid has changed since `docs/protein-app.png` was taken. Retake it --\n  \
                 cargo run --release --example ligand_binding closing.json\n  \
                 cd app && cargo run --release -- view ../closing.json --frame {FIGURE_FRAME} \
                 --snapshot ../docs/protein-app.png\n\
                 -- and write this file with `--example ligand_binding docs/protein-app.txt`. \
                 Editing the text alone would restore the green and leave the picture as stale \
                 as it is."
            );
            println!("  {:<44} {:>12}", "the app figure's caption agrees", "yes");
        }
        // A checkout without `docs/` is not this example's business, and neither is a packaged
        // crate. A *missing* file is a skip; a file that disagrees is a failure.
        _ => println!(
            "  {:<44} {:>12}",
            "no docs/protein-app.txt beside this", "skipped"
        ),
    }
}

/// How many frames the closing animation holds: out along the mode and back.
const FRAMES: usize = 48;
/// Which frame `docs/protein-app.png` is of.
///
/// Halfway round the cosine, which is the closest the walk comes to the closed structure and the
/// frame where both lids are furthest from where they started -- so it is the one frame that
/// carries the whole run in a still. Named here rather than written into the caption and the
/// retake command separately, because those two disagreeing is how a figure of one frame comes to
/// claim it is of another.
const FIGURE_FRAME: usize = FRAMES / 2;
/// How far out it goes, in thermal amplitudes of the softest mode.
///
/// Where the walk in `main` measures its closest approach to the closed structure. Past it the
/// straight line leaves the curve again, which the animation shows by going there and coming back.
const REACH: f64 = 16.0;

/// The enzyme closing along its own softest mode, frame by frame, over the molecule it closes on.
///
/// **In the closed structure's frame.** The ligand's coordinates are in that frame and nothing in
/// this crate hands out the transform that would carry them into the open one, so the network is
/// rebuilt from the open structure superposed onto the closed: a rotated network has rotated
/// modes, and the picture then has the protein and the molecule in one place.
///
/// The amplitude is not thermal. Sixteen times the excursion this mode actually has at body
/// temperature is a path being followed, not a motion being watched, and the caption says so --
/// what is physical about it is the *direction*, which is the model's and is checked in `main`.
/// Where along the mode frame `frame` sits, in thermal amplitudes: out and back as a cosine.
///
/// Shared with `main`, which walks the same schedule to check that no conformation this draws
/// folds its tube through itself. Two copies of the formula would drift, and the one in `main`
/// would then be checking shapes nothing renders.
fn excursion(frame: usize) -> f64 {
    let phase = std::f64::consts::TAU * frame as f64 / FRAMES as f64;
    REACH * (1.0 - phase.cos()) / 2.0
}

fn closing(open: &Structure, closed: &Structure, ligand: &[Atom]) -> Vec<Frame> {
    use pantometry::scene::{Panel, PanelData, Placed};

    let start = open.superposed_onto(closed).expect("the same length");
    let network = Network::new(&start, Qty::from_si(CUTOFF), Qty::from_si(SPRING));
    let modes = Modes::of(&network);
    let shape = modes.shape(0);
    let amplitude = (2.0 * 1.380_649e-23 * KELVIN / modes.eigenvalue(0).to_si()).sqrt();
    // Which way is towards the closed structure. The eigenvector's sign is arbitrary -- a mode is
    // a line, not an arrow -- so it is chosen by asking, once, and not by hoping.
    let toward = {
        let at = |s: f64| {
            let moved: Vec<[f64; 3]> = (0..start.len())
                .map(|i| {
                    let p = start.residues()[i].at;
                    [
                        p[0] + s * amplitude * shape[3 * i],
                        p[1] + s * amplitude * shape[3 * i + 1],
                        p[2] + s * amplitude * shape[3 * i + 2],
                    ]
                })
                .collect();
            Structure::from_positions(moved)
                .rmsd_to(closed)
                .expect("the same length")
        };
        if at(REACH) < at(-REACH) {
            1.0
        } else {
            -1.0
        }
    };

    // Built once. The molecule is fixed in this frame -- it is the thing the protein closes on
    // -- so rebuilding 57 spheres on each of 48 frames would be the same arithmetic 2736 times.
    let mut molecule = mesh::Mesh::default();
    for atom in ligand {
        // **Each at its own element'''s radius.** Five phosphorus atoms at 1.80 A against
        // carbon'''s 1.70 and oxygen'''s 1.52 is what a pentaphosphate looks like; one radius for
        // all of them was a picture of a molecule this is not.
        molecule.append(&mesh::sphere(
            atom.at,
            bondi(&atom.element).0,
            ATOM_LEVEL,
            1.0,
        ));
    }
    (0..FRAMES)
        .map(|f| {
            // Out and back, as a cosine, so the ends are still and the middle is quick -- which is
            // what a mode does and also what makes a loop read as a loop rather than a jump.
            // The schedule is `excursion` because `main` walks the same one to check the shapes
            // this draws, and two copies of it would drift.
            let steps = toward * excursion(f);
            let here: Vec<[f64; 3]> = (0..start.len())
                .map(|i| {
                    let p = start.residues()[i].at;
                    [
                        p[0] + steps * amplitude * shape[3 * i],
                        p[1] + steps * amplitude * shape[3 * i + 1],
                        p[2] + steps * amplitude * shape[3 * i + 2],
                    ]
                })
                .collect();
            let gone: Vec<f64> = (0..start.len())
                .map(|i| {
                    let p = start.residues()[i].at;
                    let d = [here[i][0] - p[0], here[i][1] - p[1], here[i][2] - p[2]];
                    (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() / A
                })
                .collect();
            let rmsd = Structure::from_positions(here.clone())
                .rmsd_to(closed)
                .expect("the same length");
            // **Asked of the conformation drawn, not of the one at rest.** A tube whose widest
            // ring outreaches the step it sits on has folded through itself, and it renders as a
            // bead rather than a chain while every other check still passes. `main` asks this of
            // the resting structure; the lid swings 18 A from there, so each frame asks again.
            let (curve, shaded) = mesh::smoothed(&here, &gone, TUBE_SMOOTH);
            let crowding = mesh::crowding(&curve, TUBE_RADIUS);
            assert!(
                crowding < 1.0,
                "frame {f}: the tube's widest ring reaches {crowding:.2}x the step it sits on, so \
                 it has folded through itself -- smooth the trace further or thin the tube"
            );
            let skin = mesh::tube(&curve, &shaded, TUBE_RADIUS, TUBE_SIDES);
            Frame {
                time_s: f as f64 / FRAMES as f64,
                panels: vec![
                    // **The backbone as a solid.** A line has no shape, and the thing this run
                    // is about is a lid closing over a molecule -- which a wire cannot occlude,
                    // cannot shade, and cannot show the near side of. The tube carries one value
                    // per residue to every vertex of that residue's ring, so the colour is per
                    // residue and the renderer's interpolation between rings is what makes it
                    // smooth: no shading is invented here that the numbers do not have.
                    Panel {
                        name: "backbone".into(),
                        unit: "A moved",
                        place: Placed::HERE,
                        data: PanelData::surface(skin.points, skin.faces, skin.values),
                    },
                    Panel {
                        name: "AP5A".into(),
                        unit: "the molecule",
                        place: Placed::HERE,
                        data: PanelData::surface(
                            molecule.points.clone(),
                            molecule.faces.clone(),
                            molecule.values.clone(),
                        ),
                    },
                ],
                readings: vec![
                    Reading::new("closing", "from the closed form", rmsd / A, "A"),
                    Reading::new("closing", "excursion", steps.abs() * amplitude / A, "A"),
                    Reading::new(
                        "closing",
                        "the furthest residue has moved",
                        gone.iter().fold(0.0f64, |m, v| m.max(*v)),
                        "A",
                    ),
                ],
            }
        })
        .collect()
}

/// The middle radius of curvature along a path, as a sense of how bent it is overall.
///
/// The tightest curve says whether the tube folds; the median says whether that tightest one is
/// the chain or one bad turn in it.
fn median_curve(path: &[[f64; 3]]) -> f64 {
    let mut radii: Vec<f64> = (1..path.len() - 1)
        .map(|i| {
            let u = [
                path[i][0] - path[i - 1][0],
                path[i][1] - path[i - 1][1],
                path[i][2] - path[i - 1][2],
            ];
            let v = [
                path[i + 1][0] - path[i][0],
                path[i + 1][1] - path[i][1],
                path[i + 1][2] - path[i][2],
            ];
            let c = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            let n = |a: [f64; 3]| (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
            let w = [u[0] + v[0], u[1] + v[1], u[2] + v[2]];
            if n(c) == 0.0 {
                f64::MAX
            } else {
                n(u) * n(v) * n(w) / (2.0 * n(c))
            }
        })
        .collect();
    radii.sort_by(|a, b| a.partial_cmp(b).expect("no NaN in a measured structure"));
    radii[radii.len() / 2]
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

/// The radius an element is drawn at, and whether [`BONDI`] knew it.
///
/// The second half is the point. An element with no entry is drawn at carbon's radius, which is
/// a reasonable guess and a silent one -- so the caller prints what it fell back on, and a
/// molecule containing a metal says so rather than showing it the size of a carbon.
fn bondi(element: &str) -> (f64, bool) {
    match BONDI.iter().find(|(name, _)| *name == element) {
        Some((_, r)) => (*r, true),
        None => (UNLISTED_RADIUS, false),
    }
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
