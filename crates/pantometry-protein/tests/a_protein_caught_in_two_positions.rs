//! The model against a motion a protein was photographed performing.
//!
//! This is the check `a_protein_a_crystallographer_measured.rs` cannot do, and the reason that
//! one is not enough. A B-factor is **one number per residue**, so a comparison against B-factors
//! asks only whether the model gets the mobile residues right — and measured, on all four
//! deposited structures in this crate, the distance from the centroid gets them just as right. A
//! model earns nothing by matching a predictor that is one line of geometry.
//!
//! What an elastic network model actually claims is a **direction in `3n` space**: not how much a
//! residue moves but which way, and together with what. Nothing scalar can test that, and a
//! random direction cannot fake it — in `3n = 642` dimensions a random unit vector overlaps any
//! fixed one by `1/√642 = 0.039`.
//!
//! # Adenylate kinase, open and closed
//!
//! The enzyme closes two lids over its substrates and opens again, and both ends of that motion
//! are in the Protein Data Bank: `4AKE` open, `1AKE` closed on its inhibitor. Superposed, they
//! are **7.1 Å** apart — a large conformational change, not a vibration.
//!
//! The model is given only the **open** structure. It has never seen the closed one, and the
//! closed one contributes nothing but the direction it is asked about.
//!
//! # Release only
//!
//! 214 residues is a `642 × 642` eigenproblem: about 4 s in a release build and minutes in a
//! debug one. The library gate and CI both run `cargo test --workspace --release`.

#![cfg(not(debug_assertions))]

use pantometry_protein::{Modes, Network, Structure};
use pantometry_units::Qty;

/// Ångström, in metres.
const A: f64 = 1e-10;

/// Chain A of an entry, checked to be the entry it says it is.
fn chain_a(text: &str, expect_id: &str) -> Structure {
    let header = text.lines().next().expect("a HEADER line");
    assert!(
        header.starts_with("HEADER") && header.len() >= 66 && &header[62..66] == expect_id,
        "this file's HEADER says {:?}, not {expect_id}",
        header.get(62..66)
    );
    let whole = Structure::from_pdb(text).expect("alpha carbons");
    // Both entries are dimers in the asymmetric unit. Building one network across both joins two
    // molecules that are not joined, and every mode after that is about the crystal.
    assert_eq!(whole.len(), 428, "{expect_id} is no longer a dimer");
    let one = whole.chain('A').expect("a chain A");
    assert_eq!(one.len(), 214);
    one
}

/// **The softest mode of the open enzyme is the motion it performs when it closes, and a
/// direction that knows nothing about the protein gets `1/√3n`.**
///
/// Overlap `0.799` for the single lowest non-rigid mode, against `0.039` for a random direction.
/// Ten modes reach `0.966`. Nothing in this repository was tuned to produce that: the cutoff is
/// the same 15 Å every other test here uses, the spring constant cancels out, and the model is
/// never shown the closed structure.
///
/// The comparison this is read against is the one measured below — a direction that knows nothing
/// about the protein — and not a range recalled from a paper. Every number in this crate is one
/// its own test printed.
///
/// This is the answer to "what does the structure do when something binds", and it is the one
/// claim in this crate that a scalar comparison could not have supported.
#[test]
fn the_lowest_mode_of_the_open_enzyme_points_at_the_closed_one() {
    let open = chain_a(include_str!("../structures/4AKE.pdb"), "4AKE");
    let closed = chain_a(include_str!("../structures/1AKE.pdb"), "1AKE");

    let before = open.rmsd_to(&closed).expect("the same length");
    let fitted = closed.superposed_onto(&open).expect("the same length");
    let after = open.rmsd_to(&fitted).expect("the same length");
    assert!(
        after < before,
        "superposing made it worse: {:.3} to {:.3} A",
        before / A,
        after / A
    );
    assert!(
        (6.0 * A..8.0 * A).contains(&after),
        "the two conformations are {:.3} A apart and were 7.131",
        after / A
    );

    let went = open
        .direction_to(&fitted)
        .expect("two different structures");
    let network = Network::new(&open, Qty::from_si(15.0 * A), Qty::from_si(0.695));
    assert_eq!(network.components(), 1);
    let modes = Modes::of(&network);
    assert_eq!(modes.rigid_body_modes(), 6);

    let lowest = modes.overlap_with(0, &went);
    let mut cumulative = 0.0;
    for k in 0..10 {
        let o = modes.overlap_with(k, &went);
        cumulative += o * o;
        println!("    mode {k}: {o:.3}   cumulative {:.3}", cumulative.sqrt());
    }
    let ten = cumulative.sqrt();
    // `1/√3n`: what a direction chosen without looking would give.
    let chance = 1.0 / (3.0 * open.len() as f64).sqrt();
    println!(
        "  4AKE -> 1AKE: {:.3} A apart, lowest mode {lowest:.3}, ten modes {ten:.3}, chance \
         {chance:.3}",
        after / A
    );

    assert!(
        lowest > 0.70,
        "the lowest mode overlaps the observed motion by {lowest:.3}, and it was 0.799"
    );
    assert!(
        ten > 0.90,
        "ten modes reach {ten:.3} of the observed motion, and they reached 0.966"
    );
    assert!(
        lowest > 15.0 * chance,
        "{lowest:.3} is only {:.1} times what a random direction gives",
        lowest / chance
    );

    // **And the null hypothesis, measured rather than quoted.** Without it `0.799` is a number
    // with no units. It shares this test's eigendecomposition on purpose: a second one of a
    // `642 × 642` matrix costs as much as the whole rest of the crate's suite.
    let n = 3 * open.len();
    let mut worst: f64 = 0.0;
    let mut mean = 0.0;
    let trials = 200;
    for trial in 0..trials {
        let mut rng = pantometry_core::Rng::for_index(20260913, trial);
        let mut v: Vec<f64> = (0..n).map(|_| rng.gaussian()).collect();
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        for x in v.iter_mut() {
            *x /= norm;
        }
        let o = modes.overlap_with(0, &v);
        worst = worst.max(o);
        mean += o / trials as f64;
    }
    println!(
        "  {trials} random directions against mode 0: mean {mean:.4}, worst {worst:.4}, \
         1/sqrt(3n) = {chance:.4}"
    );
    // The mean of `|cos|` between a random unit vector and a fixed one in `n` dimensions is
    // `√(2/πn)`, which is `0.798/√n` — so this is a closed form and not a rule of thumb.
    let expected = (2.0 / (std::f64::consts::PI * n as f64)).sqrt();
    assert!(
        (mean / expected - 1.0).abs() < 0.15,
        "the mean random overlap is {mean:.4} against the closed form {expected:.4}"
    );
    assert!(
        worst < lowest,
        "a random direction reached {worst:.3} against the model's {lowest:.3}"
    );
}
