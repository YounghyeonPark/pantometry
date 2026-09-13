//! The model against a measurement nobody in this repository made.
//!
//! Every other test in this crate checks the code against a theorem: a rigid motion is in the
//! null space, a chain's spectrum is `4 sin²(jπ/2N)`, a ligand cannot make a residue more mobile.
//! Those say the arithmetic is what it claims to be. **None of them says it describes a
//! protein.** A network of springs satisfies all of them and so does a network of springs that
//! has nothing to do with any molecule.
//!
//! What says otherwise is in the file: a crystallographer refined a temperature factor for every
//! atom, `B = 8π²⟨Δr²⟩/3`, and that is an independent measurement of the thing this model
//! predicts. Convention 1 of this workspace asks for a closed form and never for a second
//! implementation; this is a third kind, a number from the world.
//!
//! # It is a weak check, and this file is where that is written down
//!
//! [`the_trivial_predictors_do_as_well`] measures two predictors that are one line of geometry
//! each — the distance from the centroid, and how many neighbours a residue has — against the
//! same B-factors. **They match or beat the model on all four structures.** So agreement here is
//! necessary and is not evidence that the collective-motion physics is right: a B-factor is one
//! scalar per residue, and a scalar cannot distinguish a model that knows which way things move
//! from one that knows only which residues are on the outside.
//!
//! The check that a trivial predictor cannot pass is in
//! `a_protein_caught_in_two_positions.rs`, which compares against a **direction** rather than a
//! magnitude, and where the model beats a random direction by twenty-five to one.
//!
//! # What it can and cannot be asked for
//!
//! Not the **size** of the fluctuations. The spring constant divides all of them equally and
//! nothing here predicts it, so the amplitude is a scale and not a result. What is a result is
//! the **shape**: whether the same residues are the mobile ones. A correlation coefficient is
//! blind to scale, which is why it and not a residual is the comparison.
//!
//! # The cutoff is the published one
//!
//! 15 Å, the value the anisotropic network model was published with. Choosing it to make a
//! correlation larger would turn a prediction into a one-parameter fit, and it would work:
//! crambin goes from `0.395` at 15 Å to `0.685` at 7 Å. The number stays where the literature
//! put it.
//!
//! # Not under `wasm32`
//!
//! These structures are `include_str!`, so they compile in anywhere — but see
//! [`a_modern_and_an_old_refinement_of_the_same_protein_disagree`] for the one test here that is
//! release-only, and why.

use pantometry_protein::{correlation, Modes, Network, Structure};
use pantometry_units::{Qty, Stiffness, Temperature};

/// Ångström, in metres.
const A: f64 = 1e-10;

/// The published cutoff for the anisotropic network model.
const CUTOFF: f64 = 15.0 * A;

/// One kcal/mol/Å². It cancels out of every correlation below; it is here because a `Network`
/// needs one, not because anything depends on it.
fn gamma() -> Stiffness {
    Qty::from_si(0.695)
}

fn room() -> Temperature {
    Qty::from_si(300.0)
}

const CRAMBIN: &str = include_str!("../structures/1CRN.pdb");
const UBIQUITIN: &str = include_str!("../structures/1UBQ.pdb");

/// Predicted against deposited, for one entry.
fn agreement(text: &str, expect_id: &str, expect_residues: usize) -> f64 {
    // A fixture that has been swapped, truncated or regenerated is a fixture that is no longer
    // the entry the numbers below were measured against. The entry names itself in columns
    // 63-66 of its own HEADER.
    let header = text.lines().next().expect("a HEADER line");
    assert!(
        header.starts_with("HEADER") && header.len() >= 66 && &header[62..66] == expect_id,
        "this file's HEADER says {:?}, not {expect_id}",
        header.get(62..66)
    );

    let s = Structure::from_pdb(text).expect("alpha carbons");
    assert_eq!(s.len(), expect_residues, "{expect_id} changed size");
    let measured = s
        .experimental_fluctuations()
        .expect("a deposited structure has B-factors");

    let network = Network::new(&s, Qty::from_si(CUTOFF), gamma());
    assert_eq!(
        network.components(),
        1,
        "{expect_id} falls into {} pieces at 15 A, so its fluctuations are a drift",
        network.components()
    );
    let predicted = Modes::of(&network).fluctuations(room());
    correlation(&predicted, &measured).expect("neither series is constant")
}

/// **Two proteins the model should get right, and does.**
///
/// Ubiquitin at `0.571` and crambin at `0.395`. Both are pinned to narrow windows around what
/// was measured rather than to a wide band, because a bound loose enough to hold both would be
/// loose enough to hold noise.
///
/// Crambin being the lower of the two is what a reader would expect and is not evidence of
/// anything: it is forty-six residues stapled by three disulfides, about as rigid as a protein
/// gets, and what little its B-factors vary by is largely crystal contact. That is an explanation
/// offered after the fact, which is the weakest kind, and it is written here as one rather than as
/// a reason the number is right.
#[test]
fn the_predicted_fluctuations_agree_with_the_deposited_ones() {
    let crambin = agreement(CRAMBIN, "1CRN", 46);
    let ubiquitin = agreement(UBIQUITIN, "1UBQ", 76);
    println!("  1CRN  46 residues: {crambin:.3}");
    println!("  1UBQ  76 residues: {ubiquitin:.3}");

    assert!(
        (0.30..0.50).contains(&crambin),
        "crambin came out at {crambin:.3}, and it was 0.395"
    );
    assert!(
        (0.50..0.65).contains(&ubiquitin),
        "ubiquitin came out at {ubiquitin:.3}, and it was 0.571"
    );
    // The big one lives in the release-only test below, because it costs 18 s in a debug build.
    println!("  lysozyme is in a_modern_and_an_old_refinement_of_the_same_protein_disagree");
}

/// **The same protein, refined twice, and the model agrees with one of them.**
///
/// 193L and 4LYZ are both hen egg-white lysozyme: the same sequence, the same fold, the same 129
/// residues. The model predicts `0.659` against the 1996 refinement and `−0.404` against the 1975
/// one.
///
/// That is the most useful number in this crate, and it is the one a check written only against
/// structures that agree would have hidden. Two things follow from it. The comparison is
/// sensitive to the **measurement** rather than being a property of the model, since the model
/// gave both answers about one protein. And a correlation near zero or below is a statement about
/// a file: 4LYZ's B-factors span 0.46 to 12.35 Å² where 193L's span 9.1 to 32.9, and its
/// `REMARK 3` names its refinement program as `NULL`.
///
/// # Release only, and what that costs
///
/// A 129-residue protein is a `387 × 387` eigenproblem, which takes **18 s** in a debug build and
/// **0.5 s** in a release one — measured, both, on the machine this was written on. The library
/// gate and CI both run `cargo test --workspace --release`, so this runs on every commit; what it
/// does not do is add 36 s to the debug pass as well.
#[test]
#[cfg(not(debug_assertions))]
fn a_modern_and_an_old_refinement_of_the_same_protein_disagree() {
    // Inside the test rather than beside its neighbours: `cfg(not(debug_assertions))` deletes
    // this function in a debug build, and a `const` it was the only user of becomes dead code
    // that `-D warnings` refuses. Clippy runs in debug, so that is not a hypothetical.
    let modern = agreement(include_str!("../structures/193L.pdb"), "193L", 129);
    let old = agreement(include_str!("../structures/4LYZ.pdb"), "4LYZ", 129);
    println!("  193L 129 residues, 1996 refinement: {modern:.3}");
    println!("  4LYZ 129 residues, 1975 refinement: {old:.3}");

    assert!(
        (0.60..0.72).contains(&modern),
        "the modern refinement came out at {modern:.3}, and it was 0.659"
    );
    assert!(
        old < 0.0,
        "4LYZ came out at {old:.3}; it was -0.404, and a positive number here would mean the \
         two refinements no longer disagree -- which is a thing to read about before believing"
    );
    assert!(
        modern - old > 1.0,
        "the two refinements are {:.3} apart and were 1.063",
        modern - old
    );
}

/// **The trivial predictors do as well, and this is the honest place to say so.**
///
/// Three substitutes for the predicted fluctuations:
///
/// - **the residue index**, which is what a model that had lost the structure and kept the order
///   would give — the one that should lose, and does,
/// - **the distance from the centroid**, one line of geometry,
/// - **minus the number of neighbours within the cutoff**, one more line.
///
/// Measured, at the published 15 Å cutoff:
///
/// | | model | centroid | neighbours | index |
/// | --- | --- | --- | --- | --- |
/// | 1CRN | +0.395 | +0.402 | +0.423 | +0.375 |
/// | 1UBQ | +0.571 | **+0.804** | **+0.768** | +0.346 |
/// | 193L | +0.659 | **+0.699** | **+0.697** | +0.262 |
/// | 4LYZ | −0.404 | −0.301 | −0.267 | −0.234 |
///
/// **The model does not beat either geometric predictor on any structure here.** That is a real
/// result and it is not news in the field — a residue's B-factor is mostly a statement about how
/// buried it is — but it is the kind of thing a crate reports about itself only if somebody wrote
/// the check. The assertion below is therefore not "the model wins": it is that the ordering is
/// what was measured, so that a change which alters it has to be read rather than merely pass.
///
/// What follows from it is in this file's own documentation: the B-factor comparison is
/// necessary and is not sufficient, and the sufficient one needs a direction.
#[test]
fn the_trivial_predictors_do_as_well() {
    let s = Structure::from_pdb(UBIQUITIN).expect("ubiquitin");
    let measured = s.experimental_fluctuations().expect("B-factors");
    let network = Network::new(&s, Qty::from_si(CUTOFF), gamma());
    let predicted = Modes::of(&network).fluctuations(room());
    let model = correlation(&predicted, &measured).expect("a correlation");

    let index: Vec<f64> = (0..s.len()).map(|i| i as f64).collect();
    let mut middle = [0.0; 3];
    for r in s.residues() {
        for (k, slot) in middle.iter_mut().enumerate() {
            *slot += r.at[k] / s.len() as f64;
        }
    }
    let centroid: Vec<f64> = s
        .residues()
        .iter()
        .map(|r| {
            (0..3)
                .map(|k| (r.at[k] - middle[k]).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect();
    let neighbours: Vec<f64> = s
        .residues()
        .iter()
        .map(|a| {
            let count = s
                .residues()
                .iter()
                .filter(|b| {
                    let d: f64 = (0..3).map(|k| (a.at[k] - b.at[k]).powi(2)).sum();
                    d > 0.0 && d <= CUTOFF * CUTOFF
                })
                .count() as f64;
            -count
        })
        .collect();

    let score = |v: &[f64]| correlation(v, &measured).expect("a correlation");
    let (c, n, i) = (score(&centroid), score(&neighbours), score(&index));
    println!("  1UBQ: model {model:+.3}   centroid {c:+.3}   neighbours {n:+.3}   index {i:+.3}");

    // The one that should lose: a predictor that kept the sequence and forgot the shape.
    assert!(
        i < model,
        "the residue index scored {i:.3} against the model's {model:.3}, so the fluctuations are \
         not using the structure at all"
    );
    // And the two that do not. Asserted in the direction they were measured in, because that is
    // the claim this crate is entitled to make: a change that made the model win here would be a
    // change worth reading, not one to wave through.
    assert!(
        c > model && n > model,
        "the geometric predictors no longer beat the model (centroid {c:.3}, neighbours {n:.3}, \
         model {model:.3}); they did, at 0.804 and 0.768 against 0.571, and the documentation of \
         this crate says so in three places"
    );
}
