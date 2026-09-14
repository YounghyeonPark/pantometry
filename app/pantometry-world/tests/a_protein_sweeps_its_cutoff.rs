//! The sweep a protein has, in place of the one it cannot have.
//!
//! Every other spatial scene here states a grid, and [`Scene::refined`] halves it to ask whether
//! the answer belongs to the problem or to the discretisation. A protein has no grid. The
//! structure is where the atoms were measured to be, and there is nothing between two of them to
//! subdivide — so `refined()` refuses with the same message a scene of pure systems gets, and a
//! scene whose only check was that sweep would ship measuring nothing.
//!
//! What can quietly be doing the work is the **cutoff**: the one number in the scene that decides
//! which pairs are joined, and therefore what the network is. If the softest mode turned with it,
//! the motion the scene draws would be a statement about a parameter rather than about the fold.
//! So that is what is swept here, and the comparison is the same one the model is judged by
//! elsewhere — `|cos|` between mode shapes, absolute because a mode's sign is arbitrary.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped scene off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use pantometry::protein::{Modes, Network, Structure};
use pantometry_world::{Beside, DomainSpec, Scene, World};

/// The shipped scene, and the directory its structure file sits beside.
fn scene() -> (Scene, Beside) {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("scenes");
    let path = dir.join("31-a-protein-shaking-at-body-temperature.json");
    let text = std::fs::read_to_string(&path).expect("the shipped scene reads");
    (
        serde_json::from_str(&text).expect("the shipped scene parses"),
        Beside(dir),
    )
}

/// The scene's own protein, at whatever cutoff is asked for.
fn modes_at(cutoff_a: f64) -> (Network, Modes) {
    let (scene, beside) = scene();
    let DomainSpec::Protein { pdb, chain, .. } = &scene.domains[0] else {
        panic!("the scene's first domain is not the protein it was written to be");
    };
    let bytes = pantometry_world::Parts::bytes(&beside, pdb).expect("the structure reads");
    let whole =
        Structure::from_pdb(&String::from_utf8(bytes).expect("text")).expect("alpha carbons");
    let structure = match chain {
        Some(id) => whole.chain(*id).expect("the chain the scene names"),
        None => whole,
    };
    let network = Network::new(
        &structure,
        pantometry::units::Qty::from_si(cutoff_a * 1e-10),
        pantometry::units::Qty::from_si(0.695),
    );
    let modes = Modes::of(&network);
    (network, modes)
}

/// **The softest mode is still the softest mode from 10 Å to 20 Å, and that is asserted without
/// a threshold.**
///
/// This is the scene's sweep. A resolution sweep asks whether halving the cells changes the
/// answer; this asks whether the one parameter that decides what the network *is* changes the
/// motion. Overlap is `|cos|` between mode shapes in `3N` space, absolute because a mode's sign
/// is arbitrary, and in `3 × 46 = 138` dimensions a direction picked without looking scores
/// `1/√138 = 0.085`.
///
/// # Why the claim is an ordering and not a number
///
/// The first version of this test asked for an overlap above `0.9`. Measured against the 15 Å
/// reference: **0.815** at 10 Å, **0.927** at 12, **0.849** at 18, **0.809** at 20 — each about
/// ten times chance, and three of the four below the bar. Moving the bar to fit them would be
/// choosing a number because it made the test pass, which is the thing `CONTRIBUTING.md` names
/// first.
///
/// So the claim is the one the numbers actually support and the scene actually needs: the softest
/// mode at every cutoff is **best matched by the softest mode** at the reference, out of ten
/// candidates. That is a statement about identity rather than about distance, it has no parameter
/// in it, and it is what would break if the cutoff were choosing the motion — below percolation
/// the softest mode becomes one piece of a broken network drifting away from another, which is
/// not mode 0 of anything and would lose to a different mode instantly.
#[test]
fn the_softest_mode_does_not_turn_with_the_cutoff() {
    let (reference_net, reference) = modes_at(15.0);
    assert_eq!(reference.rigid_body_modes(), 6, "crambin is not collinear");
    println!(
        "  15 A: {} springs, {} modes (the reference)",
        reference_net.springs(),
        reference.len()
    );

    let chance = 1.0 / (3.0 * 46.0f64).sqrt();
    let mut tightest = f64::INFINITY;
    for cutoff in [10.0, 12.0, 18.0, 20.0] {
        let (net, other) = modes_at(cutoff);
        assert_eq!(
            other.rigid_body_modes(),
            6,
            "at {cutoff} A the network is not one piece"
        );
        // This cutoff's softest mode against the reference's ten softest.
        let against: Vec<f64> = (0..10).map(|k| reference.overlap(k, &other, 0)).collect();
        let (best, &score) = against
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .expect("ten candidates");
        let runner_up = against
            .iter()
            .enumerate()
            .filter(|(k, _)| *k != best)
            .map(|(_, v)| *v)
            .fold(0.0f64, f64::max);
        println!(
            "  {cutoff:4.0} A: {:4} springs, softest mode matches reference mode {best} at \
             {score:.3}, next best {runner_up:.3}",
            net.springs()
        );
        assert_eq!(
            best, 0,
            "at {cutoff} A the softest mode is reference mode {best} rather than mode 0, so the \
             cutoff is choosing the motion this scene draws"
        );
        tightest = tightest.min(score - runner_up);
    }
    println!("  closest call: the right mode won by {tightest:.3}, and chance is {chance:.3}");
    // The identification has to be *unambiguous*, not merely correct: a margin of nothing would
    // mean the loop above was reporting a coin toss it happened to win four times.
    assert!(
        tightest > chance,
        "the softest mode won by only {tightest:.3}, which is under the {chance:.3} a direction \
         that knew nothing would score, so calling it the same mode is not supported"
    );
}

/// **The scene has no resolution to halve, and says so rather than skipping quietly.**
///
/// The failure this guards is the one commit `7374d3b` was written for: a scene that silently
/// has no sweep is a scene nobody measured, and it looks exactly like a scene that passed one.
#[test]
fn the_scene_refuses_to_refine_and_names_the_reason() {
    let (scene, _) = scene();
    let why = scene
        .refined()
        .map(|_| ())
        .expect_err("a protein is not a discretisation of anything");
    assert!(
        why.contains("nothing in this scene is a discretisation"),
        "it refused for some other reason: {why}"
    );
    println!("  refined(): {why}");
}

/// **The scene builds, runs, and the protein is the size the file says.**
#[test]
fn the_shipped_scene_builds_and_moves() {
    let (scene, beside) = scene();
    let mut world = World::build_with(scene, &beside).expect("the shipped scene builds");
    let before = world.capture();
    world.run().expect("the shipped scene runs");
    let after = world.capture();

    let bodies = |f: &pantometry::scene::Frame| {
        f.panels
            .iter()
            .find(|p| p.name == "crambin")
            .map(|p| p.values().len())
            .expect("a panel for the protein")
    };
    assert_eq!(bodies(&before), 46, "crambin has 46 residues");
    assert_eq!(bodies(&after), 46);

    // It has to have moved. A structure sitting still would draw the same picture every frame and
    // every check above would still pass.
    let positions = |f: &pantometry::scene::Frame| {
        let panel = f
            .panels
            .iter()
            .find(|p| p.name == "crambin")
            .expect("a panel for the protein");
        match &panel.data {
            pantometry::scene::PanelData::Points { positions, .. } => positions.clone(),
            other => panic!("the protein draws as {other:?} rather than as points"),
        }
    };
    let (a, b) = (positions(&before), positions(&after));
    let moved = a
        .iter()
        .zip(&b)
        .map(|(p, q)| (0..3).map(|k| (p[k] - q[k]).abs()).fold(0.0f64, f64::max))
        .fold(0.0f64, f64::max);
    assert!(
        moved > 1e-12,
        "nothing moved: the largest displacement over the whole run was {moved:e} m"
    );
    println!("  46 residues, largest displacement over the run {moved:.3e} m");
}

/// **The structure beside the scenes is the one the library crate is checked against.**
///
/// `crates/pantometry-protein/structures/1CRN.pdb` is a test fixture and this is a scene's input,
/// and they are the same deposited entry. Two copies of a file are two files that can drift, and
/// the one that would drift is the one nothing compares — so this compares them, byte for byte.
/// Neither is generated: both are the Protein Data Bank's own bytes, which is what makes the
/// comparison meaningful rather than circular.
#[test]
fn the_two_copies_of_crambin_are_one_file() {
    let here = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let scene_copy = here.join("scenes/structures/1CRN.pdb");
    let crate_copy = here.join("../../crates/pantometry-protein/structures/1CRN.pdb");
    let (a, b) = (
        std::fs::read(&scene_copy).expect("the scene's copy"),
        std::fs::read(&crate_copy).expect("the crate's copy"),
    );
    assert_eq!(
        a.len(),
        b.len(),
        "{} is {} bytes and {} is {}",
        scene_copy.display(),
        a.len(),
        crate_copy.display(),
        b.len()
    );
    assert!(a == b, "the two copies of 1CRN.pdb have diverged");
    // And neither is empty, which two missing files would also satisfy.
    assert!(a.len() > 40_000, "1CRN.pdb is only {} bytes", a.len());
    println!("  both copies of 1CRN.pdb are {} bytes", a.len());
}
