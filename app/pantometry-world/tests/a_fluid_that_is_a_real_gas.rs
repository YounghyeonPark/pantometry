//! **A fluid with no substance is in reduced units, and used to say it was in metres.**
//!
//! `LennardJones::reduced()` sets `σ = ε = m = 1`. That is the right thing for reproducing a
//! paper's state point and the wrong thing for a file: every length the run reports is then a
//! number of `σ`, written into a wire format whose positions are metres. The shipped fluid scene
//! declared a box `5.0388 x 5.0388 x 5.0388 m` for 108 atoms that span 1.7 nanometres, and
//! `pantometry view` drew it under a scale bar reading **`2 M`**. The bar was right about the
//! numbers it was given.
//!
//! A scene can name its gas now. These check the two halves of that: the named one arrives in SI,
//! and the unnamed one stops claiming to be.
//!
//! # Why this is here and not in a scene
//!
//! A scene of argon needs its `duration_s` in seconds — argon's `τ = σ√(m/ε)` is 2.16 ps, so the
//! shipped `6.0` would be `10¹²` of them. A new scene is also a number that appears in four
//! documents. This goes through `World::build` instead, which is the same path a scene takes and
//! runs on every commit.

use pantometry_world::{Scene, World};

/// A fluid scene, optionally naming a gas.
fn scene_text(substance: &str) -> String {
    format!(
        r#"{{
  "title": "a small box of atoms",
  "duration_s": 1.0e-11,
  "frames": 2,
  "domains": [
    {{ "kind": "atoms", "name": "gas", "cells": 2, "density": 0.8442,
       "temperature": 1.4, "seed": 1{substance} }}
  ]
}}"#
    )
}

/// The refusal, or a panic naming what was accepted. `expect_err` needs `Debug` on the success
/// type and a `World` has no business having one -- the same helper `declared_materials` keeps.
fn refusal(substance: &str) -> String {
    match build(substance) {
        Err(e) => e,
        Ok(_) => panic!("built, and should not have: {substance}"),
    }
}

fn build(substance: &str) -> Result<World, String> {
    let scene: Scene =
        serde_json::from_str(&scene_text(substance)).map_err(|e| format!("parse: {e}"))?;
    World::build(scene)
}

/// **Argon arrives in metres, and the box is the size argon's box is.**
///
/// 32 atoms at `ρ* = 0.8442` occupy `N/ρ*` cubic `σ`, so the side is `(32/0.8442)^(1/3) = 3.3592 σ`
/// — and with argon's `σ = 3.405 Å` that is **1.145 nm**. The assertion is on the number, not on
/// the fact that it built: a substance that parsed and was then dropped would build perfectly
/// well and give a box 3.363 metres across, which is what this used to do.
#[test]
fn a_named_gas_puts_the_run_in_metres() {
    let world = build(r#", "substance": "argon""#).expect("argon builds");
    let fluid = world
        .simulation()
        .domain_as::<pantometry::molecular::Fluid>("gas")
        .expect("the fluid is there");

    let side = fluid.bounds().length;
    let expected = (32.0f64 / 0.8442).cbrt() * 3.405e-10;
    println!("  the box is {side:e} m, and argon's is {expected:e}");
    assert!(
        (side - expected).abs() < 1e-3 * expected,
        "32 argon atoms at rho* = 0.8442 make a box {expected:e} m across, and this one is {side:e}"
    );
    // A nanometre and not a metre, said plainly, because the failure this replaces was exactly a
    // factor of 10^9 that nothing in the file objected to.
    assert!(
        (1e-9..2e-9).contains(&side),
        "a box of 32 argon atoms is about a nanometre across, and this one is {side:e} m"
    );

    use pantometry::prelude::Bodies;
    assert_eq!(
        fluid.value_unit(),
        "m/s",
        "a fluid with a real sigma and epsilon reports real speeds"
    );
}

/// **And without a substance it says `sigma/tau` rather than `m/s`.**
///
/// Reduced is a legitimate thing to be — it is how a paper's state point is reproduced and how
/// this crate's own tests are written. What is not legitimate is calling it metres.
#[test]
fn an_unnamed_fluid_says_it_is_in_reduced_units() {
    let world = build("").expect("a reduced fluid builds");
    let fluid = world
        .simulation()
        .domain_as::<pantometry::molecular::Fluid>("gas")
        .expect("the fluid is there");

    use pantometry::prelude::Bodies;
    assert_eq!(
        fluid.value_unit(),
        "sigma/tau",
        "a fluid at sigma = epsilon = 1 is in reduced units and has to say so"
    );
    let side = fluid.bounds().length;
    println!("  the reduced box is {side} sigma");
    assert!(
        (3.359..3.360).contains(&side),
        "32 atoms at rho* = 0.8442 make a box 3.3592 sigma across, and this one is {side}"
    );
}

/// **A gas this potential does not describe is refused by name**, with the list beside it.
///
/// Not a default and not the nearest match: a scene asking for nitrogen and quietly getting argon
/// would be a run of the wrong gas reported in the right units, which is worse than a refusal.
#[test]
fn a_substance_this_potential_cannot_describe_is_refused() {
    let why = refusal(r#", "substance": "nitrogen""#);
    println!("  {why}");
    assert!(
        why.contains("nitrogen") && why.contains("argon"),
        "the refusal should name what was asked for and what it could have been: {why}"
    );
    // Neon and xenon are measured in `pantometry-molecular`'s own suite and are not offered here
    // — neon because the potential is wrong for it, xenon because its parameters could not be
    // sourced. A scene that names either gets the same refusal as nitrogen.
    for out in ["neon", "xenon"] {
        let why = refusal(&format!(r#", "substance": "{out}""#));
        assert!(why.contains(out), "the refusal should name {out}: {why}");
    }
}
