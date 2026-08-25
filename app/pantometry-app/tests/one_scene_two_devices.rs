//! **The same scene, on the CPU and on the device, compared.**
//!
//! This is the test the merge existed to make possible, and it found three defects on its first
//! run. Until the CLI, the viewer, the editor and the accelerator were one application there was no
//! binary that could build a scene *and* attach a device, so a `"device": "gpu"` scene had never
//! been run end to end by anything — only by unit tests that reached for `&mut` accessors and
//! therefore refreshed state the trait methods could not.
//!
//! What it found, all three in the same run:
//!
//! - **The mirror was never synced.** `ledger`, `readings` and `as_field` are all `&self`, and
//!   `Simulation::advance` never hands back a `&mut`, so the device's cells were read back once at
//!   construction and never again. Measured: `peak` read `120.0000011 °C` at every one of five
//!   frames while the CPU's spot diffused to `20.34`.
//! - **`as_field` was missing**, so `capture` made no panel: no report view, no CSV field, no glTF,
//!   no USD geometry, no viewer. The CLI said "no field and no bodies — not drawn" about a block.
//! - **The ledger's scale was its own net**, which is `Σ cᵢ(Tᵢ − T̄)` — identically zero for equal
//!   capacities. A relative tolerance against zero is not a tolerance, and once the mirror moved it
//!   fired on a correct run at `1.487e0`.
//!
//! Every test here skips with a printed reason when there is no adapter. A CI runner has none, and
//! a software rasteriser would be checking a different implementation than anyone uses.

use pantometry_gpu::OnTheGpu;
use pantometry_world::{OnDisk, Scene, World};
use viewer_core::Run;

/// A hot cell in the middle of a block, which is the cheapest scene whose answer moves.
fn scene(device: &str) -> Scene {
    let text = format!(
        r#"{{
      "title": "a hot spot, two ways",
      "duration_s": 0.05,
      "frames": 5,
      "conservation_tolerance": 1e-4,
      "domains": [
        {{ "kind": "block", "name": "part", "cells": [16, 16, 16], "cell_mm": 1.0,
           "material": "aluminium", "initial_c": 20.0{device},
           "hot_spot": {{ "at": [8, 8, 8], "above_k": 100.0 }} }}
      ]
    }}"#
    );
    serde_json::from_str(&text).expect("the scene parses")
}

/// Run it both ways. `None` when this machine has no adapter, with a reason on stdout.
fn both() -> Option<(Vec<pantometry_world::Frame>, Vec<pantometry_world::Frame>)> {
    let asked = scene(", \"device\": \"gpu\"");
    let mut on_device = match World::build_with_accelerator(asked, &OnDisk, &OnTheGpu) {
        Ok(w) => w,
        Err(why) => {
            println!("skipped: {why}");
            return None;
        }
    };
    let device = on_device
        .run()
        .expect("the run holds its books on the device");
    let mut reference = World::build(scene("")).expect("the cpu scene builds");
    let cpu = reference.run().expect("the run holds its books on the cpu");
    Some((device, cpu))
}

/// **The device's answer moves, and moves to the CPU's answer.**
#[test]
fn the_device_reports_what_it_computed_and_not_what_it_started_with() {
    let Some((device, cpu)) = both() else { return };
    assert_eq!(device.len(), cpu.len(), "the same scene, the same frames");
    assert!(device.len() >= 3, "a run worth comparing");

    let peak = |frames: &[pantometry_world::Frame], i: usize| {
        frames[i]
            .readings
            .iter()
            .find(|r| r.label == "peak")
            .map(|r| r.value)
            .expect("a block reports its peak")
    };

    // The claim the frozen mirror failed: the first frame and the last are different runs of the
    // same physics, and a domain that reports its initial condition forever passes every
    // agreement check against itself.
    let (first, last) = (peak(&device, 0), peak(&device, device.len() - 1));
    println!("  device peak {first:.4} → {last:.4} °C");
    assert!(
        (first - last).abs() > 1.0,
        "the device reported {first} at both ends of the run — the mirror is not being read back"
    );

    for i in 0..device.len() {
        let (g, c) = (peak(&device, i), peak(&cpu, i));
        // Single precision against double, over the run: about `1e-7` a step walking as `√k`. A
        // thousandth is that with room, and far tighter than a *different operator* would give.
        assert!(
            (g - c).abs() / c.abs().max(1e-9) < 1e-3,
            "frame {i}: device {g} against cpu {c}"
        );
    }
}

/// **The device draws.** A run with no panel has no report, no export and no viewer.
#[test]
fn a_block_on_the_device_still_has_a_picture() {
    let Some((device, cpu)) = both() else { return };
    let last = device.last().expect("a last frame");

    assert_eq!(
        last.panels.len(),
        cpu.last().expect("a last frame").panels.len(),
        "the device offers the same panels as the cpu"
    );
    let panel = last.panels.first().expect("a block is drawn");
    println!("  the device's panel: {} {}", panel.name, panel.unit);

    // Through the run file's own reader, which is what every consumer of this actually uses — the
    // viewer, the editor and the HTML report all go through `viewer_core::Run`.
    let json = pantometry::view::to_json("a hot spot, two ways", &device);
    let parsed: Run = serde_json::from_str(&json).expect("the run file parses");
    let field = parsed
        .frames
        .last()
        .and_then(|f| f.panels.first())
        .expect("a panel survives the round trip");
    let viewer_core::Panel::Field { nx, ny, nz, .. } = field else {
        panic!("a block is a field, not {}", field.name());
    };
    assert_eq!(
        (*nx, *ny, *nz),
        (16, 16, 16),
        "the whole grid, not a corner of it"
    );

    // And it is the run's own values, not the initial condition: the middle is hotter than a face.
    let values = field.values();
    let middle = values[8 + 16 * (8 + 16 * 8)];
    let corner = values[0];
    println!("  middle {middle:.4} against corner {corner:.4}");
    assert!(
        middle > corner + 0.01,
        "the spot is not in the picture: middle {middle}, corner {corner}"
    );
}
