//! **A scene states the device, and the answer comes from where it said.**
//!
//! End to end: the JSON carries `"device": "gpu"`, `World::build_with_accelerator` honours it
//! through `OnTheGpu`, and the run's own readings come off the device. The library alone refuses the
//! same scene, by name, because it has no device and cannot acquire one — its workspace is thirteen
//! licence-gated crates that compile to `wasm32` and to Rust 1.78.
//!
//! Nothing here picks a device. That is the whole design: an accelerator is a lower-precision
//! computation, not a faster one, so choosing it is choosing what the run may lose — and a heuristic
//! that chose by grid size would have moved half the shipped scenes onto a different physics.

use pantometry_gpu::OnTheGpu;
use pantometry_world::{OnDisk, Scene, World};

fn scene(device: &str) -> Scene {
    let text = format!(
        r#"{{
      "title": "a block that says where it runs",
      "duration_s": 0.02,
      "frames": 3,
      "conservation_tolerance": 1e-4,
      "domains": [
        {{ "kind": "block", "name": "part", "cells": [12, 12, 12], "cell_mm": 1.0,
           "material": "aluminium", "initial_c": 20.0{device} }}
      ]
    }}"#
    );
    serde_json::from_str(&text).expect("the scene parses")
}

/// **The library refuses the device it does not have, and names what to do instead.**
///
/// Not a silent fall back. A run that quietly used the CPU when the scene said otherwise would be
/// answering a question nobody asked, and in single precision against double the two answers really
/// are different — which is the whole reason the device is stated.
#[test]
fn the_library_alone_refuses_and_says_why() {
    let err = match World::build(scene(", \"device\": \"gpu\"")) {
        Err(why) => why,
        Ok(_) => panic!("the library has no device and must not pretend otherwise"),
    };
    assert!(err.contains("part"), "name the domain: {err}");
    assert!(err.contains("gpu"), "name what was asked for: {err}");
    assert!(
        err.contains("build_with_accelerator"),
        "say what would honour it: {err}"
    );

    // And a scene that says nothing still runs, because `Cpu` is the default and every scene
    // written before this key existed means it.
    World::build(scene("")).expect("a scene with no device is a cpu scene");
}

/// **An application honours it, and the run comes off the device.**
#[test]
fn the_app_runs_the_block_where_the_scene_said() {
    let asked = scene(", \"device\": \"gpu\"");
    let mut world = match World::build_with_accelerator(asked, &OnDisk, &OnTheGpu) {
        Ok(w) => w,
        // No adapter is a skip, and it says so. A software rasteriser would be checking an
        // implementation nobody runs.
        Err(why) => return println!("skipped: {why}"),
    };
    let on_device = world.run().expect("the run holds its books on the device");

    // **The same scene on the CPU, and the two compared.** Counting frames would pass on a device
    // that computed nothing; what says the device ran the scene is that it got the scene's answer.
    let mut reference = World::build(scene("")).expect("the cpu scene builds");
    let on_cpu = reference.run().expect("the run holds its books on the cpu");
    assert_eq!(
        on_device.len(),
        on_cpu.len(),
        "the same scene, the same frames"
    );

    let last_gpu = &on_device.last().expect("a last frame").readings;
    let last_cpu = &on_cpu.last().expect("a last frame").readings;
    assert_eq!(last_gpu.len(), last_cpu.len(), "the same scalars");
    for (g, c) in last_gpu.iter().zip(last_cpu) {
        assert_eq!(g.label, c.label, "the same readings in the same order");
        assert!(
            g.value.is_finite(),
            "{}: not a number on the device",
            g.label
        );
        // Single precision against double, over the run: about `1e-7` a step walking as `√k`. A
        // thousandth of the value is that with room, and is far tighter than the difference a
        // *different operator* would produce — the uniform-coefficient kernel this replaced was
        // wrong by whole per cent on a block of two materials.
        let scale = c.value.abs().max(1e-9);
        assert!(
            (g.value - c.value).abs() / scale < 1e-3,
            "{}: device {} against cpu {}",
            g.label,
            g.value,
            c.value
        );
    }

    // And the audit did not stop either of them. `conservation_tolerance` is 1e-4 in this scene
    // because single precision cannot meet the 1e-9 default, and choosing that number is choosing
    // what the run is allowed to lose — which the scene does in writing.
}

/// **A block the device cannot run is refused through the accelerator too, with the reason.**
#[test]
fn a_cooled_block_on_the_device_is_refused_by_name() {
    let text = r#"{
      "title": "a block that sheds",
      "duration_s": 0.02,
      "frames": 3,
      "conservation_tolerance": 1e-4,
      "domains": [
        { "kind": "block", "name": "part", "cells": [8, 8, 8], "cell_mm": 1.0,
          "material": "aluminium", "initial_c": 200.0, "device": "gpu",
          "cooling": [{ "face": "z-max", "ambient_c": 20.0,
                        "convection_w_per_m2_k": 7.0, "area_cm2": 0.64 }] }
      ]
    }"#;
    let asked: Scene = serde_json::from_str(text).expect("the scene parses");
    match World::build_with_accelerator(asked, &OnDisk, &OnTheGpu).map(|_| ()) {
        Err(why) => {
            assert!(why.contains("part"), "name the domain: {why}");
            assert!(
                why.contains("film") || why.contains("adapter") || why.contains("device"),
                "say what it cannot do: {why}"
            );
        }
        Ok(_) => panic!("a film has no device pass and must not be accepted silently"),
    }
}
