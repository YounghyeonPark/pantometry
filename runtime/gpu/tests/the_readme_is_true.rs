//! **The README's scene is parsed, not typed.**
//!
//! `README.md` shows the JSON that asks for a device. Its first draft put `conservation_tolerance`
//! inside the domain, where `Scene`'s `deny_unknown_fields` refuses it — and nothing would have
//! said so, because a README is not compiled. This parses the fenced block out of the file itself,
//! so the two cannot drift.

use pantometry_world::{Device, DomainSpec, Scene};

/// The first ```json fence in the README, as a `Scene`.
#[test]
fn the_readme_scene_parses() {
    let readme = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/README.md"))
        .expect("this crate's README is beside its manifest");

    // Take the fence rather than a line count: prose moves and a line number is a claim about
    // formatting.
    let after = readme
        .split_once("```json\n")
        .expect("the README shows a scene in a json fence")
        .1;
    let json = after.split_once("\n```").expect("the fence is closed").0;

    let scene: Scene = serde_json::from_str(json).unwrap_or_else(|why| {
        panic!("the README's scene does not parse: {why}\n--- what it says ---\n{json}")
    });

    // And it says the thing the section is about.
    let block = scene
        .domains
        .first()
        .expect("the scene has a domain to run");
    assert_eq!(
        block.device(),
        Device::Gpu,
        "the section is about asking for a device, so the example has to ask for one"
    );
    assert!(
        matches!(block, DomainSpec::Block { .. }),
        "a device is a block's to ask for"
    );

    // The tolerance is the point of the paragraph under it: the default cannot be met in f32.
    assert!(
        scene.conservation_tolerance > 1e-9,
        "a device scene that kept the default audit would fail it: {}",
        scene.conservation_tolerance
    );
}
