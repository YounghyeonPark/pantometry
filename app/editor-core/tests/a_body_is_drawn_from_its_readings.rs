//! **A compartment model is drawn as its body, and the drawing is the run.**
//!
//! Scene `34`, a tablet through the liver into lean tissue and fat, is the case. The diagram takes
//! its shape from the scene and every value from the run; each check here is against one of those
//! two, or against a closed form, and never against a second copy of the mapping:
//!
//! - **nothing is left out**: one disc per compartment and depot and one line per vessel, in every
//!   shipped scene that has a compartment model, and a disc's name is its compartment's;
//! - **a disc's area is its volume**, to rounding;
//! - **a disc's colour is the reading under its own label**, looked up here from the frame rather
//!   than through the diagram, on one scale for the whole run;
//! - and **the discs hold the body.** At every frame, the concentration each disc is coloured by
//!   times the volume the scene gives it, summed, is the `in the body` the domain reports on its
//!   own. That is the definition of the amount in a well-stirred compartment, `A = C·V`, and the
//!   domain sums the amounts rather than the products, so the two meet only if every disc read
//!   its own compartment.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped scenes off a disk.

#![cfg(not(target_family = "wasm"))]

use editor_core::diagram::{bodies, view};
use pantometry_world::{DomainSpec, Scene};
use viewer_core::Run;

const SCENES: &str = "../pantometry-world/scenes";
const TABLET: &str = "34-a-tablet-that-meets-the-liver-first.json";
const MOTOR: &str = "11-motor-thermal-network.json";

fn text(name: &str) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join(SCENES)
        .join(name);
    std::fs::read_to_string(path).expect("the shipped scene reads")
}

fn scene(name: &str) -> Scene {
    serde_json::from_str(&text(name)).expect("the shipped scene parses")
}

fn ran(name: &str) -> Run {
    let json = editor_core::run(&text(name), &editor_core::OnDisk).expect("the scene runs");
    Run::from_json(&json).expect("the run parses")
}

/// The reading under `domain` and `label` in frame `f`, found here and not by the diagram.
fn reading(run: &Run, f: usize, domain: &str, label: &str) -> f64 {
    run.frames[f]
        .readings
        .iter()
        .find(|r| r.domain == domain && r.label == label)
        .unwrap_or_else(|| panic!("frame {f} has no {domain} {label}"))
        .value
}

const RECT: [f64; 4] = [0.0, 0.0, 640.0, 420.0];

/// **Every compartment, vessel and depot of every shipped compartment model is drawn, once.**
#[test]
fn nothing_in_a_model_is_left_out() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(SCENES);
    let mut models = 0;
    for entry in std::fs::read_dir(&dir).expect("the scenes directory reads") {
        let path = entry.expect("an entry").path();
        if path.extension().and_then(|e| e.to_str()) != Some("json") {
            continue;
        }
        let s: Scene = serde_json::from_str(&std::fs::read_to_string(&path).expect("reads"))
            .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let laid = bodies(&s);
        let specs: Vec<_> = s
            .domains
            .iter()
            .filter_map(|d| match d {
                DomainSpec::Compartments {
                    name,
                    volumes,
                    links,
                    perfusions,
                    depots,
                    ..
                } => Some((name, volumes, links, perfusions, depots)),
                _ => None,
            })
            .collect();
        assert_eq!(laid.len(), specs.len(), "{}", path.display());
        for (b, (name, volumes, links, perfusions, depots)) in laid.iter().zip(&specs) {
            models += 1;
            assert_eq!(&b.domain, *name);
            let labels: Vec<&str> = b.organs.iter().map(|o| o.label.as_str()).collect();
            let want: Vec<&str> = volumes.iter().map(|c| c.label.as_str()).collect();
            assert_eq!(
                labels,
                want,
                "{}: the compartments, in order",
                path.display()
            );
            assert_eq!(b.vessels.len(), links.len() + perfusions.len());
            for (v, p) in b.vessels.iter().zip(perfusions.iter()) {
                assert_eq!(b.organs[v.from].label, p.blood);
                assert_eq!(b.organs[v.to].label, p.tissue);
                assert_eq!(v.flow_l_per_h, p.flow_l_per_h);
            }
            for (v, l) in b.vessels[perfusions.len()..].iter().zip(links.iter()) {
                assert_eq!(b.organs[v.from].label, l.between[0]);
                assert_eq!(b.organs[v.to].label, l.between[1]);
            }
            assert_eq!(b.depots.len(), depots.len());
            for (p, d) in b.depots.iter().zip(depots.iter()) {
                assert_eq!(p.label, d.label);
                assert_eq!(b.organs[p.into].label, d.into);
            }
        }
    }
    // Scenes 32 and 34. Fewer and this test has stopped looking at what it says it does.
    assert!(models >= 2, "only {models} compartment model(s) found");
}

/// **A disc for every compartment and depot, inside the drawing, its area its volume, its name
/// beside it.**
#[test]
fn a_disc_is_its_compartment() {
    let run = ran(TABLET);
    let b = &bodies(&scene(TABLET))[0];
    let last = run.frames.len() - 1;
    let d = b.draw(&run, last, RECT);
    assert_eq!(d.discs.len(), b.organs.len() + b.depots.len());
    for c in &d.discs {
        assert!(c.radius > 0.0);
        assert!(c.centre[0] - c.radius >= RECT[0] && c.centre[0] + c.radius <= RECT[2]);
        assert!(c.centre[1] - c.radius >= RECT[1] && c.centre[1] + c.radius <= RECT[3]);
    }
    // Area over volume is one number for every compartment; `r²/V` against the first's.
    let k = d.discs[0].radius * d.discs[0].radius / b.organs[0].volume_l;
    for (c, o) in d.discs.iter().zip(&b.organs) {
        let ratio = c.radius * c.radius / o.volume_l / k;
        assert!(
            (ratio - 1.0).abs() < 1e-12,
            "{}: area/volume is {ratio} of blood's",
            o.label
        );
    }
    for o in &b.organs {
        let named = d.texts.iter().filter(|t| t.text == o.label).count();
        assert_eq!(named, 1, "{} is named once", o.label);
    }
}

/// **A disc's colour is its own reading, on one scale for the whole run, and zero is an outline.**
#[test]
fn a_disc_is_coloured_by_its_own_reading() {
    let run = ran(TABLET);
    let b = &bodies(&scene(TABLET))[0];
    let (lo, hi) = b.scale(&run).expect("the run has positive concentrations");
    // The scale is whole decades holding every positive concentration in every frame.
    let all: Vec<f64> = (0..run.frames.len())
        .flat_map(|f| {
            b.organs
                .iter()
                .map(move |o| (f, o.label.clone()))
                .collect::<Vec<_>>()
        })
        .map(|(f, l)| reading(&run, f, &b.domain, &l))
        .filter(|v| *v > 0.0)
        .collect();
    let min = all.iter().copied().fold(f64::INFINITY, f64::min);
    let max = all.iter().copied().fold(0.0, f64::max);
    assert!(
        lo <= min && min < 10.0 * lo,
        "{lo} is the decade under {min}"
    );
    assert!(
        hi >= max && max > hi / 10.0,
        "{hi} is the decade over {max}"
    );
    let mut filled = 0;
    for f in 0..run.frames.len() {
        let d = b.draw(&run, f, RECT);
        for (c, o) in d.discs.iter().zip(&b.organs) {
            let v = reading(&run, f, &b.domain, &o.label);
            // The place on the scale worked out here, from the decades, and the colour straight
            // off the library's ramp — not through `Body::colour`, which is the thing under test.
            let want = (v > 0.0).then(|| {
                let u = (v.log10() - lo.log10()) / (hi.log10() - lo.log10());
                pantometry::view::ramp::sequential(u)
            });
            assert_eq!(c.fill, want, "frame {f} {} at {v}", o.label);
            if v > 0.0 {
                filled += 1;
                assert!(c.fill.is_some(), "frame {f} {} at {v} is filled", o.label);
            } else {
                assert!(
                    c.fill.is_none(),
                    "frame {f} {} at {v} is an outline",
                    o.label
                );
            }
        }
    }
    // Before the lag, nothing has reached any compartment; after it, every one is filled.
    let first = b.draw(&run, 0, RECT);
    assert!(
        first.discs.iter().all(|c| c.fill.is_none()),
        "nothing at t = 0"
    );
    assert!(filled > 4 * (run.frames.len() - 4), "{filled} filled discs");
}

/// **The discs hold the body.** `Σ Cᵢ·Vᵢ` over the concentrations the discs are coloured by and
/// the volumes the scene gives them is the `in the body` the domain reports, at every frame.
///
/// The bound is the file's. A run reaches an editor as JSON with seven significant figures, so
/// every reading is good to `5e-7` of itself: the sum is good to `5e-7·Σ|Cᵢ·Vᵢ|` and the reading
/// it is compared with to `5e-7` of that, and nothing else in the comparison rounds at that size.
#[test]
fn the_discs_hold_the_body() {
    let run = ran(TABLET);
    let b = &bodies(&scene(TABLET))[0];
    let mut worst: f64 = 0.0;
    let mut weighed = 0;
    for (f, frame) in run.frames.iter().enumerate() {
        let held: Vec<f64> = b
            .concentrations(frame)
            .iter()
            .zip(&b.organs)
            .map(|(c, o)| c.expect("every compartment has a reading") * o.volume_l)
            .collect();
        let sum: f64 = held.iter().sum();
        let body = reading(&run, f, &b.domain, "in the body");
        let bound = 5e-7 * (held.iter().map(|h| h.abs()).sum::<f64>() + body.abs());
        assert!(
            (sum - body).abs() <= bound,
            "frame {f}: the discs hold {sum} mg and the body {body} mg, {} apart against {bound}",
            (sum - body).abs()
        );
        if body > 1.0 {
            weighed += 1;
            worst = worst.max((sum - body).abs() / bound);
        }
    }
    // Most of four days has milligrams in the body; a test that weighed nothing proved nothing.
    assert!(
        weighed > run.frames.len() / 2,
        "only {weighed} frames weighed"
    );
    eprintln!("the discs hold the body to {worst:.3} of the bound over {weighed} frames");
}

/// **A run with no compartment model is drawn as its chart and nothing else**, the same drawing
/// the chart makes on its own.
#[test]
fn a_network_is_still_its_chart() {
    let run = ran(MOTOR);
    let laid = bodies(&scene(MOTOR));
    assert!(laid.is_empty());
    let frame = run.frames.len() / 2;
    let drawn = view(&laid, &run, false, 800.0, 600.0, frame).expect("a drawing");
    let chart = editor_core::chart::chart(&run, false)
        .expect("a chart")
        .draw(800.0, 600.0, Some(run.frames[frame].t));
    assert!(drawn.discs.is_empty());
    assert_eq!(drawn, chart);
}
