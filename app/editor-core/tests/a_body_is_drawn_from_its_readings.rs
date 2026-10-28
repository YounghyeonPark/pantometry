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

const BOLUS: &str = "32-a-dose-distributing-and-leaving.json";

/// The label the domain documents for each vessel's flow, read off the scene rather than the
/// diagram: `"blood to tissue"` for a perfusion and `"a to b"` in the declared order for a link.
/// Scenes 32 and 34 have no link beside a perfusion, so each vessel has its own.
fn flow_labels(s: &Scene) -> Vec<String> {
    s.domains
        .iter()
        .flat_map(|d| match d {
            DomainSpec::Compartments {
                perfusions, links, ..
            } => perfusions
                .iter()
                .map(|p| format!("{} to {}", p.blood, p.tissue))
                .chain(
                    links
                        .iter()
                        .map(|l| format!("{} to {}", l.between[0], l.between[1])),
                )
                .collect(),
            _ => Vec::new(),
        })
        .collect()
}

/// The heads of a drawing in the order it draws them: each three-point stroke in
/// [`editor_core::diagram::FLOW`], as its middle, the unit direction it points, and its length
/// from middle to tip.
fn heads(d: &editor_core::chart::Drawing) -> Vec<([f64; 2], [f64; 2], f64)> {
    d.strokes
        .iter()
        .filter(|s| s.colour == editor_core::diagram::FLOW && s.points.len() == 3)
        .map(|s| {
            let back = [
                0.5 * (s.points[0][0] + s.points[2][0]),
                0.5 * (s.points[0][1] + s.points[2][1]),
            ];
            let tip = s.points[1];
            let v = [tip[0] - back[0], tip[1] - back[1]];
            let n = v[0].hypot(v[1]);
            let mid = [0.5 * (tip[0] + back[0]), 0.5 * (tip[1] + back[1])];
            (mid, [v[0] / n, v[1] / n], 0.5 * n)
        })
        .collect()
}

/// **A head on every vessel that carries drug, on the vessel, pointing the way the domain says
/// it goes, and longer for more.**
///
/// At every frame of both compartment scenes: one head per vessel whose reading is not zero, in
/// the vessels' order; parallel to the line between its two discs and pointing from the one the
/// reading's label names first exactly when the reading is positive; its middle halfway between
/// the two rims. Then the length, over the whole run at once: a larger flow never has a shorter
/// head, in any two frames, which a scale renormalised per frame would break.
#[test]
fn a_vessel_points_the_way_its_drug_goes() {
    for name in [BOLUS, TABLET] {
        let s = scene(name);
        let run = ran(name);
        let b = &bodies(&s)[0];
        let labels = flow_labels(&s);
        assert_eq!(labels.len(), b.vessels.len());
        let mut sized: Vec<(f64, f64)> = Vec::new();
        for f in 0..run.frames.len() {
            let d = b.draw(&run, f, RECT);
            let flows: Vec<f64> = labels
                .iter()
                .map(|l| reading(&run, f, &b.domain, l))
                .collect();
            let carrying: Vec<usize> = (0..flows.len()).filter(|&k| flows[k] != 0.0).collect();
            let drawn = heads(&d);
            assert_eq!(drawn.len(), carrying.len(), "{name} frame {f}: {flows:?}");
            for (&k, &(mid, dir, len)) in carrying.iter().zip(&drawn) {
                let v = &b.vessels[k];
                let (from, to) = (&d.discs[v.from], &d.discs[v.to]);
                let along = [to.centre[0] - from.centre[0], to.centre[1] - from.centre[1]];
                let l = along[0].hypot(along[1]);
                let cross = dir[0] * along[1] - dir[1] * along[0];
                let dot = dir[0] * along[0] + dir[1] * along[1];
                assert!(
                    cross.abs() <= 1e-9 * l,
                    "{name} frame {f}: {} is off its vessel",
                    labels[k]
                );
                assert_eq!(
                    dot > 0.0,
                    flows[k] > 0.0,
                    "{name} frame {f}: {} reads {} and its head points the other way",
                    labels[k],
                    flows[k]
                );
                let halfway = from.radius + 0.5 * (l - from.radius - to.radius);
                let at = [mid[0] - from.centre[0], mid[1] - from.centre[1]];
                assert!(
                    (at[0].hypot(at[1]) - halfway).abs() <= 1e-9 * l,
                    "{name} frame {f}: {} is not halfway between the rims",
                    labels[k]
                );
                sized.push((flows[k].abs(), len));
            }
        }
        sized.sort_by(|a, b| a.0.total_cmp(&b.0));
        for w in sized.windows(2) {
            assert!(
                w[1].1 >= w[0].1 - 1e-9,
                "{name}: {} mg/h has a head of {} px and {} mg/h one of {}",
                w[0].0,
                w[0].1,
                w[1].0,
                w[1].1
            );
        }
        let (short, long) = (sized[0].1, sized[sized.len() - 1].1);
        assert!(long > short, "{name}: every head is {short} px");
        eprintln!("{name}: {} heads, {short:.2} to {long:.2} px", sized.len());
    }
}

/// **A tissue peaks the frame its head turns.** A compartment joined to nothing but one vessel,
/// clearing nothing and fed by no depot, changes only by what crosses that vessel; so it rises
/// while the head points in and falls once it points out, and its largest reading is at the last
/// frame its head points in or the first it does not. Those two frames bracket the step it turned
/// at, and nothing finer is in the run.
///
/// Read from the drawing, which way each head points frame by frame, against the tissue's own
/// concentration found here. Three tissues: scene 32's, and scene 34's lean tissue and fat.
#[test]
fn a_tissue_peaks_the_frame_its_head_turns() {
    let mut tissues = 0;
    for name in [BOLUS, TABLET] {
        let run = ran(name);
        let b = &bodies(&scene(name))[0];
        for (k, v) in b.vessels.iter().enumerate() {
            let t = &b.organs[v.to];
            let joined = b
                .vessels
                .iter()
                .filter(|w| w.from == v.to || w.to == v.to)
                .count();
            if joined != 1 || t.clearance_l_per_h > 0.0 || b.depots.iter().any(|p| p.into == v.to) {
                continue;
            }
            // Which way vessel `k`'s head points in each frame: in, out, or not drawn. Its place
            // among the heads is its place among the vessels that carry something.
            let inward: Vec<Option<bool>> = (0..run.frames.len())
                .map(|f| {
                    let d = b.draw(&run, f, RECT);
                    let flows = b.flows(&run.frames[f]);
                    let at = (0..b.vessels.len())
                        .filter(|&j| flows[j].is_some_and(|x| x != 0.0))
                        .position(|j| j == k)?;
                    let (_, dir, _) = heads(&d)[at];
                    let (from, to) = (&d.discs[v.from], &d.discs[v.to]);
                    let along = [to.centre[0] - from.centre[0], to.centre[1] - from.centre[1]];
                    Some(dir[0] * along[0] + dir[1] * along[1] > 0.0)
                })
                .collect();
            let filling = inward
                .iter()
                .position(|x| *x == Some(true))
                .unwrap_or_else(|| panic!("{name}: {} never fills", t.label));
            let turned = filling
                + inward[filling..]
                    .iter()
                    .position(|x| *x != Some(true))
                    .unwrap_or_else(|| panic!("{name}: {} never turns", t.label));
            let c: Vec<f64> = (0..run.frames.len())
                .map(|f| reading(&run, f, &b.domain, &t.label))
                .collect();
            let peak = (0..c.len()).max_by(|&a, &b| c[a].total_cmp(&c[b])).unwrap();
            eprintln!(
                "{name}: {} peaks at frame {peak} and its head turns at {turned}",
                t.label
            );
            // As values, not indices: seven figures keep the order of two readings and can make
            // them equal, and a tie would move an index without moving the claim.
            let bracket = c[turned - 1].max(c[turned]);
            assert_eq!(
                bracket.to_bits(),
                c[peak].to_bits(),
                "{name}: {} peaks at frame {peak}, and its head turns at frame {turned}",
                t.label
            );
            tissues += 1;
        }
    }
    assert_eq!(
        tissues, 3,
        "scene 32's tissue and scene 34's lean tissue and fat"
    );
}

/// **A link beside a perfusion is one head, pointing the way the two together carry.** No shipped
/// scene has one, so this is scene 34 with a 3 L/h link added between the blood and the lean
/// tissue: four vessels on three pairs. The domain reads the pair once, under the perfusion's
/// label and as the sum of both, so at every frame there are as many heads as pairs that carry
/// drug, and the lean tissue's points the way that one reading says.
#[test]
fn a_link_beside_a_perfusion_is_one_head() {
    let mut json: serde_json::Value = serde_json::from_str(&text(TABLET)).expect("parses");
    json["domains"][0]["links"] = serde_json::json!([
        { "between": ["blood", "lean tissue"], "clearance_l_per_h": 3.0 }
    ]);
    let text = json.to_string();
    let s: Scene = serde_json::from_str(&text).expect("still a scene");
    let run = Run::from_json(&editor_core::run(&text, &editor_core::OnDisk).expect("it runs"))
        .expect("the run parses");
    let b = &bodies(&s)[0];
    assert_eq!(b.vessels.len(), 4, "three perfusions and the link");
    let pairs = ["blood to liver", "blood to lean tissue", "blood to fat"];
    let lean = 1;
    let mut weighed = 0;
    for f in 0..run.frames.len() {
        let flows: Vec<f64> = pairs
            .iter()
            .map(|l| reading(&run, f, &b.domain, l))
            .collect();
        assert!(
            run.frames[f]
                .readings
                .iter()
                .filter(|r| r.unit == "mg/h")
                .count()
                == pairs.len(),
            "frame {f}: one reading per pair"
        );
        let d = b.draw(&run, f, RECT);
        let drawn = heads(&d);
        let carrying = flows.iter().filter(|j| **j != 0.0).count();
        assert_eq!(drawn.len(), carrying, "frame {f}: {flows:?}");
        if carrying == pairs.len() {
            let v = &b.vessels[lean];
            let (from, to) = (&d.discs[v.from], &d.discs[v.to]);
            let along = [to.centre[0] - from.centre[0], to.centre[1] - from.centre[1]];
            let (_, dir, _) = drawn[lean];
            assert_eq!(
                dir[0] * along[0] + dir[1] * along[1] > 0.0,
                flows[lean] > 0.0,
                "frame {f}: the lean tissue reads {} and its head points the other way",
                flows[lean]
            );
            weighed += 1;
        }
    }
    assert!(
        weighed > run.frames.len() / 2,
        "only {weighed} frames had every pair carrying"
    );
}
