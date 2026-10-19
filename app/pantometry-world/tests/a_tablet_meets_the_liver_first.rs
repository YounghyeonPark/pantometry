//! Hepatic first pass, the partition of every tissue, and the steady state an infusion implies —
//! each from the scene's own numbers, and none against a second implementation.
//!
//! `34-a-tablet-that-meets-the-liver-first` is one blood pool with three tissues perfused off it:
//! a liver that eliminates (the well-stirred liver, `E = CL/(Q_h + CL)`), lean tissue and fat. A
//! tablet's depot empties into the **liver**, because there is no portal vein in a one-pool model
//! and that is how first pass is written. Renal clearance acts on the blood.
//!
//! # The identities, and why they are exact at any step
//!
//! Integrate every organ's balance from the dose to the end of time. Each amount starts and ends
//! at zero, so each derivative integrates to nothing and what is left is algebra in the areas
//! `A = ∫C dt`. For a tissue that does not eliminate, `0 = Q(A_b − A_t/K_p)`, so
//! **`A_t = K_p·A_b` for any dose by any route**. For the liver and the blood, a dose `D` into the
//! blood gives `D = A_b(CL_r + Q_h E)`, and the bioavailable `F·D` into the liver gives
//! `F·D = A_b(CL_r + Q_h E)/(1 − E)`. So
//!
//! ```text
//!   AUC_b(oral) / AUC_b(IV into blood) = F·(1 − E)          whatever CL_r is
//!   AUC_liver / AUC_b                  = K_p (1 − E)          after the IV dose
//!                                      = K_p (Q_h + CL_r)/Q_h after the oral one
//! ```
//!
//! The stepped model satisfies the same algebra **exactly**, not to first order: explicit Euler is
//! `x' = (I + hM)x + b`, and summing it over every step gives `h Σ xₙ = −M⁻¹(x₀ + Σ b)` — the
//! continuous area, at any `h` inside the stability limit. So these are checked at the shipped
//! 10 s and again at 40 s, and the tolerance is rounding, not truncation.
//!
//! # What the area's tail costs, bounded rather than guessed
//!
//! A run ends, and the area does not. What is left at the end is `A_rem` milligrams in the body and
//! the gut, and renal clearance alone removes `CL_r·C_b` of it each hour, so the blood's area still
//! to come is at most `A_rem/CL_r` — exactly, since the step eliminates `CL_r·C_b·h` and nothing
//! eliminates more than is there. The test asserts that bound is below 1e-15 of the area, rather
//! than assuming 240 hours is long enough.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped scene off a disk.

#![cfg(not(target_family = "wasm"))]

use pantometry::prelude::*;
use pantometry_world::{DomainSpec, DoseSpec, InfusionSpec, OnDisk, Scene, World};

const FILE: &str = "scenes/34-a-tablet-that-meets-the-liver-first.json";

fn shipped() -> Scene {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(FILE);
    let text = std::fs::read_to_string(&path).expect("the shipped scene reads");
    serde_json::from_str(&text).expect("the shipped scene parses")
}

/// One tissue's numbers, off the file: volume, blood flow, `K_p`, intrinsic clearance, in litres
/// and hours.
#[derive(Clone, Copy, Debug)]
struct Tissue {
    volume: f64,
    flow: f64,
    partition: f64,
    clearance: f64,
}

/// The body the scene states, read by label so a reordered file reads the same.
#[derive(Debug)]
struct Body {
    blood_volume: f64,
    renal: f64,
    liver: Tissue,
    lean: Tissue,
    fat: Tissue,
    bioavailability: f64,
    dose: f64,
}

impl Body {
    fn of(scene: &Scene) -> Body {
        let DomainSpec::Compartments {
            volumes,
            perfusions,
            depots,
            oral,
            ..
        } = &scene.domains[0]
        else {
            panic!("the scene's first domain is not the compartment model it was written to be");
        };
        let volume = |label: &str| {
            volumes
                .iter()
                .find(|c| c.label == label)
                .unwrap_or_else(|| panic!("no compartment called {label}"))
        };
        let tissue = |label: &str| {
            let p = perfusions
                .iter()
                .find(|p| p.tissue == label)
                .unwrap_or_else(|| panic!("{label} is not perfused"));
            assert_eq!(p.blood, "blood", "{label} is perfused from {}", p.blood);
            Tissue {
                volume: volume(label).volume_l,
                flow: p.flow_l_per_h,
                partition: p.partition,
                clearance: volume(label).clearance_l_per_h,
            }
        };
        assert_eq!(depots.len(), 1, "one route by mouth");
        assert_eq!(oral.len(), 1, "one tablet");
        assert_eq!(
            depots[0].into, "liver",
            "the tablet's depot has to feed the liver"
        );
        Body {
            blood_volume: volume("blood").volume_l,
            renal: volume("blood").clearance_l_per_h,
            liver: tissue("liver"),
            lean: tissue("lean tissue"),
            fat: tissue("fat"),
            bioavailability: depots[0].bioavailability,
            dose: oral[0].dose_mg,
        }
    }

    /// The well-stirred liver's extraction ratio.
    fn extraction(&self) -> f64 {
        self.liver.clearance / (self.liver.flow + self.liver.clearance)
    }
}

/// The scene with its tablet replaced by the same dose given intravenously, into the blood.
fn intravenous(mut scene: Scene, dose: f64) -> Scene {
    let DomainSpec::Compartments {
        boluses,
        depots,
        oral,
        ..
    } = &mut scene.domains[0]
    else {
        unreachable!("checked by Body::of");
    };
    depots.clear();
    oral.clear();
    boluses.push(DoseSpec {
        into: "blood".into(),
        dose_mg: dose,
    });
    scene
}

/// Neumaier's compensated sum, so the area is the model's and not the summation's.
#[derive(Default, Clone, Copy)]
struct Sum {
    sum: f64,
    carry: f64,
}

impl Sum {
    fn add(&mut self, x: f64) {
        let t = self.sum + x;
        self.carry += if self.sum.abs() >= x.abs() {
            (self.sum - t) + x
        } else {
            (x - t) + self.sum
        };
        self.sum = t;
    }
    fn value(self) -> f64 {
        self.sum + self.carry
    }
}

const ORGANS: [&str; 4] = ["blood", "liver", "lean tissue", "fat"];

/// The model inside a built world.
fn model(world: &World) -> &CompartmentModel {
    world
        .simulation()
        .domain_as::<CompartmentModel>("body")
        .expect("the scene's domain is a compartment model called body")
}

/// Every organ's area under its concentration curve, mg·h/L, by the left Riemann sum over every
/// step — which is the sum the identity above is about — run for `hours` at `window_s`. Returns the
/// areas, the step count, and what was left in the body and the gut at the end, in mg.
fn areas(scene: Scene, window_s: f64, hours: f64) -> ([f64; 4], usize, f64) {
    let mut scene = scene;
    scene.duration_s = hours * 3600.0;
    scene.window_s = Some(window_s);
    let steps = scene.steps();
    let dt = scene.duration_s / steps as f64;
    let mut world = World::build_with(scene, &OnDisk).expect("the scene builds");
    let handles: Vec<Compartment> = ORGANS
        .iter()
        .map(|l| model(&world).compartment_named(l).expect("an organ"))
        .collect();
    let mut sums = [Sum::default(); 4];
    let h = dt / 3600.0;
    for _ in 0..steps {
        let m = model(&world);
        for (s, c) in sums.iter_mut().zip(&handles) {
            // kg/m³ is g/L; the area is in mg·h/L.
            s.add(h * m.concentration(*c).to_si() * 1e3);
        }
        world
            .advance(Time::from_si(dt))
            .expect("the model steps inside its limit");
    }
    let m = model(&world);
    let left = (m.body_burden().to_si() + m.depot_mass().to_si() + m.pending_dose().to_si()) * 1e6;
    (sums.map(Sum::value), steps, left)
}

/// **The hepatic first pass: `AUC_b(oral)/AUC_b(IV) = F·(1 − E)`, at the shipped step and at four
/// times it, and every tissue's area is `K_p` times the blood's.**
///
/// Measured, as departures from the closed form, 240 hours at each step:
///
/// | ratio | closed form | 10 s, 86 400 steps | 40 s, 21 600 steps |
/// | --- | --- | --- | --- |
/// | `AUC_b` oral / IV | `F(1 − E)` = 0.54 | 4.4e-16 | 0 |
/// | liver / blood, IV | `K_p(1 − E)` = 1.2 | 0 | 0 |
/// | liver / blood, oral | `K_p(Q_h + CL_r)/Q_h` = 2.1333 | 2.2e-16 | 0 |
/// | lean / blood, both | `K_p` = 0.7 | 0, 0 | 0, 0 |
/// | fat / blood, IV and oral | `K_p` = 6 | 4.4e-16, 2.2e-16 | 1.1e-16, 2.2e-16 |
/// | `AUC_b` after IV | `D/(CL_r + Q_h E)` = 9.5238 mg·h/L | 2.2e-16 | 2.2e-16 |
///
/// The bound is `16ε` plus the tail. Not `n·ε`: a rounding made at one step is a perturbation the
/// stable map carries into the area once, as a relative error of its own size, so the roundings
/// of every step together are a few ε of the area and not `n` of them — eight operations at most
/// enter one amount's update, and the compensated sum and the ratio add two more. The tail is
/// bounded above, and is 3.0e-19 of the area at 240 hours.
///
/// What it would catch: the depot feeding the blood instead of the liver (the ratio is then `F`,
/// 0.9), the liver eliminating on `C_t` rather than on its venous `C_t/K_p` (an effective 120 L/h,
/// `E` = 0.571 and a ratio of 0.386), a perfusion's `K_p` swapped between two tissues or put on the
/// blood's side (the `K_p` ratios), and an IV reference that was not intravenous (1).
#[test]
fn the_liver_takes_its_share_before_the_blood_sees_the_tablet() {
    let scene = shipped();
    let body = Body::of(&scene);
    let e = body.extraction();
    let first_pass = body.bioavailability * (1.0 - e);
    println!(
        "  E = {e:.6}, F = {}, so F(1 - E) = {first_pass:.6}",
        body.bioavailability
    );

    for window in [10.0, 40.0] {
        let hours = 240.0;
        let (oral, n, oral_left) = areas(shipped(), window, hours);
        let (iv, _, iv_left) = areas(intravenous(shipped(), body.dose), window, hours);

        // The tail, bounded by renal clearance alone, in units of the area it would add to.
        let tail = (oral_left / body.renal / oral[0]).max(iv_left / body.renal / iv[0]);
        assert!(
            tail < 1e-15,
            "{window} s: {oral_left:.3e} and {iv_left:.3e} mg are left after {hours} h, which can \
             still add {tail:.2e} of the area"
        );
        println!("  {window:4.0} s: the area still to come is at most {tail:.1e} of it");
        let bound = 16.0 * f64::EPSILON + tail;

        let check = |what: &str, got: f64, want: f64| {
            let off = (got / want - 1.0).abs();
            println!(
                "  {window:4.0} s, {n:6} steps: {what:38} {got:.15} against {want:.15}, \
                 off {off:.2e} (bound {bound:.1e})"
            );
            assert!(
                off < bound,
                "{window} s: {what} is {got:.15} and the closed form is {want:.15}, {off:.2e} \
                 apart against a bound of {bound:.2e}"
            );
        };
        check(
            "AUC_b oral / AUC_b IV = F(1 - E)",
            oral[0] / iv[0],
            first_pass,
        );

        // Every tissue's area against the blood's, by route.
        let liver_after_iv = body.liver.partition * (1.0 - e);
        let liver_after_oral =
            body.liver.partition * (body.liver.flow + body.renal) / body.liver.flow;
        check(
            "IV: AUC_liver / AUC_b = K_p (1 - E)",
            iv[1] / iv[0],
            liver_after_iv,
        );
        check(
            "oral: AUC_liver / AUC_b = K_p (Q+CL_r)/Q",
            oral[1] / oral[0],
            liver_after_oral,
        );
        for (route, a) in [("IV", &iv), ("oral", &oral)] {
            check(
                &format!("{route}: AUC_lean / AUC_b = K_p"),
                a[2] / a[0],
                body.lean.partition,
            );
            check(
                &format!("{route}: AUC_fat / AUC_b = K_p"),
                a[3] / a[0],
                body.fat.partition,
            );
        }

        // And the blood's own area after the IV dose, which is the dose over the body's
        // clearance: renal plus what the liver extracts from the blood it is sent.
        check(
            "IV: AUC_b = D / (CL_r + Q_h E)",
            iv[0],
            body.dose / (body.renal + body.liver.flow * e),
        );
    }
}

/// **Started at the steady state an infusion implies, the body stays there: every tissue at `K_p`
/// times the blood, the liver at `K_p (1 − E)` times it.**
///
/// Under an infusion `R` into the blood the closed form is `C_b = R/(CL_r + Q_h E)`, and Euler's
/// fixed point `Mx* = −R` has no `h` in it, so a run started there does not move except by
/// rounding. The scene's tablet is replaced by that infusion, and the organs are dosed with their
/// steady-state amounts as boluses — all through the scene format, so the builder's mapping of
/// every perfusion is what is tested.
///
/// The bound is `(n + 8)·ε` relative: the start is the closed form to a few roundings (milligrams
/// to kilograms, litres to cubic metres, one division), and each step adds at most one rounding
/// that the stable map does not grow. What it would catch is any perfusion's `K_p` or flow mapped
/// to the wrong tissue, the liver eliminating on `C_t` rather than `C_t/K_p`, and renal clearance
/// lost: each moves the fixed point by a fraction of the whole, not by a rounding.
#[test]
fn an_infusion_holds_every_organ_at_its_partition() {
    let body = Body::of(&shipped());
    let e = body.extraction();
    // A round blood concentration, 1 mg/L, and the rate that holds it.
    let blood = 1.0;
    let rate = blood * (body.renal + body.liver.flow * e);
    let want = [
        blood,
        body.liver.partition * blood * (1.0 - e),
        body.lean.partition * blood,
        body.fat.partition * blood,
    ];
    let volumes = [
        body.blood_volume,
        body.liver.volume,
        body.lean.volume,
        body.fat.volume,
    ];
    println!(
        "  R = {rate} mg/h holds blood at {blood} mg/L; liver {:.6}, lean {:.6}, fat {:.6}",
        want[1], want[2], want[3]
    );

    let mut scene = shipped();
    let hours = 24.0;
    let DomainSpec::Compartments {
        boluses,
        infusions,
        depots,
        oral,
        ..
    } = &mut scene.domains[0]
    else {
        unreachable!("checked by Body::of");
    };
    depots.clear();
    oral.clear();
    for (organ, (c, v)) in ORGANS.iter().zip(want.iter().zip(volumes)) {
        boluses.push(DoseSpec {
            into: (*organ).into(),
            dose_mg: c * v,
        });
    }
    // Running past the end of the run, so no step takes the syringe's remainder: the last slice
    // of an infusion is `total − delivered`, which carries every earlier slice's rounding at once
    // and is a property of the bookkeeping rather than of the fixed point.
    infusions.push(InfusionSpec {
        into: "blood".into(),
        dose_mg: rate * 2.0 * hours,
        from_h: 0.0,
        to_h: 2.0 * hours,
    });
    scene.duration_s = hours * 3600.0;
    let steps = scene.steps();
    let dt = scene.duration_s / steps as f64;
    let mut world = World::build_with(scene, &OnDisk).expect("the held scene builds");
    let handles: Vec<Compartment> = ORGANS
        .iter()
        .map(|l| model(&world).compartment_named(l).expect("an organ"))
        .collect();

    let mut worst = [0.0f64; 4];
    for _ in 0..steps {
        world
            .advance(Time::from_si(dt))
            .expect("the model steps inside its limit");
        let m = model(&world);
        for (k, c) in handles.iter().enumerate() {
            let got = m.concentration(*c).to_si() * 1e3;
            worst[k] = worst[k].max((got / want[k] - 1.0).abs());
        }
    }
    // `(n + 8)ε`: the start is the closed form to a few roundings, and each step adds at most one
    // that the stable map does not grow — the bound the library's own steady-state test holds.
    //
    // **And not the clock's**, which would be a second term of `n·ε` and is zero here. A slice is
    // `(now + h) − now`, which rounds to half an ulp of `now` and makes the rate the blood is fed
    // at wobble by up to `ε·t/h` of itself — unless every `now` is an integer, as it is at a 10 s
    // step, where the subtraction is exact. Asserted, because this bound depends on it.
    assert_eq!(dt, 10.0, "the shipped step lands every `now` on an integer");
    let bound = (steps as f64 + 8.0) * f64::EPSILON;
    for (organ, w) in ORGANS.iter().zip(worst) {
        println!(
            "  {organ:12} worst departure over {steps} steps: {w:.2e} = {:.1} eps, bound {bound:.1e}",
            w / f64::EPSILON
        );
        assert!(
            w < bound,
            "{organ} left the steady state its partition implies by {w:.2e}, against a rounding \
             bound of {bound:.2e}"
        );
    }
}

/// A small body written out in full, with `extra` spliced into the domain, so each refusal below
/// is one key away from a scene that builds.
fn small(extra: &str) -> String {
    format!(
        r#"{{
  "title": "a tablet in a small body",
  "duration_s": 3600.0,
  "frames": 4,
  "window_s": 10.0,
  "domains": [
    {{
      "kind": "compartments",
      "name": "body",
      "volumes": [
        {{ "label": "blood", "volume_l": 5.0, "clearance_l_per_h": 6.0 }},
        {{ "label": "liver", "volume_l": 1.5, "clearance_l_per_h": 60.0 }}
      ]{extra}
    }}
  ]
}}"#
    )
}

fn built(extra: &str) -> Result<World, String> {
    let scene: Scene = serde_json::from_str(&small(extra)).map_err(|e| format!("parse: {e}"))?;
    World::build(scene)
}

fn refused(extra: &str) -> String {
    match built(extra) {
        Ok(_) => panic!("built a scene that should have been refused: {extra}"),
        Err(e) => e,
    }
}

const PERFUSED: &str = r#",
      "perfusions": [ { "blood": "blood", "tissue": "liver", "flow_l_per_h": 90.0, "partition": 2.0 } ]"#;
const GUT: &str = r#",
      "depots": [ { "label": "gut", "into": "liver", "absorption_per_h": 1.0, "bioavailability": 0.9 } ]"#;
const TABLET: &str = r#",
      "oral": [ { "depot": "gut", "dose_mg": 400.0, "at_h": 0.0 } ]"#;

/// **A scene dosed only by mouth is administered, its milligrams are in the note, and each new
/// key's mistakes are refused by the field that made them.**
///
/// The guard `nothing was administered` used to look at boluses and infusions alone, so the first
/// scene with only a tablet would have been refused as saying nothing; the note's sum would have
/// read 0.0 mg beside a 400 mg tablet. Every refusal names the field and the file's own number,
/// because the library's message carries SI — 90 L/h arrives there as 2.5e-5 m³/s.
#[test]
fn a_tablet_alone_is_a_dose_and_every_new_key_says_what_it_refused() {
    let world = built(&format!("{PERFUSED}{GUT}{TABLET}")).expect("a tablet alone builds");
    let notes = world.notes().join("\n");
    println!("  {notes}");
    assert!(
        notes.contains("400.0 mg administered, 400.0 mg of it by mouth"),
        "the note does not count the tablet: {notes}"
    );

    let cases: [(String, &[&str]); 10] = [
        // A perfusion's own refusals, by index, with the file's numbers.
        (
            PERFUSED.replace("90.0", "-90.0") + GUT + TABLET,
            &["perfusions[0]", "\"liver\" from \"blood\"", "flow_l_per_h -90"],
        ),
        (
            PERFUSED.replace("2.0 }", "0.0 }") + GUT + TABLET,
            &["perfusions[0]", "partition 0", "partition coefficient"],
        ),
        (
            PERFUSED.replace("\"tissue\": \"liver\"", "\"tissue\": \"kidney\"") + GUT + TABLET,
            &["no compartment called \"kidney\"", "in perfusions[0].tissue"],
        ),
        // A depot's.
        (
            format!("{PERFUSED}{}{TABLET}", GUT.replace("\"into\": \"liver\"", "\"into\": \"portal vein\"")),
            &["no compartment called \"portal vein\"", "in depots[0].into"],
        ),
        (
            format!("{PERFUSED}{}{TABLET}", GUT.replace("0.9", "1.5")),
            &["depots[0] \"gut\"", "bioavailability 1.5", "(0, 1]"],
        ),
        (
            format!("{PERFUSED}{}{TABLET}", GUT.replace("\"absorption_per_h\": 1.0", "\"absorption_per_h\": 0.0")),
            &["depots[0] \"gut\"", "absorption_per_h 0"],
        ),
        (
            format!(
                "{PERFUSED}{}{TABLET}",
                GUT.replace(
                    "\"bioavailability\": 0.9 }",
                    "\"bioavailability\": 0.9 }, { \"label\": \"gut\", \"into\": \"blood\", \"absorption_per_h\": 1.0, \"bioavailability\": 1.0 }"
                )
            ),
            &["depots[1]", "two depots are called \"gut\""],
        ),
        // A dose's, and a depot nothing goes through.
        (
            format!("{PERFUSED}{GUT}{}", TABLET.replace("\"depot\": \"gut\"", "\"depot\": \"stomach\"")),
            &["oral[0]", "no depot called \"stomach\"", "\"gut\""],
        ),
        (
            format!("{PERFUSED}{GUT}{}", TABLET.replace("400.0", "-400.0")),
            &["oral[0] through \"gut\"", "dose_mg -400"],
        ),
        (
            format!("{PERFUSED}{GUT}"),
            &["the depot \"gut\" has no dose in `oral`"],
        ),
    ];
    for (extra, wanted) in &cases {
        let err = refused(extra);
        println!("  refused: {err}");
        for w in *wanted {
            assert!(err.contains(w), "{err:?} does not say {w:?}");
        }
    }

    // A body with organs and a bolus, and no tablet, says it gave nothing by mouth — not `-0.0`.
    let world = built(&format!(
        "{PERFUSED},
      \"boluses\": [ {{ \"into\": \"blood\", \"dose_mg\": 100.0 }} ]"
    ))
    .expect("a bolus into a perfused body builds");
    let notes = world.notes().join(
        "
",
    );
    assert!(
        notes.contains("100.0 mg administered, 0.0 mg of it by mouth") && !notes.contains("-0.0"),
        "{notes}"
    );

    // And nothing at all is still nothing.
    let err = refused(PERFUSED);
    assert!(err.contains("nothing was administered"), "{err}");
    // A key the format does not have is refused, not dropped, inside each new entry.
    let err = refused(&format!(
        "{PERFUSED}{GUT}{}",
        TABLET.replace("\"at_h\"", "\"at_hour\"")
    ));
    assert!(
        err.starts_with("parse:") && err.contains("at_hour"),
        "{err}"
    );
}
