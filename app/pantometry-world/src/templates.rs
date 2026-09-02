//! A minimal, valid example of every domain the scene format defines.
//!
//! The editor needs these to offer "add a domain" at all, and where they live is the whole
//! question: a list of nineteen things that has to stay in step with an enum nobody can enumerate
//! at runtime.
//!
//! # Why here, and how the list is kept honest
//!
//! Beside the format, because that is what they are facts about — and because the check that
//! keeps them in step reads the shipped scenes off a disk, which `editor-core` cannot do: it
//! compiles to `wasm32`, where there is no repository. That is the lesson `counts_in_prose`
//! learned when it lived in the wrong crate and five `test` jobs went red.
//!
//! Rust cannot list an enum's variants, so no compiler check can prove this list is complete.
//! What can be proved is that it agrees with a set maintained somewhere else entirely, and one
//! exists: **every kind the format defines appears in at least one shipped scene.** Measured, not
//! assumed — nineteen variants, nineteen distinct kinds across the thirty scenes. So the
//! test compares the two sets **in both directions**, which is the same shape as
//! `counts_in_prose`: a template with no scene fires one half, a twentieth domain with a scene
//! fires the other.
//!
//! # Text, not a serialised `DomainSpec`
//!
//! These are spliced into a file whose formatting is the thing the editor's whole write path
//! exists to preserve — see `editor_core::edit`. Serialising a value would insert an object with
//! alphabetised keys into a document written `kind` and `name` first, which is the reflow that
//! path refuses to cause. So they are written the way the scenes are written, and the test that
//! parses each one is what keeps them valid rather than merely well-indented.
//!
//! # Two of them name another domain, and cannot not
//!
//! `beam` states what it shines `onto` and `structure` states the block it `follows`. Both are
//! required, so there is no such thing as a template for them that stands alone: whatever name is
//! written here is a name the receiving scene probably does not have. The build says so by name,
//! which is the right place for that to be said — and `a_template_either_builds_or_names_the_
//! reference_it_needs` pins that those two are the only ones, so a third arriving unnoticed is a
//! test failure rather than a surprise in the editor.
/// One kind of domain, with everything needed to offer it and to start a scene from it.
///
/// Grown from a `(kind, json)` pair when the editor gained a chooser: a person picking what to
/// simulate needs to be told what each kind *is*, and a scene made of one needs a duration and a
/// frame count that the kind can actually run under. Those two numbers span fourteen orders of
/// magnitude across this table — `atoms` settles in picoseconds, a thermal `network` in half an
/// hour — which is a fact about the physics and the reason a chooser has to say so rather than
/// pick a default and hope.
pub struct Template {
    /// The `kind` the scene format spells.
    pub kind: &'static str,
    /// One line, for somebody choosing rather than somebody maintaining.
    pub about: &'static str,
    /// A duration this kind is known to run under, from the first shipped scene that uses it.
    pub duration_s: f64,
    /// The frame count from that same scene.
    pub frames: usize,
    /// The domain itself, as text. See the module documentation for why it is not a value.
    pub json: &'static str,
}

impl Template {
    /// Whether this kind names another domain it cannot do without.
    ///
    /// Two do: `beam` states what it shines `onto` and `structure` states the block it `follows`,
    /// and both keys are required. Derived from the template's own text rather than listed, so a
    /// third arriving does not need this to be remembered — and
    /// `exactly_two_kinds_name_a_partner` pins the set so a fourth key spelled differently is a
    /// test failure rather than a row in a chooser that quietly stops saying so.
    pub fn needs_a_partner(&self) -> bool {
        self.json.contains("\"onto\"") || self.json.contains("\"follows\"")
    }
}

/// The examples, keyed by the `kind` the scene format spells.
///
/// Sorted by kind, so a menu built from this is in a stable order.
pub const TEMPLATES: [Template; 19] = [
    Template {
        kind: "atoms",
        about: "A Lennard-Jones fluid in a periodic box",
        duration_s: 6e-12,
        frames: 12,
        json: r#"{ "kind": "atoms", "name": "atoms", "cells": 3, "density": 0.8442, "temperature": 1.4,
          "thermostat_t": 1.4, "seed": 20260808 }"#,
    },
    Template {
        kind: "bar",
        about: "Heat along a one-dimensional conducting bar",
        duration_s: 4.0,
        frames: 9,
        json: r#"{ "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61, "area_mm2": 100.0,
          "initial_c": 20.0 }"#,
    },
    Template {
        kind: "beam",
        about: "A beam that heats where it lands, over a shared boundary",
        duration_s: 0.2,
        frames: 13,
        json: r#"{ "kind": "beam", "name": "beam", "onto": "another domain in this scene", "faces": 61,
          "face_area_mm2": 0.5, "watts": 4.0, "reserve_j": 0.4, "waist_fraction": 0.08 }"#,
    },
    Template {
        kind: "block",
        about: "Heat through a three-dimensional conducting block",
        duration_s: 0.006,
        frames: 11,
        json: r#"{ "kind": "block", "name": "block", "cells": [11, 11, 11], "cell_mm": 1.0,
          "initial_c": -10.0, "material": "ice" }"#,
    },
    Template {
        kind: "bounce",
        about: "A ball bouncing on a floor, losing energy to its dashpot",
        duration_s: 1.0,
        frames: 24,
        json: r#"{ "kind": "bounce", "name": "bounce", "drop_m": 0.5, "mass_kg": 0.0113,
          "stiffness": 200000.0, "damping": 20.0 }"#,
    },
    Template {
        kind: "cavity",
        about: "Maxwell's equations on a Yee grid: fields in a resonant box",
        duration_s: 4e-9,
        frames: 200,
        json: r#"{ "kind": "cavity", "name": "cavity", "cells": [24, 4, 24], "cell_mm": 5.0,
          "medium": "vacuum", "mode": [1, 1], "amplitude_v_per_m": 1000.0 }"#,
    },
    Template {
        kind: "channel",
        about: "Incompressible flow on a staggered grid",
        duration_s: 4.0,
        frames: 16,
        json: r#"{ "kind": "channel", "name": "channel", "cells": [4, 16, 4], "cell_mm": 0.125,
          "fluid": "water", "walls": {"lower_m_per_s": 0.0, "upper_m_per_s": 0.0},
          "drive_m_per_s2": [0.02, 0.0, 0.0] }"#,
    },
    Template {
        kind: "conductor",
        about: "A block of conductor with two electrodes, solved for its potential",
        duration_s: 1.0,
        frames: 6,
        json: r#"{ "kind": "conductor", "name": "conductor", "cells": [12, 5, 5], "cell_mm": 1.0,
          "resistivity_ohm_m": 1.724e-08, "volts": 0.001,
          "blocked": [[6, 0, 0], [6, 0, 1], [6, 0, 2], [6, 0, 3], [6, 0, 4], [6, 1, 0], [6, 1, 1], [6, 1, 2], [6, 1, 3], [6, 1, 4], [6, 2, 0], [6, 2, 1], [6, 2, 2], [6, 2, 3], [6, 2, 4]] }"#,
    },
    Template {
        kind: "hall",
        about: "A room with a ceiling: the wave equation in three dimensions",
        duration_s: 0.02,
        frames: 11,
        json: r#"{ "kind": "hall", "name": "hall", "width_m": 4.4, "height_m": 3.1, "depth_m": 2.4,
          "nodes_across": 23, "mode": [1, 1, 1], "amplitude_pa": 1.0 }"#,
    },
    Template {
        kind: "heater",
        about: "A heat source with a finite tank of energy to spend",
        duration_s: 4.0,
        frames: 9,
        json: r#"{ "kind": "heater", "name": "heater", "watts": 2.0, "reserve_j": 6.0 }"#,
    },
    Template {
        kind: "light",
        about: "A lamp on a coated surface: real spectra deciding how much becomes heat",
        duration_s: 3.0,
        frames: 9,
        json: r#"{ "kind": "light", "name": "light", "watts": 40.0, "colour_k": 3200.0,
          "finish": "aluminium", "reserve_j": 2.0 }"#,
    },
    Template {
        kind: "lump",
        about: "One thermal mass at one temperature, losing heat to still air",
        duration_s: 1.0,
        frames: 24,
        json: r#"{ "kind": "lump", "name": "lump", "volume_cm3": 2.0, "thickness_mm": 6.0,
          "initial_c": 20.0, "ambient_c": 20.0, "area_cm2": 12.0 }"#,
    },
    Template {
        kind: "network",
        about: "Lumped bodies joined by conductances: junction, case, ambient",
        duration_s: 1800.0,
        frames: 7,
        json: r#"{ "kind": "network", "name": "network", "absorbing": "winding",
          "nodes": [{"name": "winding", "material": "copper", "volume_cm3": 18.0, "thickness_mm": 2.0, "initial_c": 25.0}, {"name": "stator", "material": "electrical_steel", "volume_cm3": 140.0, "thickness_mm": 8.0, "initial_c": 25.0}, {"name": "housing", "material": "aluminium", "volume_cm3": 220.0, "thickness_mm": 4.0, "initial_c": 25.0, "loses_to": {"ambient_c": 25.0, "area_cm2": 420.0}}],
          "links": [{"from": "winding", "to": "stator", "w_per_k": 0.9}, {"from": "stator", "to": "housing", "w_per_k": 2.4}] }"#,
    },
    Template {
        kind: "orbit",
        about: "Bodies under their own gravity: a central mass with satellites",
        duration_s: 7200.0,
        frames: 15,
        json: r#"{ "kind": "orbit", "name": "orbit", "central_kg": 1.989e+30,
          "radii_m": [57900000000.0, 108200000000.0, 149600000000.0],
          "inclinations_deg": [7.0, 3.4, 0.0], "satellite_kg": 3.3e+23 }"#,
    },
    Template {
        kind: "puck",
        about: "A basket of packed grounds with liquid driven through it",
        duration_s: 8.0,
        frames: 9,
        json: r#"{ "kind": "puck", "name": "puck", "cells": [19, 10, 19], "cell_mm": 2.0,
          "radius_mm": 15.0, "grind_um": 250.0, "porosity": 0.45, "bar": 9.0, "brew_c": 93.0 }"#,
    },
    Template {
        kind: "room",
        about: "A two-dimensional box of air with rigid walls, released in a standing mode",
        duration_s: 0.02,
        frames: 11,
        json: r#"{ "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1, "cells_across": 81,
          "release": {"as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0} }"#,
    },
    Template {
        kind: "structure",
        about: "A solid body under load, solved for its displacement",
        duration_s: 60.0,
        frames: 12,
        json: r#"{ "kind": "structure", "name": "structure", "cells": [8, 8, 8], "cell_mm": 1.5,
          "material": "copper",
          "held": [ { "face": "x-min", "as": "roller" }, { "face": "y-min", "as": "roller" },
                    { "face": "z-min", "as": "roller" } ],
          "follows": "another domain in this scene", "reference_c": 20.0 }"#,
    },
    Template {
        kind: "well",
        about: "One particle in a hard-walled well",
        duration_s: 2e-13,
        frames: 20,
        json: r#"{ "kind": "well", "name": "well", "cells": 200, "width_nm": 10.0,
          "start": {"eigenstate": 3} }"#,
    },
    Template {
        kind: "winding",
        about: "A copper winding dissipating I²R, which is where a motor's heat comes from",
        duration_s: 1800.0,
        frames: 7,
        json: r#"{ "kind": "winding", "name": "winding", "length_m": 62.0, "cross_section_mm2": 0.35,
          "amps": 1.75, "at_c": 90.0, "reserve_j": 200000.0 }"#,
    },
];

/// A scene made of one domain of each named kind, ready to open in the editor.
///
/// # The duration is the shortest of them, and that is a decision with a reason
///
/// A scene has **one** `duration_s`, and these span fourteen orders of magnitude: `atoms` settles
/// in picoseconds and a thermal `network` in half an hour. There is no duration at which both are
/// a simulation — one of them does nothing at any number you pick. So the shortest is taken,
/// together with the frame count from the same shipped scene, and the chooser says on screen when
/// the span is wide rather than quietly producing a scene in which half the domains are frozen.
///
/// # References are left as the templates write them
///
/// `beam` names what it shines `onto` and `structure` names the block it `follows`, and neither is
/// satisfied by pointing at another domain's *name*: a beam's `faces` has to equal the bar's
/// `cells` and its `onto` has to match the bar's `exposes`, which is three fields agreeing. Wiring
/// that here would be guessing. Left alone, the format says exactly what is missing — `build`
/// refuses `structure`, and the conservation audit stops `beam` at the first step naming the
/// domain that is not there — which is a better sentence than any this could invent.
pub fn scene(kinds: &[&str]) -> String {
    let chosen: Vec<&Template> = TEMPLATES
        .iter()
        .filter(|t| kinds.contains(&t.kind))
        .collect();
    let Some(shortest) = chosen
        .iter()
        .min_by(|a, b| a.duration_s.total_cmp(&b.duration_s))
    else {
        return String::from("{\n  \"title\": \"an empty scene\",\n  \"duration_s\": 1.0,\n  \"frames\": 2,\n  \"domains\": []\n}\n");
    };
    let names: Vec<&str> = chosen.iter().map(|t| t.kind).collect();
    let domains: Vec<String> = chosen
        .iter()
        .map(|t| format!("    {}", t.json.replace('\n', "\n    ")))
        .collect();
    format!(
        "{{\n  \"title\": \"a new scene: {}\",\n  \"schedule\": \"multirate\",\n  \"duration_s\": {},\n  \"frames\": {},\n  \"domains\": [\n{}\n  ]\n}}\n",
        names.join(", "),
        shortest.duration_s,
        shortest.frames,
        domains.join(",\n")
    )
}

/// How far apart the chosen kinds' natural durations are, as a ratio.
///
/// One is a set that runs together. `1e6` is a set where the slowest domain advances by a
/// millionth of what it needs while the fastest finishes — which is not an error the format can
/// refuse, because every one of those scenes is well formed. It is a thing to say out loud.
pub fn timescale_span(kinds: &[&str]) -> f64 {
    let secs: Vec<f64> = TEMPLATES
        .iter()
        .filter(|t| kinds.contains(&t.kind))
        .map(|t| t.duration_s)
        .collect();
    match (
        secs.iter().copied().reduce(f64::min),
        secs.iter().copied().reduce(f64::max),
    ) {
        (Some(lo), Some(hi)) if lo > 0.0 => hi / lo,
        _ => 1.0,
    }
}
