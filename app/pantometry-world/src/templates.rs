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
//! assumed — twenty-two variants, twenty-two distinct kinds across the thirty-four scenes. So the
//! test compares the two sets **in both directions**, which is the same shape as
//! `counts_in_prose`: a template with no scene fires one half, a twenty-third domain with a scene
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
/// frame count that the kind can actually run under. Those two numbers span **seventeen** orders
/// of magnitude across this table — `well` settles in 2e-13 s and a `compartments` model runs for
/// twelve hours, a span of 2.16e17 — which is a fact about the physics and the reason a chooser has
/// to say so rather than pick a default and hope. It was 3.6e16 and `orbit` held the far end until
/// a dose took half a day to leave a body.
///
/// **It said fourteen orders and named `atoms` against a thermal `network`**, and both halves were
/// wrong. The extremes have always been `well` and `orbit`, so the count was stale whatever
/// `atoms` held; and `atoms` itself was 6e-12 against a domain whose time unit is one second,
/// which is the units error the doc below records. Measured from this table rather than
/// remembered: `well` 2e-13, `orbit` 7200.
///
/// **`atoms` said picoseconds and asked for 6e-12 s, and that was a units error.** The domain is
/// built with `LennardJones::reduced()`, where σ = ε = m = 1 in SI and the only time scale is
/// `τ = σ√(m/ε) = 1 second`. Six picoseconds is 6e-12 τ, ten orders of magnitude short of the
/// 0.01 τ the domain suggests for one step: measured, the mean square displacement grew as `t²`
/// with a ratio of exactly **4.000** between the last frame and the half-way one, which is free
/// flight — not one collision in the whole run. Argon is where the picoseconds came from: at
/// σ = 3.405 Å, ε/k = 119.8 K and m = 39.948 u, τ = 2.16 ps and six of those is a real run. The
/// number here is in the units the domain actually uses.
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
pub const TEMPLATES: [Template; 22] = [
    Template {
        kind: "atoms",
        about: "A Lennard-Jones fluid in a periodic box",
        duration_s: 6.0,
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
        frames: 52,
        json: r#"{ "kind": "beam", "name": "beam", "onto": "another domain in this scene", "faces": 61,
          "face_area_mm2": 0.5, "watts": 4.0, "reserve_j": 0.4, "waist_fraction": 0.08 }"#,
    },
    Template {
        kind: "block",
        about: "Heat through a three-dimensional conducting block",
        duration_s: 0.006,
        frames: 88,
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
        frames: 1600,
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
        kind: "compartments",
        about: "A dose distributing between compartments and leaving the body",
        duration_s: 43200.0,
        frames: 13,
        json: r#"{ "kind": "compartments", "name": "compartments",
          "volumes": [ { "label": "plasma", "volume_l": 5.0, "clearance_l_per_h": 3.0 },
                       { "label": "tissue", "volume_l": 20.0 } ],
          "links": [ { "between": ["plasma", "tissue"], "clearance_l_per_h": 8.0 } ],
          "boluses": [ { "into": "plasma", "dose_mg": 500.0 } ] }"#,
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
        kind: "molecule",
        about: "A small molecule, every atom, relaxing under a force field",
        // One minimiser iteration a frame, and aspirin's minimum takes 454 of them: a template
        // that stopped short of that would open on a molecule that has not finished relaxing.
        // The duration is a count of iterations, because a step of this kind is not a time.
        duration_s: 500.0,
        frames: 500,
        json: r#"{ "kind": "molecule", "name": "molecule", "ccd": "structures/AIN.cif" }"#,
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
        kind: "protein",
        about: "A protein moving along the collective modes of its own fold",
        duration_s: 2e-11,
        frames: 16,
        json: r#"{ "kind": "protein", "name": "protein", "pdb": "structures/1CRN.pdb",
          "cutoff_a": 15.0, "stiffness_n_per_m": 0.695, "temperature_c": 37.0,
          "seed": 20260914 }"#,
    },
    Template {
        kind: "puck",
        // **Eight seconds was a quarter of a shot**, and this suggested it because
        // `18-an-espresso-shot` ran it, which is what these numbers are taken from. That scene is
        // 25 s now — 19.90% extraction, the middle of what a barista pulls to — and a starting
        // point should open on a shot. `the_schedule_a_template_suggests_is_one_a_scene_uses`
        // caught this: the template stayed at 8 s after the scene moved, which is exactly the
        // silence that test exists for.
        about: "A basket of packed grounds with liquid driven through it",
        duration_s: 25.0,
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
/// A scene has **one** `duration_s`, and these span **sixteen** orders of magnitude: `well` settles
/// in 2e-13 s and an `orbit` takes 7200. There is no duration at which both are
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

/// One entry of a list a domain holds — a cooled face, a heat source, a region, a dose — as text.
///
/// The editor's model tree shows each list a domain states and offers to add to it, and an entry
/// cannot be made up: `deny_unknown_fields` refuses a key the entry does not take, and a missing
/// required one is refused as well. So the starting point lives here beside the format, for the
/// reason [`TEMPLATES`] does, and is text for the same reason too: it is spliced into a file whose
/// formatting the editor's whole write path exists to keep.
///
/// Each is one line, so it can be put into a list at any indentation without reflowing anything.
pub struct Entry {
    /// The `kind` of the domain that holds the list.
    pub kind: &'static str,
    /// The list's key, as the format spells it.
    pub key: &'static str,
    /// The entry itself, as text.
    pub json: &'static str,
}

/// The entries, sorted by kind and then key.
///
/// **Every list in the format, not only the ones a scene uses.** Seventeen lists of objects across
/// four kinds of domain, which is what the format's `Vec<…Spec>` fields come to; the shipped scenes
/// use fifteen of them, and `infusions` and `pressed` appear in none. A test holds the shipped
/// scenes' lists against this table, and another adds each entry to its kind's template and builds
/// the result.
///
/// The names an entry refers to — a compartment's `into`, a link's `from` and `to` — are the
/// template's own, so each builds there. In another scene they may name nothing, and the build
/// says which name, which is the sentence a person can act on. What an entry says about *where in
/// a grid* is different, and [`entry_for`] fits it to the grid the entry lands in.
pub const ENTRIES: [Entry; 17] = [
    Entry {
        kind: "block",
        key: "contact",
        json: r#"{ "axis": "z", "at": 5, "w_per_m2_k": 5000.0 }"#,
    },
    Entry {
        kind: "block",
        key: "cooling",
        json: r#"{ "face": "z-max", "ambient_c": 20.0, "convection_w_per_m2_k": 10.0, "area_cm2": 1.21 }"#,
    },
    Entry {
        kind: "block",
        key: "dissipation",
        json: r#"{ "watts": 1.0, "from": [3, 3, 3], "to": [7, 7, 7] }"#,
    },
    Entry {
        kind: "block",
        key: "parts",
        json: r#"{ "stl": "parts/rod.stl", "material": "aluminium" }"#,
    },
    Entry {
        kind: "block",
        key: "regions",
        json: r#"{ "material": "aluminium", "from": [0, 0, 0], "to": [5, 5, 5] }"#,
    },
    Entry {
        kind: "compartments",
        key: "boluses",
        json: r#"{ "into": "plasma", "dose_mg": 100.0 }"#,
    },
    Entry {
        kind: "compartments",
        key: "depots",
        json: r#"{ "label": "gut", "into": "plasma", "absorption_per_h": 1.0, "bioavailability": 0.8 }"#,
    },
    Entry {
        kind: "compartments",
        key: "infusions",
        json: r#"{ "into": "plasma", "dose_mg": 100.0, "from_h": 0.0, "to_h": 1.0 }"#,
    },
    Entry {
        kind: "compartments",
        key: "links",
        json: r#"{ "between": ["plasma", "tissue"], "clearance_l_per_h": 1.0 }"#,
    },
    Entry {
        kind: "compartments",
        key: "oral",
        json: r#"{ "depot": "gut", "dose_mg": 100.0, "at_h": 0.0 }"#,
    },
    Entry {
        kind: "compartments",
        key: "perfusions",
        json: r#"{ "blood": "plasma", "tissue": "tissue", "flow_l_per_h": 10.0, "partition": 1.0 }"#,
    },
    Entry {
        kind: "compartments",
        key: "volumes",
        json: r#"{ "label": "compartment", "volume_l": 1.0 }"#,
    },
    Entry {
        kind: "network",
        key: "links",
        json: r#"{ "from": "winding", "to": "housing", "w_per_k": 0.5 }"#,
    },
    Entry {
        kind: "network",
        key: "nodes",
        json: r#"{ "name": "node", "material": "aluminium", "volume_cm3": 10.0, "thickness_mm": 2.0, "initial_c": 25.0 }"#,
    },
    Entry {
        kind: "structure",
        key: "held",
        json: r#"{ "face": "z-max", "as": "clamp" }"#,
    },
    Entry {
        kind: "structure",
        key: "pressed",
        json: r#"{ "face": "z-max", "mpa": 1.0 }"#,
    },
    Entry {
        kind: "structure",
        key: "regions",
        json: r#"{ "material": "aluminium", "from": [0, 0, 0], "to": [4, 4, 4] }"#,
    },
];

/// The starting entry for `kind`'s list `key`, if the format has one.
pub fn entry(kind: &str, key: &str) -> Option<&'static Entry> {
    ENTRIES.iter().find(|e| e.kind == kind && e.key == key)
}

/// How an entry says where on its domain's grid it is, which is what a click in a viewport can
/// set for it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Placing {
    /// One whole face, by its `face`: a held or a pressed face.
    Face,
    /// A face, and a box of cells on it by `from` and `to` in the face's two axes, absent for the
    /// whole face: a cooled face, whose `area_cm2` is the box's when there is one.
    FacePatch,
    /// A box of cells by `from` and `to`, three each: a heat source, a region of another material.
    Cells,
}

/// One list whose entries say where on a grid they are.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Site {
    /// The `kind` of the domain that holds the list.
    pub kind: &'static str,
    /// The list's key.
    pub key: &'static str,
    /// How its entries say where.
    pub placing: Placing,
}

/// **The lists a viewport can point at**, by how each says where.
///
/// Read off the format's types, and held to them by
/// `every_list_that_says_where_is_one_a_click_can_set`: an [`ENTRIES`] entry naming a `face` is a
/// face site, one naming three-cell `from` and `to` is a cells site, and a face site takes a patch
/// exactly when its type accepts `from` and `to` — `CoolingSpec` does, and a held or pressed face
/// is refused one by `deny_unknown_fields`. A contact is a plane between cells rather than a face
/// or a box, and is not here.
pub const SITES: [Site; 6] = [
    Site {
        kind: "block",
        key: "cooling",
        placing: Placing::FacePatch,
    },
    Site {
        kind: "block",
        key: "dissipation",
        placing: Placing::Cells,
    },
    Site {
        kind: "block",
        key: "regions",
        placing: Placing::Cells,
    },
    Site {
        kind: "structure",
        key: "held",
        placing: Placing::Face,
    },
    Site {
        kind: "structure",
        key: "pressed",
        placing: Placing::Face,
    },
    Site {
        kind: "structure",
        key: "regions",
        placing: Placing::Cells,
    },
];

/// How `kind`'s list `key` says where its entries are, if it says.
pub fn placing(kind: &str, key: &str) -> Option<Placing> {
    SITES
        .iter()
        .find(|s| s.kind == kind && s.key == key)
        .map(|s| s.placing)
}

/// What [`entry_for`] fits an entry to: the domain it is being added to, as far as an entry cares.
#[derive(Clone, Copy, Debug, Default)]
pub struct Fit<'a> {
    /// The domain's `cells` and `cell_mm`, when it states both.
    pub grid: Option<([usize; 3], f64)>,
    /// The faces the domain's lists already name — every list, not only the one being added to:
    /// a face a structure holds is not a face to press, and pressing a clamped face moves nothing.
    pub faces: &'a [&'a str],
    /// The domain's own `material`, when it states one.
    pub material: Option<&'a str>,
}

/// The starting entry for `kind`'s list `key`, fitted to the domain it is added to.
///
/// [`ENTRIES`] is right for its kind's template and for no domain in particular, and for a block's
/// cooled face that is not a refusal but a wrong answer. `area_cm2` is the face's **whole** area,
/// which the domain shares out among its boundary cells, so the template's `1.21` added to scene
/// 21's wax, 11 x 11 cells of 2 mm, would state its 4.84 cm² top as a quarter of that — and the
/// build cannot know, because a stated area is the scene's to state. A box of cells beyond the
/// grid is refused by name, which is better, and still a button that only ever produces an error.
///
/// So what an entry says about where in a grid — a face's area, a box of cells, the plane a
/// contact sits on — is taken from [`Fit::grid`]: the whole of the first face no list holds yet,
/// the middle third of the grid, the lower half of it, the plane halfway up. A region is made of
/// something the domain is not, since a region of the block's own material changes nothing and
/// says nothing about having changed nothing. Everything else is [`ENTRIES`]' as written, and for
/// the template the two agree exactly, which is what `an_entry_fits_the_grid_it_is_added_to` pins
/// so they cannot drift apart.
///
/// `None` when the format keeps no entry for that list.
pub fn entry_for(kind: &str, key: &str, fit: &Fit) -> Option<String> {
    let start = entry(kind, key)?;
    // **A face nothing holds yet.** A block refuses a face cooled twice — "a face cannot lose heat
    // to two different airs" — so a second cooled face added to the top would be refused, and a
    // domain that has the top already wants the next face. The order is the top, the bottom, then
    // the sides; with all six taken it is the top, and the build says why.
    let face = FACES
        .iter()
        .find(|f| !fit.faces.contains(f))
        .copied()
        .unwrap_or("z-max");
    let substance = if fit.material == Some("aluminium") {
        "copper"
    } else {
        "aluminium"
    };
    match (kind, key) {
        ("structure", "held") => {
            return Some(format!(r#"{{ "face": "{face}", "as": "clamp" }}"#));
        }
        ("structure", "pressed") => {
            return Some(format!(r#"{{ "face": "{face}", "mpa": 1.0 }}"#));
        }
        _ => {}
    }
    let Some(([a, b, c], mm)) = fit.grid else {
        return Some(start.json.to_string());
    };
    let n = [a, b, c];
    let lo = n.map(|k| k / 3);
    let hi = [0, 1, 2].map(|i| (2 * n[i] / 3).max(lo[i] + 1));
    let half = n.map(|k| (k / 2).max(1));
    // The cells a face spans: the two counts that lie in it.
    let across = match &face[..1] {
        "x" => b * c,
        "y" => a * c,
        _ => a * b,
    };
    Some(match (kind, key) {
        ("block", "cooling") => format!(
            r#"{{ "face": "{face}", "ambient_c": 20.0, "convection_w_per_m2_k": 10.0, "area_cm2": {} }}"#,
            decimal(across as f64 * mm * mm / 100.0)
        ),
        ("block", "dissipation") => format!(
            r#"{{ "watts": 1.0, "from": [{}, {}, {}], "to": [{}, {}, {}] }}"#,
            lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]
        ),
        ("block" | "structure", "regions") => format!(
            r#"{{ "material": "{substance}", "from": [0, 0, 0], "to": [{}, {}, {}] }}"#,
            half[0], half[1], half[2]
        ),
        ("block", "contact") => format!(
            r#"{{ "axis": "z", "at": {}, "w_per_m2_k": 5000.0 }}"#,
            half[2]
        ),
        _ => start.json.to_string(),
    })
}

/// The six faces, in the order a new entry takes the first free one of: the top, the bottom, then
/// the sides. Spelled as the format spells them.
const FACES: [&str; 6] = ["z-max", "z-min", "x-max", "x-min", "y-max", "y-min"];

/// A float that stays a float in the file: `4.84`, and `4.0` rather than `4`.
fn decimal(v: f64) -> String {
    let s = format!("{v}");
    if s.contains(['.', 'e', 'E']) {
        s
    } else {
        format!("{s}.0")
    }
}
