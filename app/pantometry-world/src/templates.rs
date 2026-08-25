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
//! assumed — nineteen variants, nineteen distinct kinds across the twenty-eight scenes. So the
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
/// The examples, keyed by the `kind` the scene format spells.
///
/// Sorted by kind, so a menu built from this is in a stable order.
pub const TEMPLATES: [(&str, &str); 19] = [
    (
        "atoms",
        r#"{ "kind": "atoms", "name": "atoms", "cells": 3, "density": 0.8442, "temperature": 1.4,
          "thermostat_t": 1.4, "seed": 20260808 }"#,
    ),
    (
        "bar",
        r#"{ "kind": "bar", "name": "bar", "length_mm": 20.0, "cells": 61, "area_mm2": 100.0,
          "initial_c": 20.0 }"#,
    ),
    (
        "beam",
        r#"{ "kind": "beam", "name": "beam", "onto": "another domain in this scene", "faces": 61,
          "face_area_mm2": 0.5, "watts": 4.0, "reserve_j": 0.4, "waist_fraction": 0.08 }"#,
    ),
    (
        "block",
        r#"{ "kind": "block", "name": "block", "cells": [11, 11, 11], "cell_mm": 1.0,
          "initial_c": -10.0, "material": "ice" }"#,
    ),
    (
        "bounce",
        r#"{ "kind": "bounce", "name": "bounce", "drop_m": 0.5, "mass_kg": 0.0113,
          "stiffness": 200000.0, "damping": 20.0 }"#,
    ),
    (
        "cavity",
        r#"{ "kind": "cavity", "name": "cavity", "cells": [24, 4, 24], "cell_mm": 5.0,
          "medium": "vacuum", "mode": [1, 1], "amplitude_v_per_m": 1000.0 }"#,
    ),
    (
        "channel",
        r#"{ "kind": "channel", "name": "channel", "cells": [4, 16, 4], "cell_mm": 0.125,
          "fluid": "water", "walls": {"lower_m_per_s": 0.0, "upper_m_per_s": 0.0},
          "drive_m_per_s2": [0.02, 0.0, 0.0] }"#,
    ),
    (
        "conductor",
        r#"{ "kind": "conductor", "name": "conductor", "cells": [12, 5, 5], "cell_mm": 1.0,
          "resistivity_ohm_m": 1.724e-08, "volts": 0.001,
          "blocked": [[6, 0, 0], [6, 0, 1], [6, 0, 2], [6, 0, 3], [6, 0, 4], [6, 1, 0], [6, 1, 1], [6, 1, 2], [6, 1, 3], [6, 1, 4], [6, 2, 0], [6, 2, 1], [6, 2, 2], [6, 2, 3], [6, 2, 4]] }"#,
    ),
    (
        "hall",
        r#"{ "kind": "hall", "name": "hall", "width_m": 4.4, "height_m": 3.1, "depth_m": 2.4,
          "nodes_across": 23, "mode": [1, 1, 1], "amplitude_pa": 1.0 }"#,
    ),
    (
        "heater",
        r#"{ "kind": "heater", "name": "heater", "watts": 2.0, "reserve_j": 6.0 }"#,
    ),
    (
        "light",
        r#"{ "kind": "light", "name": "light", "watts": 40.0, "colour_k": 3200.0,
          "finish": "aluminium", "reserve_j": 2.0 }"#,
    ),
    (
        "lump",
        r#"{ "kind": "lump", "name": "lump", "volume_cm3": 2.0, "thickness_mm": 6.0,
          "initial_c": 20.0, "ambient_c": 20.0, "area_cm2": 12.0 }"#,
    ),
    (
        "network",
        r#"{ "kind": "network", "name": "network", "absorbing": "winding",
          "nodes": [{"name": "winding", "material": "copper", "volume_cm3": 18.0, "thickness_mm": 2.0, "initial_c": 25.0}, {"name": "stator", "material": "electrical_steel", "volume_cm3": 140.0, "thickness_mm": 8.0, "initial_c": 25.0}, {"name": "housing", "material": "aluminium", "volume_cm3": 220.0, "thickness_mm": 4.0, "initial_c": 25.0, "loses_to": {"ambient_c": 25.0, "area_cm2": 420.0}}],
          "links": [{"from": "winding", "to": "stator", "w_per_k": 0.9}, {"from": "stator", "to": "housing", "w_per_k": 2.4}] }"#,
    ),
    (
        "orbit",
        r#"{ "kind": "orbit", "name": "orbit", "central_kg": 1.989e+30,
          "radii_m": [57900000000.0, 108200000000.0, 149600000000.0],
          "inclinations_deg": [7.0, 3.4, 0.0], "satellite_kg": 3.3e+23 }"#,
    ),
    (
        "puck",
        r#"{ "kind": "puck", "name": "puck", "cells": [19, 10, 19], "cell_mm": 2.0,
          "radius_mm": 15.0, "grind_um": 250.0, "porosity": 0.45, "bar": 9.0, "brew_c": 93.0 }"#,
    ),
    (
        "room",
        r#"{ "kind": "room", "name": "room", "width_m": 4.4, "height_m": 3.1, "cells_across": 81,
          "release": {"as": "mode", "nx": 1, "ny": 1, "amplitude_pa": 1.0} }"#,
    ),
    (
        "structure",
        r#"{ "kind": "structure", "name": "structure", "cells": [8, 8, 8], "cell_mm": 1.5,
          "material": "copper",
          "held": [ { "face": "x-min", "as": "roller" }, { "face": "y-min", "as": "roller" },
                    { "face": "z-min", "as": "roller" } ],
          "follows": "another domain in this scene", "reference_c": 20.0 }"#,
    ),
    (
        "well",
        r#"{ "kind": "well", "name": "well", "cells": 200, "width_nm": 10.0,
          "start": {"eigenstate": 3} }"#,
    ),
    (
        "winding",
        r#"{ "kind": "winding", "name": "winding", "length_m": 62.0, "cross_section_mm2": 0.35,
          "amps": 1.75, "at_c": 90.0, "reserve_j": 200000.0 }"#,
    ),
];
