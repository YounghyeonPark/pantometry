//! The thirty shipped scenes, offered as starting points.
//!
//! **A chooser asked "what are you simulating?" and answered with `bar`, `block` and `hall`.**
//! That is the scene format's vocabulary, which is the right vocabulary for a file and the wrong
//! one for somebody deciding what to make: the reader has to already know what a `puck` is to
//! discover that it is an espresso basket.
//!
//! So the first screen offers **what these are simulations of**, and the things it offers are the
//! scenes this repository already ships. Every one of them runs in CI on every commit and checks
//! itself against a closed form, so a starting point that works is guaranteed by tests that exist
//! rather than by a new set of examples somebody would have to keep true.
//!
//! # The areas are an assignment, and the set of scenes is not
//!
//! Which area a scene belongs in is a judgement — `24-a-power-module` is heat in a solid and is
//! filed under power electronics because that is what somebody opening it is doing. What is *not*
//! a judgement is the set: `every_shipped_scene_is_offered` holds these against the directory in
//! both directions, so a new scene has to be placed and a removed one has to go.
//!
//! # Embedded, because a binary does not know where the repository is
//!
//! `include_str!`, for the reason `templates` gives: `editor-core` compiles to `wasm32`, where
//! there is no disk to read a scene off. Twenty-five kilobytes for all thirty.

/// One scene, offered as a starting point.
pub struct Preset {
    /// The file it came from, which is also how it is held against the directory.
    pub file: &'static str,
    /// Which [`AREAS`] entry it is filed under.
    pub area: &'static str,
    /// The scene's own title, read from the file rather than written again here.
    pub title: &'static str,
    /// The kinds of domain it holds, so a reader can see what physics they are getting.
    pub kinds: &'static [&'static str],
    /// Whether it names a file beside itself — an STL for a designed part.
    ///
    /// **One does.** A preset opened as an unsaved project has nowhere to resolve that from, so
    /// the chooser says so rather than opening a scene that will refuse on its first check.
    pub needs_a_part: bool,
    /// The scene, as text.
    pub json: &'static str,
    /// A picture of its last frame, as PNG, or `None` when there is nothing to draw.
    ///
    /// **Three have none.** A `network` and two `winding`s report readings rather than places, so
    /// a run of them carries no panel and there is nothing for the shell to render. They are
    /// offered with their words and no tile, which is what they are.
    pub thumb: Option<&'static [u8]>,
}

/// What the areas are called, in the order they are offered.
///
/// `(key, name, what it is)`. The key is what a [`Preset`] files itself under.
pub const AREAS: [(&str, &str, &str); AREA_COUNT] = [
    (
        "heat",
        "Heat in solids",
        "conduction, melting, a coating, a part radiating across a gap",
    ),
    (
        "power",
        "Power electronics",
        "a busbar, a power module, and what cold does to its solder",
    ),
    (
        "motors",
        "Motors",
        "a winding, the steel around it, and the housing outside that",
    ),
    (
        "rooms",
        "Rooms and sound",
        "standing modes in two dimensions and in three",
    ),
    (
        "optics",
        "Light on a surface",
        "a beam, a lamp with a real spectrum, and a coating",
    ),
    (
        "flow",
        "Flow",
        "pressure-driven water, and liquid through a packed bed",
    ),
    (
        "matter",
        "Matter",
        "atoms in a lattice and out of one, and a particle in a well",
    ),
    (
        "radio",
        "Resonant cavities",
        "Maxwell's equations on a Yee grid",
    ),
    ("orbits", "Orbits", "bodies under their own gravity"),
    (
        "contact",
        "Contact",
        "a ball on a floor, losing energy to its dashpot",
    ),
    (
        "everything",
        "All of it at once",
        "five domains, four crates, one clock and one audit",
    ),
];

/// How many areas there are.
pub const AREA_COUNT: usize = 11;

/// Every shipped scene, as a starting point.
pub const PRESETS: [Preset; 30] = [
    Preset {
        file: "01-room-mode.json",
        area: "rooms",
        title: "a small room ringing in its (1,1) mode",
        kinds: &["room"],
        needs_a_part: false,
        json: include_str!("../scenes/01-room-mode.json"),
        thumb: Some(include_bytes!("../thumbnails/01-room-mode.png")),
    },
    Preset {
        file: "02-room-higher-mode.json",
        area: "rooms",
        title: "the (3,2) mode: more nodal lines, and a higher note",
        kinds: &["room"],
        needs_a_part: false,
        json: include_str!("../scenes/02-room-higher-mode.json"),
        thumb: Some(include_bytes!("../thumbnails/02-room-higher-mode.png")),
    },
    Preset {
        file: "03-room-pulse.json",
        area: "rooms",
        title: "a clap in the corner, spreading and coming back",
        kinds: &["room"],
        needs_a_part: false,
        json: include_str!("../scenes/03-room-pulse.json"),
        thumb: Some(include_bytes!("../thumbnails/03-room-pulse.png")),
    },
    Preset {
        file: "04-heater-and-bar.json",
        area: "heat",
        title: "a heater pays joules onto the bus and a bar takes them",
        kinds: &["bar", "heater"],
        needs_a_part: false,
        json: include_str!("../scenes/04-heater-and-bar.json"),
        thumb: Some(include_bytes!("../thumbnails/04-heater-and-bar.png")),
    },
    Preset {
        file: "05-beam-on-bar.json",
        area: "optics",
        title: "a beam that heats where it lands, and the heat spreading afterwards",
        kinds: &["bar", "beam"],
        needs_a_part: false,
        json: include_str!("../scenes/05-beam-on-bar.json"),
        thumb: Some(include_bytes!("../thumbnails/05-beam-on-bar.png")),
    },
    Preset {
        file: "06-orbits.json",
        area: "orbits",
        title: "four satellites, tilted out of one plane, and Kepler's third law",
        kinds: &["orbit"],
        needs_a_part: false,
        json: include_str!("../scenes/06-orbits.json"),
        thumb: Some(include_bytes!("../thumbnails/06-orbits.png")),
    },
    Preset {
        file: "07-bouncing-ball.json",
        area: "contact",
        title: "a ball bouncing on a penalty contact, and the heat its dashpot makes",
        kinds: &["bounce", "lump"],
        needs_a_part: false,
        json: include_str!("../scenes/07-bouncing-ball.json"),
        thumb: Some(include_bytes!("../thumbnails/07-bouncing-ball.png")),
    },
    Preset {
        file: "08-atoms-crystal.json",
        area: "matter",
        title: "a Lennard-Jones crystal at T* = 0.15: atoms rattling in place",
        kinds: &["atoms"],
        needs_a_part: false,
        json: include_str!("../scenes/08-atoms-crystal.json"),
        thumb: Some(include_bytes!("../thumbnails/08-atoms-crystal.png")),
    },
    Preset {
        file: "09-atoms-liquid.json",
        area: "matter",
        title: "the same atoms at T* = 1.4: the lattice is gone and they wander",
        kinds: &["atoms"],
        needs_a_part: false,
        json: include_str!("../scenes/09-atoms-liquid.json"),
        thumb: Some(include_bytes!("../thumbnails/09-atoms-liquid.png")),
    },
    Preset {
        file: "10-lamp-on-a-mirror.json",
        area: "optics",
        title: "a tungsten lamp on an aluminium mirror, and the heat the blue end leaves",
        kinds: &["bar", "light"],
        needs_a_part: false,
        json: include_str!("../scenes/10-lamp-on-a-mirror.json"),
        thumb: Some(include_bytes!("../thumbnails/10-lamp-on-a-mirror.png")),
    },
    Preset {
        file: "11-motor-thermal-network.json",
        area: "motors",
        title: "a winding, the steel around it and the housing outside that",
        kinds: &["heater", "network"],
        needs_a_part: false,
        json: include_str!("../scenes/11-motor-thermal-network.json"),
        thumb: None,
    },
    Preset {
        file: "12-winding-heats-a-motor.json",
        area: "motors",
        title: "the same motor, with the heat computed instead of stated",
        kinds: &["network", "winding"],
        needs_a_part: false,
        json: include_str!("../scenes/12-winding-heats-a-motor.json"),
        thumb: None,
    },
    Preset {
        file: "13-winding-that-heats-itself.json",
        area: "motors",
        title: "a winding whose resistance follows its own temperature",
        kinds: &["network", "winding"],
        needs_a_part: false,
        json: include_str!("../scenes/13-winding-that-heats-itself.json"),
        thumb: None,
    },
    Preset {
        file: "14-a-world.json",
        area: "everything",
        title: "a world: five domains, four crates, one clock and one audit",
        kinds: &["bar", "beam", "light", "orbit", "room"],
        needs_a_part: false,
        json: include_str!("../scenes/14-a-world.json"),
        thumb: Some(include_bytes!("../thumbnails/14-a-world.png")),
    },
    Preset {
        file: "15-a-hot-spot-in-a-block.json",
        area: "heat",
        title: "a hot spot in a block of aluminium, spreading in three dimensions",
        kinds: &["block"],
        needs_a_part: false,
        json: include_str!("../scenes/15-a-hot-spot-in-a-block.json"),
        thumb: Some(include_bytes!("../thumbnails/15-a-hot-spot-in-a-block.png")),
    },
    Preset {
        file: "16-a-room-with-a-ceiling.json",
        area: "rooms",
        title: "the oblique (1,1,1) mode of a room with a ceiling",
        kinds: &["hall"],
        needs_a_part: false,
        json: include_str!("../scenes/16-a-room-with-a-ceiling.json"),
        thumb: Some(include_bytes!("../thumbnails/16-a-room-with-a-ceiling.png")),
    },
    Preset {
        file: "17-a-busbar-with-a-notch.json",
        area: "power",
        title: "a busbar with a notch, and the resistance the shape actually has",
        kinds: &["conductor", "lump"],
        needs_a_part: false,
        json: include_str!("../scenes/17-a-busbar-with-a-notch.json"),
        thumb: Some(include_bytes!("../thumbnails/17-a-busbar-with-a-notch.png")),
    },
    Preset {
        file: "18-an-espresso-shot.json",
        area: "flow",
        title: "an espresso shot, and the same basket with a gap at its wall",
        kinds: &["puck"],
        needs_a_part: false,
        json: include_str!("../scenes/18-an-espresso-shot.json"),
        thumb: Some(include_bytes!("../thumbnails/18-an-espresso-shot.png")),
    },
    Preset {
        file: "19-a-coating-stops-the-heat.json",
        area: "optics",
        title: "a hot spot in aluminium, meeting a wall of borosilicate halfway",
        kinds: &["block"],
        needs_a_part: false,
        json: include_str!("../scenes/19-a-coating-stops-the-heat.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/19-a-coating-stops-the-heat.png"
        )),
    },
    Preset {
        file: "20-melting-a-block-of-ice.json",
        area: "heat",
        title: "a heater melting a block of ice, and the plateau it holds at while it does",
        kinds: &["block", "heater"],
        needs_a_part: false,
        json: include_str!("../scenes/20-melting-a-block-of-ice.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/20-melting-a-block-of-ice.png"
        )),
    },
    Preset {
        file: "21-a-wax-thermal-buffer.json",
        area: "heat",
        title: "a wax thermal buffer holding a plateau, in a substance the library does not ship",
        kinds: &["block", "heater"],
        needs_a_part: false,
        json: include_str!("../scenes/21-a-wax-thermal-buffer.json"),
        thumb: Some(include_bytes!("../thumbnails/21-a-wax-thermal-buffer.png")),
    },
    Preset {
        file: "22-wax-in-an-aluminium-matrix.json",
        area: "heat",
        title: "the same wax buffer with an aluminium matrix through it, declared as a composite",
        kinds: &["block", "heater"],
        needs_a_part: false,
        json: include_str!("../scenes/22-wax-in-an-aluminium-matrix.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/22-wax-in-an-aluminium-matrix.png"
        )),
    },
    Preset {
        file: "23-a-part-radiating-to-its-lid.json",
        area: "heat",
        title: "a hot part radiating across its clearance to a cooled lid",
        kinds: &["block"],
        needs_a_part: false,
        json: include_str!("../scenes/23-a-part-radiating-to-its-lid.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/23-a-part-radiating-to-its-lid.png"
        )),
    },
    Preset {
        file: "24-a-power-module-junction-to-ambient.json",
        area: "power",
        title: "a power module: 45 W from the die, junction to ambient through the stack",
        kinds: &["block"],
        needs_a_part: false,
        json: include_str!("../scenes/24-a-power-module-junction-to-ambient.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/24-a-power-module-junction-to-ambient.png"
        )),
    },
    Preset {
        file: "25-what-140-kelvin-does-to-the-solder.json",
        area: "power",
        title: "the same power module, and what 140 kelvin does to its solder",
        kinds: &["block", "structure"],
        needs_a_part: false,
        json: include_str!("../scenes/25-what-140-kelvin-does-to-the-solder.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/25-what-140-kelvin-does-to-the-solder.png"
        )),
    },
    Preset {
        file: "26-poiseuille-in-a-cooling-channel.json",
        area: "flow",
        title: "pressure-driven water in a 2 mm channel, against Poiseuille",
        kinds: &["channel"],
        needs_a_part: false,
        json: include_str!("../scenes/26-poiseuille-in-a-cooling-channel.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/26-poiseuille-in-a-cooling-channel.png"
        )),
    },
    Preset {
        file: "27-a-cavity-ringing-at-its-own-frequency.json",
        area: "radio",
        title: "a 120 x 120 mm vacuum cavity ringing in its (1,0,1) mode",
        kinds: &["cavity"],
        needs_a_part: false,
        json: include_str!("../scenes/27-a-cavity-ringing-at-its-own-frequency.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/27-a-cavity-ringing-at-its-own-frequency.png"
        )),
    },
    Preset {
        file: "28-an-eigenstate-that-does-not-move.json",
        area: "matter",
        title: "an electron in the third state of a 10 nm well, which does not move",
        kinds: &["well"],
        needs_a_part: false,
        json: include_str!("../scenes/28-an-eigenstate-that-does-not-move.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/28-an-eigenstate-that-does-not-move.png"
        )),
    },
    Preset {
        file: "29-a-designed-bracket-becomes-cells.json",
        area: "heat",
        title: "an aluminium bracket, designed as a mesh and solved as cells, cooling to still air",
        kinds: &["block"],
        needs_a_part: true,
        json: include_str!("../scenes/29-a-designed-bracket-becomes-cells.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/29-a-designed-bracket-becomes-cells.png"
        )),
    },
    Preset {
        file: "30-two-phases-crossing-at-a-clearance.json",
        area: "power",
        title:
            "two blackened busbars crossing at a clearance, one at twice the current, in still air",
        kinds: &["block"],
        needs_a_part: false,
        json: include_str!("../scenes/30-two-phases-crossing-at-a-clearance.json"),
        thumb: Some(include_bytes!(
            "../thumbnails/30-two-phases-crossing-at-a-clearance.png"
        )),
    },
];
