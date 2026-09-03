"""Write `presets.rs` and the tiles it embeds.

Run from anywhere, with a release build of the binary:

    cargo build --release --bin pantometry --manifest-path app/Cargo.toml
    python tools/presets/make.py

It runs every shipped scene, renders its last frame as a tile, and writes
`app/pantometry-world/src/presets.rs` embedding both the scene text and the tile. The README
beside this file says why that output is committed rather than built.

**Nothing here decides what a tile looks like.** The cropping, the magnification cap and the
shrink are `pantometry view --thumbnail`, in `app/pantometry-app/src/view.rs`, because they are
also what the editor's own snapshot does and a second copy of a rendering rule is a second thing
to get out of step. This script drives the binary and writes Rust. It had its own copy of the
cropping for an afternoon, and the two agreed only because one was transcribed from the other.
"""
import json
import pathlib
import shutil
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
APP = ROOT / "app"
EXE = APP / "target" / "release" / ("pantometry.exe" if sys.platform == "win32" else "pantometry")
SCENES = APP / "pantometry-world" / "scenes"
THUMBS = APP / "pantometry-world" / "thumbnails"
OUT = APP / "pantometry-world" / "src" / "presets.rs"

# Which application each scene is one of. Assigned by reading them; the *set* is held against the
# directory in both directions by a test, so a new scene has to be placed and a removed one has to
# go, but which area a scene went into is a judgement a reader can disagree with.
AREA = {
    "01": "rooms", "02": "rooms", "03": "rooms", "16": "rooms",
    "05": "optics", "10": "optics", "19": "optics",
    "11": "motors", "12": "motors", "13": "motors",
    "17": "power", "24": "power", "25": "power", "30": "power",
    "04": "heat", "15": "heat", "20": "heat", "21": "heat", "22": "heat", "23": "heat", "29": "heat",
    "18": "flow", "26": "flow",
    "08": "matter", "09": "matter", "28": "matter",
    "27": "radio",
    "06": "orbits",
    "07": "contact",
    "14": "everything",
}

# What each area is called on the screen, and one line about it. In the order they are offered.
AREAS = [
    ("heat", "Heat in solids", "conduction, melting, a coating, a part radiating across a gap"),
    ("power", "Power electronics", "a busbar, a power module, and what cold does to its solder"),
    ("motors", "Motors", "a winding, the steel around it, and the housing outside that"),
    ("rooms", "Rooms and sound", "standing modes in two dimensions and in three"),
    ("optics", "Light on a surface", "a beam, a lamp with a real spectrum, and a coating"),
    ("flow", "Flow", "pressure-driven water, and liquid through a packed bed"),
    ("matter", "Matter", "atoms in a lattice and out of one, and a particle in a well"),
    ("radio", "Resonant cavities", "Maxwell's equations on a Yee grid"),
    ("orbits", "Orbits", "bodies under their own gravity"),
    ("contact", "Contact", "a ball on a floor, losing energy to its dashpot"),
    ("everything", "All of it at once", "five domains, four crates, one clock and one audit"),
]

HEADER = r"""//! The thirty shipped scenes, offered as starting points.
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
pub const AREAS: [(&str, &str, &str); AREA_COUNT] = ["""


def tiles(work):
    """Render one tile per scene. Returns the stems that have one, and the stems that do not.

    A scene with no panel at all gets none, and the chooser says so in as many words. **That is
    the only reason a tile may be missing**: any other refusal raises here rather than quietly
    producing a preset that claims to report readings rather than places.
    `a_preset_has_a_picture_unless_its_scene_has_nothing_to_draw` holds the other end of it.
    """
    made, without = [], []
    for scene in sorted(SCENES.glob("*.json")):
        stem = scene.stem
        run = work / (stem + ".json")
        done = subprocess.run(
            [str(EXE), "run", str(scene), str(run)],
            cwd=APP, capture_output=True, text=True,
        )
        if done.returncode != 0:
            raise SystemExit(stem + ": the run refused:\n" + done.stdout + done.stderr)
        frames = len(json.loads(run.read_text(encoding="utf-8")).get("frames", []))
        if frames == 0:
            raise SystemExit(stem + ": the run wrote no frames")
        tile = THUMBS / (stem + ".png")
        shot = subprocess.run(
            [
                str(EXE), "view", str(run),
                "--snapshot", str(tile),
                "--frame", str(frames - 1),
                "--thumbnail",
            ],
            cwd=APP, capture_output=True, text=True,
        )
        said = shot.stdout + shot.stderr
        if "no panels at all" in said:
            without.append(stem)
            continue
        if shot.returncode != 0 or not tile.is_file():
            raise SystemExit(stem + ": the tile was not written:\n" + said)
        made.append(stem)

    # A scene that was removed leaves its tile behind, and a tile no scene names is ten kilobytes
    # of repository nothing reads. Swept, so the directory is the set rather than a superset.
    keep = set(s + ".png" for s in made)
    for stale in sorted(THUMBS.glob("*.png")):
        if stale.name not in keep:
            stale.unlink()
            print("  removed " + stale.name + ", which no scene names any more")
    return made, without


def rust(made):
    """Write `presets.rs`."""
    rows = []
    for f in sorted(SCENES.glob("*.json")):
        area = AREA.get(f.name[:2])
        if area is None:
            raise SystemExit(f.name + " has no area - add it to AREA in this file")
        s = json.loads(f.read_text(encoding="utf-8"))
        kinds = sorted(set(d.get("kind", "?") for d in s.get("domains", [])))
        needs = any("stl" in json.dumps(d) for d in s.get("domains", []))
        rows.append((f.name, area, s.get("title", ""), kinds, needs, f.stem in made))

    seen = set(r[1] for r in rows)
    named = set(a for a, _, _ in AREAS)
    if seen != named:
        raise SystemExit("areas assigned %s but named %s" % (sorted(seen), sorted(named)))

    out = [HEADER]
    for key, name, about in AREAS:
        out.append('    ("%s", "%s", "%s"),' % (key, name, about))
    out.append("];\n")
    out.append("/// How many areas there are.\npub const AREA_COUNT: usize = %d;\n" % len(AREAS))
    out.append(
        "/// Every shipped scene, as a starting point.\npub const PRESETS: [Preset; %d] = ["
        % len(rows)
    )
    for name, area, title, kinds, needs, has_tile in rows:
        ks = ", ".join('"%s"' % k for k in kinds)
        esc = title.replace("\\", "\\\\").replace('"', '\\"')
        stem = name[:-5]
        tile = ('Some(include_bytes!("../thumbnails/%s.png"))' % stem) if has_tile else "None"
        out += [
            "    Preset {",
            '        file: "%s",' % name,
            '        area: "%s",' % area,
            '        title: "%s",' % esc,
            "        kinds: &[%s]," % ks,
            "        needs_a_part: %s," % ("true" if needs else "false"),
            '        json: include_str!("../scenes/%s"),' % name,
            "        thumb: %s," % tile,
            "    },",
        ]
    out.append("];")
    OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
    # **rustfmt, because the gate's first step is `cargo fmt --all --check`.** This writes each
    # `AREAS` entry on one line at up to 104 characters and rustfmt reflows them across four, so
    # the committed file was *not* what this script produces. The README said "byte for byte",
    # which read as true only because `cargo fmt --all` happened to run after the comparison that
    # checked it. Formatting here makes the sentence true rather than deleting it.
    done = subprocess.run(
        ["rustfmt", "--edition", "2021", str(OUT)], capture_output=True, text=True
    )
    if done.returncode != 0:
        raise SystemExit("rustfmt refused:\n" + done.stdout + done.stderr)
    return rows


def main():
    if not EXE.is_file():
        raise SystemExit(
            "no binary at %s\n"
            "  cargo build --release --bin pantometry --manifest-path app/Cargo.toml" % EXE
        )
    THUMBS.mkdir(exist_ok=True)
    work = pathlib.Path(tempfile.mkdtemp(prefix="pantometry-presets-"))
    try:
        made, without = tiles(work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    rows = rust(made)
    size = sum(p.stat().st_size for p in THUMBS.glob("*.png"))
    print("%d presets across %d areas" % (len(rows), len(AREAS)))
    print(
        "%d tiles, %.0f KiB in total, %.1f KiB each"
        % (len(made), size / 1024, size / max(len(made), 1) / 1024)
    )
    print("%d scenes have no panel and are offered without one: %s" % (len(without), without))
    counts = {}
    for r in rows:
        counts[r[1]] = counts.get(r[1], 0) + 1
    for key, name, _ in AREAS:
        print("  %2d  %s" % (counts.get(key, 0), name))


main()
