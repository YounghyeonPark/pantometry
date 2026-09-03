# `tools/presets` — the starting points the editor offers

`make.py` writes two things and both are **committed**:

| what | where | why it exists |
| --- | --- | --- |
| `presets.rs` | `app/pantometry-world/src/` | the thirty shipped scenes, grouped by what they are a simulation of, with each scene's text embedded |
| 27 PNG tiles | `app/pantometry-world/thumbnails/` | a picture of each scene's last frame, 240×156, 282 KiB in total and 10.4 KiB each |

```sh
cargo build --release --bin pantometry --manifest-path app/Cargo.toml
python tools/presets/make.py
```

It runs every scene, renders the last frame with `pantometry view --thumbnail`, sweeps away any
tile no scene names any more, writes the Rust and runs `rustfmt` over it. Running it against an
unchanged tree then reproduces `presets.rs` byte for byte and rewrites the tiles identically.

The `rustfmt` step is not tidiness. This script writes each `AREAS` entry on one line at up to 104
characters and rustfmt reflows them across four, so without it the committed file is **not** what
the script produces and the tree fails the app gate's first step. That sentence about reproducing
byte for byte was here before the step was, and it read as true only because `cargo fmt --all`
happened to run after the comparison that checked it.

## This is the exception to "nothing generated is committed"

That rule is stated twice — `EXAMPLES.md` and `app/pantometry-world/scenes/README.md` — and it is
about **what a run writes**: the SVG an example plots, the glTF a scene exports. Those are outputs.
Nobody has to have them for the repository to build, and a stale one is worse than none.

These are the other thing: **build inputs**. `presets.rs` is `include_str!`'d and the tiles are
`include_bytes!`'d, so they have to exist before `cargo build` starts. The alternatives were
considered and each fails on something concrete:

* **A `build.rs` that renders them.** It would need a GPU adapter at build time. CI's runners have
  none, `editor-core` compiles to `wasm32` where there is no adapter and no disk, and the whole
  point of the tiles is to be *in* the wasm binary.
* **Render them at startup.** The editor would need the scenes directory, which a shipped binary
  does not have — the reason `templates` and the scenes are embedded in the first place — and it
  would run thirty simulations to draw a menu.
* **Ship without pictures.** That is what the screen was, and it is the thing the user said was
  hard to choose from twice.

So the cost is 282 KiB in git, once per change to a scene, and the guard against the usual price
of a committed artefact — that it drifts from its source and nobody notices — is that the two
generated things are checked against the scenes rather than against themselves:

| check | where | what it would catch |
| --- | --- | --- |
| `every_shipped_scene_is_offered` | `pantometry-world` | a scene added or removed and not regenerated |
| `a_presets_title_and_kinds_are_the_scenes_own` | `pantometry-world` | a scene's title edited after the last run of this script |
| `a_preset_has_a_picture_unless_its_scene_has_nothing_to_draw` | `pantometry-world` | a tile that failed to render, which the screen would otherwise explain as "this one reports readings, not places" |
| `every_tile_decodes_to_something` | `pantometry-world` | a tile of the wrong size, or of nothing but background |
| `the_tiles_tell_the_scenes_apart_or_say_which_they_do_not` | `pantometry-world` | two scenes whose tiles are the same picture |
| `every_shipped_scene_puts_something_on_the_canvas` | `pantometry-app` | a scene that stopped producing a panel — the live half, run on every commit |
| `an_open_area_shows_its_scenes_and_a_shut_one_does_not` | `pantometry-app` | a tile that stopped decoding *in the app*, counted through the sentence the screen says about it |

**The third row is a pin on a generated file and not corroboration**, and it was described here as
"the other end of the same fact, from a GPU render rather than from the files", which was wrong
twice. Both it and the `pantometry-app` row trace to one `eprintln!` — `this run has no panels at
all`, which this script greps for and which the test greps for. And it is not a GPU render either
way: that branch runs before an adapter is requested, which is why the live half still reports the
same three on CI, where `ci.yml` says the runner has none.

## What is a judgement here and what is not

The `AREA` map and the `AREAS` table in `make.py` are editorial: which of eleven applications a
scene is filed under, and what that application is called on screen. A reader can disagree with
`24-a-power-module` being under power electronics rather than heat, and no test asserts the
assignment — one that did would only repeat it. What *is* asserted is that every area is named and
non-empty and that every scene is placed.

## Two tiles are nearly empty, and that is the renderer being honest

`07-bouncing-ball` lights 5 pixels of 37 440 and `06-orbits` lights 29, against 10 704 for a room
mode. Both are a handful of point bodies, and a point is drawn as a cross about 1.2% of the frame
across whatever the camera is fitted to — so cropping harder does not rescue them. Measured: at a
twentyfold magnification cap the bouncing ball is a small cross in the middle of a grey field,
which is a picture of the *marker* rather than of the run. `crop_to_content` caps at three for that
reason, and `every_tile_decodes_to_something` pins those two **by rank** — the two lowest, with
their counts — so a third one arriving is a failure rather than something to notice on the screen.
By rank because the threshold it replaced sat at 90 and the third-lowest tile is 94: four pixels of
37 440, which is not a margin.

## And four pairs of tiles are the same picture

27 tiles, **23 distinct images**. `08-atoms-crystal` and `09-atoms-liquid` draw the identical
picture; so do `24` and `25`, and so do `20`, `21` and `22`.

This is explicable rather than a fault. A tile's shade is the sample's value normalised over the
run's own range, and the frame is the last one — so three heat scenes in the same block, driven to
the same steady state by the same heater, produce the same *normalised* field whatever the material
between them is. `--frame` is not the culprit: `20-melting-a-block-of-ice` at frames 0, 6 and 11
gives three different pictures, and 21 and 22 reach the shared one by frame 6 while 20 does not
until 11.

Honest, and still bad on the screen: a chooser whose pictures cannot tell two scenes apart is a
chooser with no pictures, for those two. `the_tiles_tell_the_scenes_apart_or_say_which_they_do_not`
pins the collision set so a new one fails, and choosing a better frame — or drawing the run rather
than one instant of it — is the open question.
