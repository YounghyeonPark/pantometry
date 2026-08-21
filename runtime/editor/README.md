# A scene editor beside a 3D view

```sh
cd runtime/editor
cargo run --release                 # opens on the built-in room
cargo run --release -- scene.json   # opens on a file
```

The viewport is an instrument, not a preview: a colour bar with numbers on it, a scale bar in
model metres, a probe that names whatever the cursor is over, and a transport — play, step either
way, and the space bar, the arrow keys, Home and End. What it draws is still chosen by the shape
of the data and never by a domain's name.

The left pane is the scene's JSON, checked **as you type** with the same two steps
`pantometry-world --check` runs — parse, then build — with parse errors carried as `line:column`,
which is what that error format was designed for. The viewport draws every placed extent as a
wireframe, live from the text before anything runs. **Run** streams the run in as it computes
— each frame appears when it is captured, the slider grows, **stop** ends a long run between
frames — and **Verify** runs the battery from `pantometry-world verify` and shows the report the
CLI prints, with the findings count in the window title. Drag to rotate, scroll to zoom.

## The live loop

Leave the window open and let a script do the editing. With **watch file** on (the default)
the editor polls the file's modified time and reloads when something else writes it; with
**run on change** on as well, the reload runs the scene — so
`script writes → editor rechecks → runs → draws` closes with no hand on the window, which is
the loop an agent-driven workflow needs. A change arriving mid-run stops the in-flight run at
the next frame boundary and starts a fresh one, so the picture converges on the latest text
rather than queueing history.

Two honesty rules hold it together. **Unsaved edits are never clobbered**: if the pane is
dirty and the disk changes, the status line says so and the disk's version waits for an
explicit `load` — which of two writers meant it is not the editor's call. And **a prefix is
never dressed as the run**: a streaming or stopped run is labelled on the canvas, because a
partial run that looks complete is a picture of something that did not happen.

The streaming itself is `editor-core::run_streaming`, built on `World::advance` — one
iteration of `World::run`'s loop, made public when the editor became its second consumer —
and its tests pin that the final streamed payload is byte-identical to the batch run's.

## Three things it was doing itself, and one of them wrong

`editor-core` exists so that anything which could be got wrong the same way twice is written once.
Three things had escaped it, and the failure they produced is the argument for the rule.

**The paint order.** `Camera::project` returns a depth — *distance from the eye, larger is
further* — and this shell threw it away, then recovered a stand-in by projecting each point a
second time one millimetre along world z and taking the reciprocal of how far the two landed apart
on screen. That is proportional to how far **off the view axis** a point is, not to how far away
it is. On a 6×6×6 lattice of splats at the camera the app opens at, it ordered **6,026 of 23,220
pairs backwards — 26%** — so a quarter of every translucent volume composited in the wrong order,
worst at the centre of the screen, which is where the object is. The browser shell had always
sorted by the real depth and was right. Both call `editor_core::far_to_near` now.

**The colour.** Both shells carried a copy of a straight blue-to-red line in sRGB, which made four
spellings of "cool to warm" in this workspace. Measured over 256 steps it covered **16.6 L\*** —
44.1 to 60.7 — against the 74 the library's scale covers, and seventy-five of its 255 steps ran
backwards. Lightness is what survives a greyscale print, a projector and the eight per cent of men
with a colour-vision deficiency, and blue against red is the one pair those readers cannot
separate. `editor_core::value_colour` is `pantometry::view::ramp`, whose properties are pinned by
tests in the library, and it picks the diverging scale for a range that straddles zero with the
neutral at **the value zero** rather than at the middle of the range.

**The numbers.** Readings printed `{:.4}`, so a cavity holding 3.19e-10 J showed `0.0000` beside a
field the same run had reported at 921 V/m — the third writer in this workspace to erase a
magnitude that way, after `readings_csv` and `to_json`.

## Where a field's box comes from

The run says. `PanelData::Field` carries `extent_m` — the box it was sampled over — so a run opened
without the file that produced it still draws in the right place and at the right size. The
scene's placed box is the fallback for a run written before the format carried it.

That key is also why the viewer's workspace needed a change on the same day: its reader is
`deny_unknown_fields`, so the moment the library began writing `extent_m` it refused every run file
the library produced — and the library's twenty-step gate could not see it, because these are
separate workspaces with separate CI jobs and nothing in `crates/` reads that format back. Adding a
key to the wire format is a coordinated change across three workspaces, and the gate is not what
tells you so.

## The colour is computed, not chosen

A temperature field is drawn in the colour a body at that temperature **actually is** — Planck's
spectral exitance through the CIE 1931 colour matching functions to sRGB, from
`pantometry::view::colour`. A block at 1473 K is that orange and no other, and nothing here picked
it; in Blender the same block's colour is a value somebody types.

Which is also the honest limit. Below about 900 K a body emits no visible light at all —
`glow_fraction` says so with the visible share of its radiated power rather than a threshold
somebody chose — and this workspace holds no visible *reflectance* to draw instead, since
`Substance` carries a broadband infrared emissivity and no colour. So a cool field falls back to
the conventional ramp and the canvas says **"false colour — nothing here is hot enough to
glow"**. A false colour mistaken for a real one is a wrong answer that looks right, and the
label is what stops that.

The physics is checked against closed forms and published values, not against a second
renderer: Wien's displacement law fixes where the Planck curve peaks, the curve integrates to
`σT⁴`, equal-energy white lands on `x = y = 1/3`, and the locus passes through CIE Illuminant A
at `(0.44757, 0.40745)` — a coordinate this code had no hand in producing.

## The same editor, in a browser

`editor-wasm` compiles this machinery to `wasm32-unknown-unknown` and `site/` is a page over it:
edit, check as you type, run, scrub, verify — with **no server**, because the whole library
builds for the browser and a backend would only be a slower copy of what is already in the tab.
See [site/README.md](site/README.md), and `node editor-wasm/selftest.mjs` for the browser path
exercised without a browser.

One core, two shells. `editor-core` owns checking, placement geometry, running, verifying and
the colour rule; `editor` is a native window and `editor-wasm` is a page. The camera is
`viewer-core`'s in both. A thing that could be got wrong the same way twice is written once.

## Why this is a third workspace

The measured table, one row longer:

| | external crates |
| --- | --- |
| the library, all seventeen published crates | **12** |
| `bindings/python` | 15 |
| `runtime/viewer` | 86 |
| here | **371** |

`deny.toml`, `--locked`, `wasm32` and Rust 1.78 are the library's promises and none of them
can carry a window toolkit. Nothing in `crates/` depends on this; the arrow points one way.

**Unlike the viewer, this workspace links `pantometry`** — deliberately, and the difference is the
point. The viewer proves the wire format is complete by never linking the library. The editor
exists to build, run and verify scenes, which cannot be done from a file alone; it is the
first consumer of the platform verbs from a GUI. Where the two overlap, they share: the
camera, the framing, the fit and the projection are `viewer-core`'s, imported rather than
rewritten, because that arithmetic has been wrong here before and twice is enough.

## The two halves, kept apart by crate

`ARCHITECTURE.md`'s platform rules: the **authoring half** is the composition root and may
name domains; the **inspection half** must dispatch on the shape of the data, so a new physics
costs no editor edit.

- `editor-core` is the authoring half's machinery, GUI-free and tested headlessly: checking
  (parse + build), placement geometry (posed corners of every extent, through
  `Pose::point_to_world` even though no scene can state a pose yet), and thin passes over run
  and verify. It knows domains only through `pantometry-world`'s public API, where that knowledge
  already legitimately lives.
- `editor` is the shell: a text box, a canvas and an event loop. Everything it paints, it
  paints from a shape — a box, points, paths, a reading — with the domain name used only as a
  label.

## What the first version does not do

No structured editing — the text *is* the model, and `deny_unknown_fields` plus the version
check are what stand between a typo and a silently different scene. No field rendering in the
viewport — a field's box carries a note pointing at the HTML report, which is stated on the
canvas rather than quietly drawn as nothing. No file dialogs, no undo beyond the text box's
own. Each of those is worth adding only after using this one reports back what actually
chafes; that method has a name here, and `FRICTION.md` is where its findings go.
