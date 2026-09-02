# A scene editor beside a 3D view

```sh
cd app
cargo run --release                          # opens on the start screen
cargo run --release -- scene.json            # opens on a file
cargo run --release -- scene.json --run      # and runs it at once
```

With no file named you arrive at **New project**, **Open a scene** and whatever you had open
before — not at the built-in room, which used to be the answer whatever you meant by starting the
editor. New project asks what you are simulating: every kind the scene format defines, with a line
about each. Tick one for a starting point or several for a scene of several domains, which is what
a custom simulation is here; there is no second flow for that case, because it is the same list.

It says two things while you choose that nothing downstream can. **These kinds do not share a
timescale** — `atoms` settles in 6e-12 s and a thermal `network` in 1800 s, and a scene carries one
`duration_s` — so a set that spans orders of magnitude is a *well-formed* scene in which the slowest
domain barely moves, and the chooser puts the ratio on the screen. And the two kinds that name
another domain, `beam` and `structure`, say so rather than being wired up: a beam's `faces` has to
equal the bar's `cells` and its `onto` has to match the bar's `exposes`, so what opens instead is
the format's own complaint, by name.

The duration and frame count a new scene starts from are not chosen: each comes from a shipped
scene that uses that kind, and a test holds every one of them
against the thirty scenes on disk.

The viewport is an instrument, not a preview: a colour bar with numbers on it, a scale bar in
model metres, a probe that names whatever the cursor is over, and a transport — play, step either
way, and the space bar, the arrow keys, Home and End. What it draws is still chosen by the shape
of the data and never by a domain's name.

The left pane is the scene's JSON, checked **as you type** with the same two steps
`pantometry-world --check` runs — parse, then build — with parse errors carried as `line:column`,
which is what that error format was designed for.

The right pane is the inspector, and it **writes**. Selecting a row that names a domain — its
extent before a run, its field or its readings after one — lists every value the scene states
about it, at whatever depth it sits: a number gets a drag, a `material` gets a menu of the
catalogue plus whatever the scene declared, any other string gets a field. A change goes back into
the text and re-checks on the same frame. `kind` and `name` are the two it holds back, both
because they identify rather than describe — the format refuses unknown fields when a `kind`
changes, and five other keys refer to a domain by `name`.

**Domain → Add** puts a new domain in, from a starting example of each of the nineteen kinds the
format defines, named so it does not collide with anything already there; **Domain → Delete**
takes the selected one out along with the comma that held it. The list lives in
`pantometry_world::templates`, beside the format rather than beside the editor, and a test
compares it against the kinds the shipped scenes use **in both directions** — so a twentieth
domain appears in this menu without the shell learning about it, and cannot quietly not appear.
Two kinds arrive incomplete because the format will not let them stand alone: a `beam` states
what it shines `onto` and a `structure` states the block it `follows`, and neither can be guessed.

Below those, **placement**: where the domain's origin sits in the world, and the turn — an axis
and degrees — applied about it. It is a separate control
because `poses` is a map beside `materials` rather than a field inside the domain, and it is the
one write that **creates** instead of replacing — no shipped scene states `poses` at all, so
moving anything means writing a key the file does not have. Three levels can be missing and each
is filled in where it is found; a scene that gains a pose keeps every other byte.

Every one of those is a byte splice rather than a re-serialise, so the file keeps its own
formatting and key order. `editor_core::set_number`, `set_text`, `set_pose`, `add_domain` and
`remove_domain` are where that lives and where the tests for it are — including one that adds each
of the nineteen kinds and removes it again, and asserts the file is byte-identical to what it
started as.

Nothing in any of it enumerates domain kinds, so a domain written out of tree is as editable as a
shipped one. The read-only rows below say what came *out*, and the two halves are kept apart on
purpose: a peak temperature is not a thing to drag.

A selected domain also gets **three translate handles** in the viewport — x red, y green, z blue.
Dragging one moves the domain along that axis and writes the pose; the camera does not turn while
a handle is held, because to egui the two gestures are the same gesture. A handle pointing at the
camera refuses to move rather than dividing by nearly nothing, and letting go over empty space
does not clear the selection — a drag that ends on nothing reads as a click.

The arithmetic is `editor_core::drag_along_axis`, and its error has an **order** rather than a
tolerance: halving a drag quarters the miss. That is how the first version was caught moving
things 77% of the way at every drag size, which no "close enough" assertion would have noticed.

The viewport draws every placed extent as a wireframe, live from the text before anything runs.
**Run** streams the run in as it computes
— each frame appears when it is captured, the slider grows, **stop** ends a long run between
frames — and **Verify** runs the battery from `pantometry-world verify` and shows the report the
CLI prints, with the findings count in the window title. Drag to rotate, scroll to zoom.

## The viewport is shaded, on the GPU, with a depth buffer

**View → Shaded surfaces** is the default. A field is drawn as the *boundary of the cells that
hold a value* — one quad per cell face whose neighbour is absent — lit by a key and a fill, depth
tested against everything else in the scene, with the placed extents as wire boxes that are hidden
where they are behind a solid. A body is a sphere at the radius the exports use. Paths are lines.

**A scene's `parts` are drawn as the meshes they are**, flat-shaded, in the place their voxels
will occupy — which is the only thing in this viewport that is a picture of the *scene* rather
than of a run, and so the only thing visible before anything has been run. It goes away for a
domain once that domain has a panel in the current frame: the field's surface is in the same
space, and two surfaces in one place z-fight into a picture of neither.

The pose is `Scene::placements`', the same one `Voxels::onto` rasterises against. That is the
whole of the correctness here — a part drawn a millimetre from its own voxels is two pictures of
one object, and the viewport offers nothing that would say which of them moved. Pinned by
`a_part_is_drawn_where_it_is_rasterised.rs`, whose quarter turn is worked out by hand rather than
by the code under test.

**View → Cells (see inside)** is the other picture, and the menu says so: every sample as a
translucent splat composited far to near. Neither is a rendering of the other. A surface says what
the shape is and a splat cloud says what is *inside* it, and the block whose hot spot is in the
middle looks uniform from outside because it **is** uniform from outside.

Three things this is careful about, each because of a defect this workspace has already shipped:

| | |
| --- | --- |
| the geometry is `pantometry_view::mesh` | the same function that writes glTF and USD. A 40 mm cube once exported 80 mm across — a field is sampled corner to corner, so an end node owns half a cell — and a second copy of that arithmetic here would be a second chance at it |
| the colour is `editor_core::Colouring` | Planck's law where the field is hot enough to glow, `pantometry_view::ramp` where it is not, converted to **linear** for the shader. The glTF exporter shipped sRGB into a linear slot and every export was about 2.3× too bright in the midtones |
| vertices are framing-local | `Framing::local` subtracts the centre in `f64` and only the result narrows. Folding the centre into the matrix instead was measured disagreeing with `Camera::project` by `4.4e-6` on a 9 mm block sitting 200 mm out — `f32` keeping two digits of seven out of a subtraction |

`Camera::matrix` and `Camera::project` are **one camera**: the flat painter walks points through the
second, the GPU multiplies every vertex by the first, and `one_camera_two_paths` holds them to
`2.4e-7` of each other over 23 000 points. They have to be, because the shaded pass is underneath
and every caption, the colour bar and the probe are egui on top of it. A camera with two ideas of
which way is up puts the labels where the geometry is not.

### Reading the viewport back, because a screenshot is not evidence

```sh
PANTOMETRY_VIEWPORT=1 PANTOMETRY_VIEWPORT_SHOT=out.ppm cargo run --release -- scene.json --run
```

`PANTOMETRY_VIEWPORT` prints the callback's rect and what each paint drew; `PANTOMETRY_VIEWPORT_SHOT`
reads the pixels back out of the framebuffer with `glReadPixels` and writes them as a PPM, with a
count of how many are not the background.

That exists because **capturing the window from outside does not work and does not fail loudly**.
`PrintWindow` returned a convincing image of the panels with the 3D content simply absent, and a
plain screen grab returned a *different application's* window, because a background process cannot
raise a window on Windows. Both looked exactly like a renderer that drew nothing. What settled it
was the pass reporting `1620 triangles, 12 lines` and then handing over its own pixels.

Two real defects were behind that hunt, and neither was the renderer:

- The default window is 533 points wide. Three side panels have minimum widths, and a `SidePanel`
  is served before the `CentralPanel` is — so the viewport was handed a rect of **zero width**. It
  drew nothing, correctly. The window opens at 1500 × 950 now and the panels have a ceiling.
- One outliner row is a scene's title. On a 1706-pixel window that row took the outliner to
  **940 px** and left the viewport nothing: a 3D editor showing no 3D because of a caption. The rows
  truncate now, with the full name on hover.

## One home per command

The strip under the menus carries four things: **run**, **verify**, **deep** and **fit view**.
It carried nine, and seven of those nine were also in a menu — a toolbar and a menu bar offering
the same command is two answers to one question, and the screen did not say which was the real
one. What stayed is what you reach for without thinking; what went back to the menus is what you
look up. `deep` stayed because it changes what **verify** does when you press it, so it belongs
beside that button; **Watch file** and **Run on change** went because they change what happens
when nobody presses anything, and that is a thing to *read* — the status bar says so on the right
— rather than a control to keep permanent room for.

The file's name sits at the left of that strip, with the full path on hover, and says **edited**
when the pane is dirty. It used to be an editable 220-point path field, which is what it was for:
typing a path was the only way into a file on disk. `Open…` and `Save as…` are that now.

Seven menus, and each is named for one subject: **File** opens and writes, **Edit** takes edits
back and chooses a handle, **View** is about the picture, **Window** about which panels are on
screen, **Domain** adds and removes one, **Run** runs and verifies, **Watch** follows the file.
`View` used to hold eleven items doing three of those jobs at once — it was the only menu whose
name did not predict its contents, and the panels were the part that was not about the picture at
all. A seventh name to read is the price; a menu you can guess the contents of is what it buys.

`Solo the selection` came out of that menu without going anywhere: the outliner's header already
carried the checkbox, beside the selection it filters and with its state visible. Hiding the
outliner would then have left the viewport drawing one domain out of five with nothing saying why,
so the status bar says `solo — drawing only the selection` exactly when the control is not on
screen to say it itself.

The same line applies to menus: a menu is where a thing is turned on and where a command is issued,
and it is not where a continuous value is set. **Isosurface at a value** is a menu item; the level
it opens at is a slider on the strip along the bottom of the viewport, under the colour bar whose
scale it shares. It was inside the menu, which is drawn over the picture the level moves in.

## The live loop

Leave the window open and let a script do the editing. With **Watch file** on (the default —
it is in the Watch menu, and the status bar says when it is on) the editor polls the file's
modified time and reloads when something else writes it; with
**Run on change** on as well, the reload runs the scene — so
`script writes → editor rechecks → runs → draws` closes with no hand on the window, which is
the loop an agent-driven workflow needs. A change arriving mid-run stops the in-flight run at
the next frame boundary and starts a fresh one, so the picture converges on the latest text
rather than queueing history.

Two honesty rules hold it together. **Unsaved edits are never clobbered**: if the pane is
dirty and the disk changes, the status line says so and the disk's version waits for an
explicit **revert** — which of two writers meant it is not the editor's call. That button was
called `load` and re-read the path already in the box, which made it look like the way to open a
different file when it was the way to throw away your edits; `Open…` and `Save as…` are beside it
now, through the platform's own dialog. And **a prefix is
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
the library produced — and the library's gate could not see it, because these are
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

## It links `pantometry`, and `viewer-core` beside it does not

This used to be its own workspace — the third of three above the library — and the measured
argument for that boundary is in [app/README.md](../README.md), along with why there is one
boundary now instead of four.

**The editor links `pantometry` deliberately, and the difference from the viewer is the point.**
The viewer proves the wire format is complete by never linking the library; the editor exists to
build, run and verify scenes, which cannot be done from a file alone, and is the first consumer of
the platform verbs from a GUI. The two sit in one workspace sharing one `[workspace.dependencies]`
table now, so `viewer-core`'s independence is held by `the_wire_format_is_enough` rather than by a
directory.

Where the two overlap, they share: the camera, the framing, the fit and the projection are
`viewer-core`'s, imported rather than rewritten, because that arithmetic has been wrong here before
and twice is enough.

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

## An edit can be taken back

Six ways to change a scene arrived before any way to change one back, and two of them destroy
something: a removed domain takes its whole block of JSON, and a drag on an axis rewrites a
number nobody is necessarily looking at.

`History` is a list of texts and an index into it, which is cheap here because every edit in
`edit` already returns a **whole new string** rather than mutating one. `Edit > Undo`, `Ctrl+Z`,
redo on `Ctrl+Shift+Z`, sixty-four states deep with the oldest dropped.

**The buffer changes two ways and only one of them passes through this crate.** The shell binds
the text pane straight to the scene text, so typing never reaches a splice. A history that
recorded only the splices would hold a snapshot from *before* the typing, and undoing a deleted
domain would discard whatever had been typed since — silently, which is the failure mode this
repository keeps finding.

So the current text is committed on **both sides** of every edit. The first call folds any
typing into a state of its own and is a no-op when there was none; the second records the edit.
Undo then walks back through the typing rather than over it.

Committing before a **redo** as well. Typing leaves the buffer holding a state the history has
never seen, and stepping forward would overwrite it; folding in first makes the typing the
newest state, so redo says there is nothing ahead — which is true — instead of reaching a future
by discarding the present. The first draft committed only before undo, and redo after typing
lost it.

egui's own fine-grained undo is left alone for whatever has focus, so `Ctrl+Z` is handled only
when nothing does — the same guard the frame-stepping keys use, and for the same reason. Two
undos on one key would take back different things depending on where the caret happened to be.
The menu item works either way.

Opening or reverting from disk **resets**. A file is a different document, and undoing from one
scene back into another would offer a text belonging to no path.

## What the first version does not do

No structured editing — the text *is* the model, and `deny_unknown_fields` plus the version
check are what stand between a typo and a silently different scene. No shadows, no ambient
occlusion and no anti-aliasing beyond whatever the context gives: the lighting is a key, a fill and
an ambient, which is the least that shows a shape. No file dialogs. Undo covers the editor's own
edits and leaves typing to the text box, which is the division described above rather than a gap.
Each of those is worth adding only after using this one reports back what actually
chafes; that method has a name here, and `FRICTION.md` is where its findings go.
