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
timescale** — a quantum `well` settles in 2e-13 s and an `orbit` takes 7200, and a scene carries one
`duration_s` — so a set that spans orders of magnitude is a *well-formed* scene in which the slowest
domain barely moves, and the chooser puts the ratio on the screen. And the two kinds that name
another domain, `beam` and `structure`, say so rather than being wired up: a beam's `faces` has to
equal the bar's `cells` and its `onto` has to match the bar's `exposes`, so what opens instead is
the format's own complaint, by name.

The duration and frame count a new scene starts from are not chosen: each comes from a shipped
scene that uses that kind, and a test holds every one of them
against the thirty-four scenes on disk.

The viewport is an instrument, not a preview: a colour bar with numbers on it, a scale bar in
model metres, a probe that names whatever the cursor is over, and a transport — play, step either
way, and the space bar, the arrow keys, Home and End. What it draws is still chosen by the shape
of the data and never by a domain's name.

## A model tree, and a form for each row of it

The outliner is the scene's **model tree**. The scene; every domain it states, placed or not, with
its kind in the type column; under each domain the lists it holds — `cooling`, `regions`, `nodes`,
`perfusions` — with a row per entry; and the scene-level maps of things, `materials` or `poses`,
as groups of their own. After a run, the run's fields and readings follow. It used to list placed
*extents*, so a thermal network or a compartment model was not in the tree at all until a run had
reported a scalar about it — and a domain that is not in the tree cannot be selected to be changed.

It is read from the JSON, not from the parsed scene, so a scene the format refuses still has a tree
to select from while it is being fixed, and from the **shape** of the file rather than from a list
of kinds: a domain written out of tree gets the same rows. `editor_core::model` is all of it.

The pane beside the tree has two tabs. **Settings** is a form for the selected row — every value
the row states, in the file's order, with the unit its key names: a number gets a drag, an array of
numbers one line of drags, a string the format reads as one of a set a menu, any other string a
field. With nothing selected it is the scene's own form: its duration, its frames, its schedule.
**Source** is the scene's JSON, checked **as you type** with the same two steps
`pantometry-world --check` runs — parse, then build — with parse errors carried as `line:column`.

**The text is still the model.** The form reads it on every frame and writes it by splice, so the
two cannot disagree; every change can be undone, and the check runs on the same frame, its verdict
shown at the top of the form as well as on the canvas — an add can produce a scene the build
refuses, and the sentence that says so belongs beside the button that caused it.

**Every value has one home.** A domain's form holds what it states outside its lists, and each
entry's form holds the entry; together, over every row `tree` makes, they are every value in the
file but `kind` and `name`. `tests/the_model_tree_sets_up_a_scene.rs` checks that over every shipped
scene against a walk of the JSON written in the test, which is the property a form can lose
quietly: a value no row shows is part of the file nothing but the source can change. A name with
a slash in it — a material called `FR-4/Cu` — is one row, not two; split as a path, its form was
empty and its values had no home anywhere.

**A menu's words are the format's own.** The value is replaced with a word no enum has and the scene
parsed: serde's refusal names every variant it expected. Each is then tried, and only those the
scene still parses with are offered. So a cooled face offers the six faces, a hold `clamp` and
`roller`, the schedule its three, and a room's `"as": "mode"` no menu at all — `pulse` is a variant
of the same enum, and one that would make the keys beside it unknown. A `material` keeps the
catalogue plus whatever the scene declared. Nothing in the editor lists a set of words.

**An entry is added, copied and taken out.** A list's row has `+ <list>`, which adds the starting
entry the format keeps for it — `pantometry_world::templates::ENTRIES`, one for each of the
seventeen lists of objects the format has — and a domain's row offers the lists it does not hold
yet. An entry's row has `duplicate` and `delete`. A copy of a named node is renamed, `winding 2`,
for the reason a second domain of a kind is.

The starting entry is **fitted to the domain it lands in**, by `templates::entry_for`, and the one
reason is worth stating. A block's `area_cm2` is the whole cooled face's area, which the domain
shares out among its boundary cells, so the template's `1.21` cm² added to scene 21's wax — 11 ×
11 cells of 2 mm — would have stated its 4.84 cm² top as a quarter of itself, and the build cannot
know: a stated area is the scene's to state. So a cooled face takes the first face its list does not hold yet and
states that face's area, and a box of cells or a contact plane lands inside the grid. The test that
found the other half of this was adding a cooled face to every shipped block: scene 23's housing
already cools its top, and a second entry on one face is refused — "a face cannot lose heat to two
different airs". A new face avoids the faces *every* list names, so a structure is not pressed
where it is held, and a region is made of something the domain is not.

**When a menu is a text box.** A scene that fails to parse at an earlier byte has its refusal
about that byte, so the menus below it are text boxes until the file is fixed, with the error at
the top of the form. A named fluid or medium, `"water"` or `"vacuum"`, is free text, because the
format reads any string there and resolves the name at build.

`kind` and `name` are the two keys no form offers, both because they identify rather than describe
— the format refuses unknown fields when a `kind` changes, and five other keys refer to a domain by
`name`. Both stay jobs for the Source tab.

**The inspector says what came out.** It holds a row's box, its range over the run and this frame's
spread, and no inputs: those are on the form beside the tree. The two used to share the inspector,
inputs above results, and the line between them is now a line between two panels.

**Domain → Add** puts a new domain in, from a starting example of each of the twenty-two kinds the
format defines, named so it does not collide with anything already there; **Domain → Delete**
takes the selected one out along with the comma that held it. The list lives in
`pantometry_world::templates`, beside the format rather than beside the editor, and a test
compares it against the kinds the shipped scenes use **in both directions** — so a twentieth
domain appears in this menu without the shell learning about it, and cannot quietly not appear.
Two kinds arrive incomplete because the format will not let them stand alone: a `beam` states
what it shines `onto` and a `structure` states the block it `follows`, and neither can be guessed.

A placed domain's form ends with **placement**: where the domain's origin sits in the world, and
the turn — an axis and degrees — applied about it. Only a placed one: a thermal network or a
compartment model is not at anywhere, the build takes a pose for one and does nothing with it, and
dragging its position wrote the file and moved nothing. It is a separate control
because `poses` is a map beside `materials` rather than a field inside the domain, and it is the
one write that **creates** instead of replacing — no shipped scene states `poses` at all, so
moving anything means writing a key the file does not have. Three levels can be missing and each
is filled in where it is found; a scene that gains a pose keeps every other byte.

Every one of those is a byte splice rather than a re-serialise, so the file keeps its own
formatting and key order. `editor_core::set_number`, `set_text`, `set_pose`, `add_domain`,
`remove_domain`, `add_entry`, `duplicate_entry` and `remove_entry` are where that lives and where
the tests for it are — including one that adds each of the twenty-two kinds and removes it again,
and one that adds to, copies within and removes from every list in every shipped scene, and both
assert the file is byte-identical to what it started as.

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

**View → Frame selection** (and **frame this** in the inspector) re-frames the viewport on the
selection: its box becomes the framing every projection, handle, label, probe, shaded mesh and the
scale bar share, so nothing moves relative to anything else and what is far from it goes off
screen. **Fit view** returns to the whole scene. It used to move only the focal length and keep the
whole scene's framing, on the argument that re-centring would move everything else — which is not
true when every reader takes the same framing, and the old way could not show scene 14's 20 mm bar
at all: it is `7.5e-14` of a framing spanning 2.7e11 m, both its ends narrowed to one `f32`, and
the fit, measured from a centre 2.0e10 m away, put that one point on the top edge of the viewport.
Framed on itself, its outer end lands at 0.85 of the half-frame and it spans 0.76 of the viewport.

## A run with nothing to place is drawn as its readings

Five shipped scenes — three thermal networks and two compartment models — have no panel: their
domains have nothing that is *at* anywhere. The viewport said so in a sentence and stopped. It now
draws the run's readings against time, one plot per unit and stacked on one time axis, with a
cursor at the frame on screen; **View → Readings on log axes** puts each plot on decades, where a
compartment's terminal phase is a straight line.

`editor_core::chart` is all of it: the lines, the axes, the ticks and every pixel coordinate. The
desktop paints its strokes with egui and the browser with a 2D canvas, and neither computes where
anything goes. `tests/readings_have_a_chart.rs` holds it to three things:

- every point is a frame's reading, bit for bit;
- a log axis leaves out a zero rather than pinning it to the floor;
- the slope the log chart *draws* for scene `34`'s blood, read back off two of its vertices a day
  apart, is the body's slowest mode — an eigenvalue of the scene's own perfusion matrix found by
  bisection — to `1e-6/Δt`, the rounding of the seven significant figures the run's JSON carries.
  It was written at 1e-9 across two adjacent frames first and failed at 3.2e-6; the file, not the
  chart.

### And a compartment model, as its body above the chart

A chart says how fast; it does not say where the drug is. `editor_core::diagram` draws each
`compartments` domain as a disc per compartment — area its volume, colour its concentration on one
logarithmic scale for the whole run — joined by a vessel per perfusion and link, with the depots a
dose arrives by and an arrow out of every compartment that clears. The diagram takes the top of the
viewport and the chart the rest, and playing the run moves both.

**Its shape is read from the scene, not from `Placement::marker`.** A marker is one position per
domain, which is the wrong unit: scene `34`'s blood, liver, lean tissue and fat are one domain. The
run carries no topology either, so the layout is computed here from `DomainSpec::Compartments` —
which this crate, as the editor's authoring half, may read — and every value comes from the
readings. `tests/a_body_is_drawn_from_its_readings.rs` holds it to the scene and to `A = C·V`:
at every frame, the discs' concentrations times the scene's volumes add to the `in the body` the
domain reports, to the seven figures the run's JSON carries.

**A head on each vessel points the way its drug is going**, and its length goes with the flow. The
flow is the domain's own reading — `"blood to fat"` in mg/h, `Q·(C_b − C_t/K_p)` — and not
recomputed here, so there is one copy of the formula and the drawing is checked against it rather
than against a second. The closed form it answers to is the tissue's: one joined only to its
blood and clearing nothing rises while its head points in and falls once it points out, so the
frame a head turns is the frame the tissue peaks, and the test reads both off the run.

## Results under the view

A run on screen was a picture, a slider and a readout under the pointer. What a person running a
simulation tool asks next is a number, so once a run is on screen a panel opens under everything
but the tree, with five tabs:

- **Region**: the largest and smallest value and where they are, the mean by volume, the integral
  and the volume. It reads the domain or a box of it, and charts the max, mean and min over the run.
- **Slice**: the field on a plane across x, y or z, in the view's own colours and with a colour
  bar. Under the pointer it gives the value read off the field, not off the pixel.
- **Line**: the field along a straight line, against distance.
- **Point**: the value at a point now, and at every frame.
- **Readings**: every reading the frame reports as a table, and the one picked charted against
  time.

What a tab reads is drawn in the viewport too: the plane, the line, the probe or the box.
Positions are in the domain's own frame, in the unit its size is said in.

`editor_core::post` computes all of it and the shell computes nothing, so the browser editor can
draw the same numbers from the same calls. What it promises, and what holds each promise:

- **Between samples, trilinear; at a sample, the sample to the bit.** A neighbour whose weight is
  zero is not multiplied in, so a probe on a sample next to a hole reads the sample.
- **A hole is a hole.** A point whose neighbours include one has no value, and a line through it
  breaks. A probe on the edge of a part would otherwise read the average of the neighbours that
  exist, which is a number no cell has.
- **Past the last cell centre, that cell.** A cell-centred field has no sample on its walls, and
  extending a line past the data can overshoot what any cell holds.
- **Each lattice integrates by its own rule.** Cell averages use the midpoint rule and samples on
  the faces the trapezoid rule. The linear test cannot tell the two apart, so
  `each_lattice_is_integrated_by_its_own_rule` holds each to its exact error on `x^2`:
  `-L^3/(12 n^2)` for the midpoint and `+L^3/(6 m^2)` for the trapezoid.
- **An insulated block keeps its mean, and a leaking one does not.** Scene 15 is a 9 mm aluminium
  cube at 20 °C with one cell of 729 at 60 K above it and no face cooled, so its mean is
  293.15 + 60/729 K in every one of its 89 frames. Its hottest cell can only cool and its coldest
  only warm. The test holds the mean to 5e-5 K, which is the seven figures the run file writes
  and not a property of the solver. It holds the other two with no slack, because rounding to a
  fixed grid keeps order. The wrong rule fails it by a distance: trapezoid weights on these cell
  averages put the mean 0.035 K high.

  As shipped, the hot cell is in the middle, and over 6 ms its heat barely reaches the walls.
  Cooling all six faces at 25 W/m²/K left the written mean where it was. So the test also puts
  the cell in a corner, where the mean holds to 3.2e-6 K, and cools that corner's face at
  1e4 W/m²/K, where the mean falls 1.6e-3 K. The test asks for both.

The panel is held to the domains' own numbers rather than to `post` alone. On scene 29's bracket
the region's max, min and mean match the block's own `peak`, `coldest` and `mean` readings. The
domain computes those over the 4100 cells the STL filled, and the panel counts the same 4100 with
the 4012 holes left out. On scene 15 the probe opens at the middle cell and reads the region's
maximum.

**What sweeping it found.** Every shipped scene on every tab, and the built-in scene at nine widths
with and without an isosurface, was checked for strings drawn over each other or cut off by the
window. That sweep found:

- **A crash.** Scene 24's line through its module is 488.5735 K to within 1.7e-13: a field flat in
  the file, read between its samples. The chart fitted its axis to that, the first tick's count
  was 9.77e15 (past 2^53, where adding one changes nothing), and the tick loop pushed one tick until
  it asked for 34 GB and aborted. An axis no wider than its numbers' rounding now has one tick, and
  `graph` draws a range under a billionth of the value as flat.
- **Unreadable ticks.** Scene 25's ticks read `435.018200000000`, 106 points wide, across the
  outliner. The protein's ticks read `2.18e11` four times over and its time axis began
  `0.000000000000`. Tick labels now carry as many figures as tell two ticks apart, and zero is `0`.
- **Overlaps the dump counted and nobody would see.** Sixty were reported on the line tab, each a
  wrapping label's box covering the controls before it. The panel's rows lay their labels out whole.
- **Layout faults.**
  - A table of three readings left the panel's row 64 points tall, under the chart's 100-point
    floor, so the chart was not drawn.
  - A panel with only a sentence in it shrank to its floor.
  - At a 700-point window the region's chart beside its statistics laid its time axis's five
    labels across 64 points.
- **Collisions in the viewport, not the panel.** Two domains on one corner had their names drawn on
  one spot: scene 25's structure on its block's own grid, scene 18's two pucks. With the panel
  taking the bottom of the window, the built-in room's corner note fell on the draw count. Names now
  stack, a note that repeats one already there is drawn once, and both stay above the readouts'
  bands and on the view.
- **Height.** At short heights the panel left a view whose readouts overlapped eight times, so it
  keeps the view 320 points tall. When even its own floor does not fit, it goes, and the status
  bar says so.

**What three reviewers found**, each reproduced before it was fixed.

- **A pass that could not fail.**
  - The insulated test, as above.
  - A hole probed where it was not one of the eight corners.
  - A width sweep that passed over a panel drawing nothing.
  - The line tab's `0` and `9` read from its own controls, not its plot.
  - An integral checked against the mean it defines.
  - A label count that always found the axis title.
- **Absent and said nothing.**
  - A chart with no room for it.
  - An unreadable field dropped, so the panel said the run had none.
  - A plane along a bar.
  - Names held above a short view and clipped.
- **Wrong, quietly.**
  - A pick in the domain menu undone the next frame by the outliner's selection.
  - A lone tick printing 0.5 s as `0`.
  - A placed domain's positions not saying they are its own.
  - Two axes at the ends of the doubles that never returned.

## Pointing at the model

A cooled face, a held one, a pad where a part is bolted and a box of cells that dissipates were
set by typing a word or six indices into a form, and checked by counting cells in one's head. A
simulation tool is set up by pointing at the model, so now the view does what the selection asks:

- **An entry that says where.** Select a cooled face, a held or pressed face, a heat source or a
  region in the tree, and its domain's faces toward the eye are named on them. The face, patch or
  box the entry states is drawn on the model. A click on a face sets the entry's `face`. A
  shift-drag across one draws a patch on it for a cooled face, and the footprint of a heat source
  or a region, which keep their depth.
- **The results panel.** On the Point, Slice, Line or Region tab, the domain the panel reads is
  what the pointer is over. A click puts the probe there, takes the plane through the point, or
  sets a line's start and then its end. A shift-drag draws the line, or the region's box over the
  rectangle it spans.

Under the pointer, the face and cell or the coordinates it is over, and what a click and a
shift-drag would do there. Plain drags still turn the view.

`editor_core::pick` does the arithmetic and the shell only paints, so the browser editor can make
the same calls. What it promises, and what holds each promise:

- **What is drawn is what is hit.** The ray is `viewer_core::Camera::ray`, the projection read
  backwards through the same rotation. Points along it project back onto the pixel it left, to
  5.7e-14 over 129 204 points and a sweep of cameras.
- **A face is the domain's own, however it is turned.** The ray is taken into the domain's frame
  by undoing its pose, so `x-max` is the face the file means. A quarter turn about z sends the
  box's `-y` to the world's `+x`, and a third of a turn about (1, 1, 1) sends `y` to `z`. Both are
  held, with the point found put back on the ray.
- **A cooled face states the area it covers.** The block divides `area_cm2` among the **solid**
  cells on the face or patch, so that is what a pick writes: solid cells times a cell's face. On a
  plain block it is the format's own new entry to the bit. On scene 29's bracket it is the 410
  solid cells of 676 that the solver's own account measured on `z-min`, read through the void map
  the build hands back. A face or patch with no solid cell is refused, because it cools nothing.
- **Pointing again changes nothing.** A click on the face an entry already has writes nothing, not
  even an area a person stated by hand. A patch drawn and taken off gives back the file byte for
  byte, on a one-line entry and on scene 29's pad written across two lines.
- **A held face is a whole face.** `from` and `to` on one is an unknown field, so a box drawn on
  it is dropped and the status line says so. Which lists take a face, a patch or a box is the
  format's answer, `templates::SITES`. It is held to the format's types: a patch added to a cooled
  face builds, and added to a held face is refused.
- **Only the aimed-at domain is hit.** Scene 25's structure is its block's own grid. A click meant
  for the cooled face is not taken by the box in front of it, and a click beside it keeps the
  entry rather than clearing the selection.

`tests/a_face_is_picked_in_the_view.rs` drives all of it through `--ui-dump --click`,
`--shift-drag` and `--hover`. It clicks the names the view draws and reads back the file. Each
face of scene 25's module is 64 cells of 2.25 mm², 1.44 cm². A drag between two points writes the
cells the readout named under each. The probe, the line and the region land on the coordinates the
readout gave.

**What three reviewers found**, each reproduced before it was fixed.

- **Passes that could not fail.**
  - The cell under the pointer was checked against the readout, and both come from the same code.
    It is now held to division: 4.6 mm on 1.5 mm cells is cell 3.
  - No test dragged without Shift. A plain drag over a selected face writing a patch would have
    passed.
  - A click "beside" the box could not tell beside from on, or from not reaching the view at all.
    The test now hovers there first to show nothing is under it. A control with a domain row
    selected shows the same click clears a selection.
  - A void map read with `i` and `j` swapped passed every test, because each part the tests used
    is the same with them swapped. Scene 30's crossing busbars are not, and their cells are counted
    by hand from the scene's regions.
  - A highlight's corners were checked at the grid's ends only. One patch inside it is now checked,
    in order.
  - The plane's status line said the coordinate clicked rather than the one stored, so a plane set
    along the wrong axis read back right. Each sentence now says what was stored.
- **Absent, and said nothing.**
  - A part's void map comes from the build, and with no build every cell was counted solid.
    Pointing scene 29's pad at its top in a scene that did not build wrote the grid face's
    27.04 cm². A block with a part or a `void` region is now refused an area until it builds.
  - An entry in a hidden domain took clicks and did nothing, and a click that missed said nothing.
    Both say so now.
- **Scene 29's top film cools nothing, and has not since it shipped.** The bracket fills layers 0
  to 9 of a 12-layer grid, 410 cells each, and leaves layers 10 and 11 empty. The block cools its
  `z-max` face at layer 11, which holds no solid cell, so the 16.5 cm² still-air film acts on
  nothing. The format refuses a patch of void but not a whole face of it. The scene is left as it
  is here, since its numbers and the figure drawn from it would move. The editor says so instead:
  the entry's form reads "no solid cell of bracket is on z-max, so its film cools nothing".

Not yet: the browser editor's view does not point. A contact's plane between cells is not
something a click sets. A click on a domain the results panel is not reading selects it, and
the probe takes a second click.

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

**It sends one frame at a time, then the settled whole.** It used to send the whole run so far
after every frame, which is quadratic in the frame count: the cavity scene's 1600 frames
serialised about 31 GB to reach a 38.5 MB run, on the editor's run thread every time somebody
pressed run, and the test that streams every scene took twenty-six minutes. A frame's JSON
depends only on that frame, so the editor appends each one, and replaces the lot with the
`Streamed::Whole` at the end.

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

The forms are structured editing *over* the text, not instead of it: a rename, a reference
followed from one domain to another, a key no form offers — those are the Source tab's, and
`deny_unknown_fields` plus the version check are still what stand between a typo and a silently
different scene. A string that names another entry — a link's `from`, a dose's `into` — is free
text, not a menu of the names the scene has; that is the next thing this wants. No shadows, no ambient
occlusion and no anti-aliasing beyond whatever the context gives: the lighting is a key, a fill and
an ambient, which is the least that shows a shape. No file dialogs. Undo covers the editor's own
edits and leaves typing to the text box, which is the division described above rather than a gap.
Each of those is worth adding only after using this one reports back what actually
chafes; that method has a name here, and `FRICTION.md` is where its findings go.
