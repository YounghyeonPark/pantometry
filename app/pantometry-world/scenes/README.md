# Scenes

Thirty worlds described as data, covering **all eleven** of the library's domains — twenty-nine of them
one question at a time, and one that is actually a world. The count above this line used to read
"twenty-seven of them ... and one", which is twenty-eight and was never the number of files in this
directory; **seventeen** state a single kind of domain and thirteen state more, so no split of them
into "one physics" and "a world" was ever going to be a count of anything. Nothing here is Rust: the
physics, the resolution, the coupling and the run length are all in the file, and the same
binary runs all of them.

```sh
cargo run --release -p pantometry-world -- scenes/03-room-pulse.json
cargo run --release -p pantometry-world -- scenes/03-room-pulse.json out.svg
cargo run --release -p pantometry-world -- --check scenes/03-room-pulse.json
```

`--check` parses and builds without running: the format version, the domain names, a `tracks`
pointing at a node the scene defines. It reports a parse failure as `file:line:column` with the
keys that were expected, which is what an editor puts a squiggle under. CI runs it over every
scene, because it would otherwise be the one entry point nothing exercises.

Every file carries a `format` number, and **absence means 1** — which is what all thirty here
are: nothing has yet changed what an existing key means. A version this build cannot read is refused
rather than half-run: `deny_unknown_fields` catches a key that was *added*, but not one whose
meaning changed, and that is what the number is for.

No second argument prints the numbers and checks them. A second argument writes an asset, and
**the extension chooses which** — because a run has several shapes and only one is a picture:

| | |
| --- | --- |
| `out.html` | A report that **picks a view per domain from its shape**, and opens in a browser. A 3D field gets two: a rotatable render and every slice |
| `out.svg` | A filmstrip: every frame on one page, one colour scale throughout so frames compare |
| `out.csv` | Every domain's scalars over time, one row per frame, units in the header |
| `out.json` | The frames themselves — fields as grids, bodies as positions in space, readings beside them |
| `out.gltf` | The last frame as **surfaces** — a field's boundary, a body's sphere, with normals — for Blender, three.js or any glTF reader |
| `out.usda` | The **whole run** as USD: geometry, colour and every domain's scalars, animated on a timeline, for usdview, Omniverse, Houdini or Maya |

`.csv` is the one that reaches the domains a picture cannot. Twelve of these thirty scenes have
a domain with no field and no bodies, and for several the scalar *is* the result: `13` is about a
winding whose resistance follows its own temperature, and it drew nothing at all. As a table it
shows the feedback directly — 12.46 W at 25 °C rising to 16.01 W at 99 °C, with the resistance
going 3.11 to 4.00 Ω beside it.

`.html` is for someone who can state a simulation and does not want to decide how to draw it.
The shape of the data makes the choice, not the domain's name: scalars over time become a line
chart, a 1D field a profile over a ghost of the whole run, a 2D field an animated heatmap, and
bodies a rotatable depth-sorted scene. A new domain gets a sensible picture without the reporter
learning about it — which is the same reason `Domain::as_field` exists.

One scale per panel, fixed across every frame, in all four views. A frame that rescales makes a
quantity *look* constant while it changes by orders of magnitude, and that is the one thing a
picture of a simulation must never do.

And the picture has numbers on it. Every spatial axis is in metres — a room is 4.4 m by 3.1 m and
not "61 x 43" — because a field carries the box it was sampled over; hovering reads back the
sample under the cursor, its index, its position and its value; the scalar chart groups series by
unit so two temperatures forty kelvin apart no longer draw the same line; and the colour scale is
built in CIE LCh so a larger value is never darker than a smaller one. The scale that shipped
before this ran backwards for 89 of its 255 steps.

**All four writers are `pantometry-view`, and none of them is this application's.** This crate is
`publish = false`, so while they lived here they were unreachable — a consumer who could state a
simulation could not draw it. `pantometry_view::{html, svg, readings_csv, to_json}` takes the frames
`pantometry_scene::capture` produces, and everything the table above describes is available to any
program without going near a scene file.

`.gltf` is the one that leaves this workspace. **Twenty-three** of the thirty scenes have geometry
to export — bodies, ray paths, a 3D field as its cell centres — and the other **seven** are
**refused with a reason** rather than written as an empty scene: a 1D or 2D field is a graph, not
something to put in a 3D viewer, and the message says which panel and why.

Those two were nineteen and ten, and had been wrong by three in each direction for long enough that
nothing records when they stopped being right. Nothing counts them: `counts_in_prose.rs` guards the
*total* half of this sentence and says in its own doc why it leaves the numerator alone — a
numerator is a different count with a different source, and this one needs every scene exported to
establish. So it is a **release-time** count now, with the command in `RELEASING.md`, rather than a
number somebody remembers. Measured 2026-09-01.

The recommendation on record is to export into the rendering tools rather than rebuild them. This
is that, and it cost no dependency: glTF is JSON with the binary base64'd inside it, which is the
same reason SVG was chosen over a raster format.

Nothing generated is committed — with one exception, and it is named so that this sentence stays true: `tools/presets` writes `presets.rs` and the chooser's tiles, which are build *inputs* rather than outputs and cannot be produced without a GPU.

`cargo run -p pantometry-world -- --emit-default scene.json` writes a starting point.

## Sound — `pantometry-acoustic`

| Scene | What it shows |
| --- | --- |
| `01-room-mode` | The (1,1) mode of a 4.4 × 3.1 m room: the whole field breathing in and out together, one nodal line each way |
| `02-room-higher-mode` | The (3,2) mode of the same room — more nodal lines, a higher note, and a shorter run to catch it |
| `03-room-pulse` | A clap near a corner. No standing shape, so it travels, reflects off all four walls and interferes with itself |
| `16-a-room-with-a-ceiling` | The same room with a **ceiling**, released in its oblique (1,1,1) mode. 97.46 Hz, and the peak rides `\|cos(2πft)\|` to within 0.3% at 23 nodes across |

`16` is the scene `01`–`03` cannot be. A floor plan does not have the floor-to-ceiling mode
inaccurately; it does not have it. At 2.4 m that mode is 71 Hz, and the number of modes below a
given frequency grows as `f³` rather than `f²` — which is why a real room's resonances merge into
a hiss where a two-dimensional model keeps them separable much further up.

It costs what the third dimension costs: 23 × 17 × 13 nodes rather than 23 × 17, and `√3` in the
Courant limit rather than `√2`. That is the trade, stated rather than hidden, and `Room` remains
the right model for a floor plan.

## Heat — `pantometry-thermal`, and the two ways domains meet

| Scene | What it shows |
| --- | --- |
| `04-heater-and-bar` | A plain channel. The heater's joules have no *place*, so they land in the bar's first cell and conduct along it |
| `05-beam-on-bar` | A shared boundary. The beam's joules do have a place, so a spot forms where it points and then spreads |

Those two are the same physics told twice, and the difference between them is the whole
argument for `Interface` and `Flux`. A conservation audit passes either way; only the picture
tells them apart.

| Scene | What it shows |
| --- | --- |
| `15-a-hot-spot-in-a-block` | One cell of a 9×9×9 aluminium block starting 60 K hot, and the spot spreading in **three** dimensions. The scene that a one-dimensional bar cannot express: heat going sideways is the whole job of a spreader plate |
| `19-a-coating-stops-the-heat` | A **heated face** in a 9×9×18 block that is aluminium for half its depth and **borosilicate** for the other half. 150× the conductivity leaves the metal nearly isothermal and the glass barely touched, so the entire temperature drop sits on one face — and the test asserts the largest cell-to-cell step along z is exactly the interface one, which is what a block with the coating quietly not applied would fail. Measured: the metal levels to **3.03%** across its own half, the interface step is **9.4656 K of a 19.1152 K rise** with the next steepest at 7.3315, and four cells into the glass is 0.14% of the metal. **It was one cell at +60 K** — 0.145 J over 1 458 cells — so every one of those claims was true and was a claim about 0.08 K; heating the whole face is eighty-one times the energy. The pulse is a `regions` entry rather than a `hot_spot` now, and that is what lets the resolution sweep run at all: a hot spot is one cell, so refining halves its physical size, and `verify` refuses. It measures 0.106% on the peak |
| `20-melting-a-block-of-ice` | A hundred watts into 1331 mm³ of ice at −10 °C. The mean holds at **exactly 0.000 °C** for 3.6 of the 5 s while `melted` climbs a straight line at 326.9 mm³/s — which is `P/ρL` and has no rate constant in it — and then the leftover 67.9 J warms what is now all water by 27.126 K, against 27.126 predicted. The scene a domain with no latent heat cannot express: without it the block would sail through zero and finish near 190 °C |

| `21-a-wax-thermal-buffer` | Twenty watts into 10 648 mm³ of **n-octadecane**, which is not in the catalogue and is written out in the scene's own `materials` block. It holds 28.150 °C from 11 s to the end while `melted` climbs 100.6968 mm³/s against a `P/ρL` of 100.6968, agreeing to **8.3e-15** — the same closed form `20` measures, now with the density and latent heat coming out of the *file*. The 1200 J of reserve runs out at 60 s and the melting stops dead at 5354 mm³, 50.3% of the block. The energy accounts to `1.1e-14`: 136.6176 J to warm 8.667 g of wax by 8.15 K and 1063.3824 J to melt half of it, which is all 1200 J and not most of it |

`20` is the one that needs latent heat, and its plateau is the whole of what latent heat *is*: a
temperature that stops moving while energy keeps arriving. Ice is the only entry in the catalogue that
melts, and a block that cannot melt does not report a `melted` column at all — so the column's
presence is itself the check that `"material": "ice"` reached the domain.

| `22-wax-in-an-aluminium-matrix` | The same wax as `21`, now four fifths of a **composite** whose other fifth is aluminium — the scene's `composites` block mixes a declared substance with a catalogue one. It holds the same 28.150 °C, and `melted` climbs **125.8710 mm³/s against 100.6968** for the pure wax: exactly `1/0.8`, a ratio with no material property left in it. 1200.0000 J = 151.2823 warming + 1048.7177 melting, to `1.7e-14` |
| `23-a-part-radiating-to-its-lid` | A hot part sealed in a housing, surrounded on five sides by **nothing** and facing a cooled lid across a two-cell clearance. The part has no conducting face to anywhere, so its heat leaves by the parallel-plate exchange across the gap and by the air around it — a two-body lumped system with a closed form. It sheds **101.93 K in 600 s against 102.05 K** from an RK4 integration of `σA(T₁⁴−T₂⁴)/(1/ε₁+1/ε₂−1)`, `hA(T−T∞)` and the air on each body’s wetted faces, agreeing to **0.12%**. **It was a vacuum and called itself a housing.** `cooling` reaches a block’s six outer faces and nothing else, so a part rasterised inside one had no convective path at all; the scene was honest about being radiation-only and was 65 K wrong for anyone reading it as a housing — 471.22 K with still air in the cavity against 536.22 without. The answer is still a **surface-finish** answer: the same geometry in polished copper, ε = 0.04, sheds far less across the gap |
| `24-a-power-module-junction-to-ambient` | A silicon die dissipating **45 W** through the stack under it — a 100 um solder joint, a 0.63 mm alumina DBC, copper, an interface material and a 3000 W/m^2K cold plate. The scene the format could not state at all before `dissipation`: every other source here hands watts to the *bus*, which carries an amount and no location, so a die's heat would have spread over the baseplate as fast as over the die and there would be no junction temperature to read. At steady state it is a resistance stack, and the test writes every term of it out from the geometry: **227.15 C measured against 227.17 C** for a 4.1594 K/W path, agreeing to `1.1e-4` of the 187 K rise. The scene **names** that place now: it read the block’s `peak` — the hottest cell anywhere, which is the junction here only because the die spans the whole cross-section — and states `probes` for the junction and the baseplate instead. The two agree to the bit, which is a check of its own: a `peak` comes from the cells the domain holds and a probe from `ScalarField::at`, which interpolates. **Three of those terms could not be written until `contact` existed.** The solder is 100 um and the ceramic 0.63 mm against a 1.5 mm cell, and a region is a box of cells — so 1.5 mm was the thinnest either could be written as, and the solder carried fifteen times its own resistance. That is a **floor** and not a discretisation error: the resolution sweep doubles a region's bounds with the counts, so refining does not approach anything, and no grid this scene can hold reaches 100 um either — that is 120 cells a side. And the **largest** resistance in a real junction-to-ambient path was absent altogether: a baseplate is bolted to its cold plate through an interface material, 1.3889 K/W here against the 3.1379 the whole stack used to be. The two errors pushed opposite ways and left something plausible — the scene answered **32.5% low** on the number it is named after. The mounting is 33% of the path now and the film 56%; the ceramic is 0.18 K/W and the solder 0.01 |
| `25-what-140-kelvin-does-to-the-solder` | The same module as `24`, with a **structure** on the same grid whose stress-free strain is the block's temperature — the first scene that couples two physics on one mesh, and the first to reach `pantometry-elastic` at all. The failure mode of a real module is not its temperature: it is the solder, fatigued by silicon at 2.6e-6 per kelvin sitting on solder at 2.15e-5. Assembled at its reflow temperature of 217 C, so it is **already strained before it is switched on** — cold, the solder wants `2.15e-5 x (40 - 217)` = `-3.8055e-3` exactly — and it *relaxes* as it heats back towards where it was built: 0.2314 J of strain energy falling to 0.0274 J, **8.4x**, monotonically. **And it is outside the model it is solved with**, which nothing said until `verify` learned to check. SAC305 yields at 30 MPa against 41 GPa — a strain of 0.0732% — and this asks for 0.3806%: **5.201×**, worst in the solder layer at `z = 5`. `pantometry-elastic` is linear with no plasticity and documents what that costs: past yield it returns a displacement that is arithmetically correct and physically meaningless, and nothing in the answer says which. The scene is kept because the **coupling** is what it demonstrates and the coupling is right; a real joint would have yielded and crept instead |
| `26-poiseuille-in-a-cooling-channel` | Water driven down a 2 mm channel by a body force, the first scene to reach `pantometry-fluid` — the domain this workspace's own docs call the hardest to trust, because *it looks like a fluid* is the easiest wrong answer in computational physics to accept. So it is written around an exact solution, and against the **discrete** parabola rather than the continuum one: a no-slip wall imposed by reflecting the first cell makes the linear interpolation vanish there, and a parabola is not its own linear interpolation, so the settled mean is `(gh²/12v)(1 + 2/n²)` exactly. Measured `6.691655e-3` against `6.691982e-3` m/s. The tolerance is the **startup transient**, predicted rather than chosen: 4 s is 9.9 time constants of `h²/(pi² v)`, leaving `5.0e-5` of the answer, and the run is off by `4.9e-5` |
| `27-a-cavity-ringing-at-its-own-frequency` | A 120 x 120 mm vacuum box seeded in its `(1,0,1)` mode, the first scene to reach `pantometry-em`. A resonance is a property of the box, so `f = (c/2)sqrt((l/a)^2 + (n/d)^2)` is what it is checked against: **1.766304 GHz measured against 1.766544**, off `1.4e-4` — and the bound is **Yee's own dispersion**, `(k dx)^2/24` = 1.4e-3, with the measurement required to come out *under* the continuum because a wave on this grid travels slow. The energy the audit watches is **not** `1/2 eE^2 + 1/2 uH^2`: E and H are half a step apart, so the naive sum swings 7.4% about the quantity leapfrog actually conserves, which holds to the bit across every frame. `div B` stays at 1.9e-13 — an identity of the discrete curl, not a convergence |
| `28-an-eigenstate-that-does-not-move` | An electron in the third state of a 10 nm hard-walled well — the last of the eleven domains to reach a scene, and the strongest claim a solver can be given: **an eigenstate is stationary.** Energy, position expectation and norm are asserted frame by frame and none of them moves at all. Checked against the **discrete** Hamiltonian's exact eigenvalue `(2hbar^2/m dx^2) sin^2(n pi/2(N+1))`, and the gap to the continuum `n^2 pi^2 hbar^2/2mL^2` is `theta^2/3` exactly — 1.8320e-4 measured against 1.8322e-4 predicted, which is the grid being measured as the grid instead of as the physics |
| `29-a-designed-bracket-becomes-cells` | The first scene whose geometry comes from a **file** rather than from numbers in the JSON: an L-bracket named as an ASCII STL, rasterised onto 2 mm cells, carrying a module’s **20 W** from the tip of one arm, round the corner, into a bolt pad at the tip of the other. The check is the **outline**, not the mesh — the shoelace area of the seven-point profile times the 20 mm extrusion is 33 000 mm3, which 4 125 cells would hold and 4 100 do, **-0.61%**. Comparing against `Mesh::volume` instead would only say the crate agrees with itself, since that is the divergence theorem over the same triangles the rasteriser read. **It cooled over its whole footprint into still air** until `verify` learned to say when a grid is idle: every cell shed where it stood, so 4 100 cells rasterised from a shape all held the same number and the field spanned **0.067 K** at a Biot number of 8.6e-4. The route is the only reason the shape is in the file, and there was no route. It holds **52.4420 °C at the module against 30.5449 at the bolts** now — 21.90 K across a path the outline decides. It starts at 44 °C rather than 20: from cold, most of the run charged the bracket’s mass rather than establishing the gradient, and nine test binaries walk every shipped scene, so 900 s cost nine times over. 200 s from near the answer is within 0.09% of the settled peak. That needed `cooling` to name a *box* on a face: stating a smaller area does not do it, because the area is divided among the cells on the face, so a tenth of the area is a tenth of the conductance spread over all of it |
| `30-two-phases-crossing-at-a-clearance` | Two blackened copper busbars, identical but for their current, one turned a quarter turn and lifted to a 4 mm clearance — an arrangement where *not conducting* is the design requirement. **The first scene to state a `poses` entry**, and it exists because the other twenty-nine could not fail: under the identity a domain’s own coordinates and the world’s are the same thing, so three separate consumers dropped the placement in turn and every scene agreed with all of them. **It was two domains and they could not see each other.** Two surfaces at 335.6 K and 319.0 K, ε = 0.9, four millimetres apart, exchange **6.5711 mW** across the 8 × 8 mm patch where they cross — 7.64% of what the cooler bar dissipates — and two `Solid3D` domains exchange nothing but bus totals, which carry an amount and no location. As one block it is `find_gaps`, the same pairing `23` is checked on. The closed form is a **pair**: each bar sheds convectively and radiatively *and* trades with the other, solved by iteration from the constants, giving 22.0246 and 6.2523 K against a measured **22.0218 and 6.2493**. Uncoupled they would sit at 22.4333 and 5.8139, so the exchange is forty times the tolerance the balance is asserted to. Each bar is a **lump** and that is the physics: with every watt leaving at the ends, a 32 mm copper bar 8 mm square at 0.344 W varies by `P·L/(8kA)` = 53.6 mK, and it measures 0.585 |

`22` is the pair to read with `21`, and the ratio between them is the whole point. Diluting the wax to four
fifths by volume makes a cubic millimetre of buffer hold a fifth less latent heat, so the same twenty watts
clear it a quarter faster — and `1.25` is `1/0.8` with the density and the latent heat cancelled out. The
*engineering* reason for the aluminium is the other number: 5 W/m·K against the wax's 0.358, fourteen times
better at getting heat in and out of the thing that stores it, which is what a metal-matrix phase-change
buffer is for.

**That conductivity is inert in this scene, and the sentence above was worth checking rather than
believing.** The block has no cooling, one material and a spatially uniform source, so every face flux is
identically zero and a uniform field is an exact fixed point of the update — for *any* conductivity.
Measured: 0.5, 5.0 and 30.0 W/m·K, a sixtyfold range inside the bounds the format enforces, produce
**byte-identical** output in every frame, and so do the same three with the reserve cut so the plateau
cannot mask them. What the scene does exercise is the mixture's density and latent heat, and those are
live — `melted` differs from `21` by 23% and matches the first-principles rule `ρL = φ_wax ρ_wax L_wax` to
8e-9. The conductivity is live too, checked separately by giving the same composite a half-block step and
fitting the slowest Neumann mode: `α = k/ρc_p` predicts a decay rate of 5.8487e-2 s⁻¹ and the solver gives
5.8503e-2, **+0.03%**. It is the *scene* that does not ask the question, not the solver that cannot answer
it. A cooled face would make it ask; that is a change to a shipped scene's physics and has not been made.

That 5 is **the caller's choice and the format checks it**. No single conductivity exists for a composite
without knowing its microstructure, so the scene has to say which, and a value outside the Voigt and Reuss
bounds — 0.4473 to 33.6864 for this pair — is refused with both bounds and the tighter Hashin–Shtrikman
pair in the message, because a file has nowhere else to learn them. 5 sits inside HS's 0.6244 to 24.185,
which is where an aluminium foam at a fifth by volume belongs.

The latent heat dilutes by **mass** and not by volume, and that is the trap the scene exists to make
visible: four fifths of the volume is 54.67% of the mass, so the composite stores 133.4 kJ/kg and not
195.2. Using the volume fraction there would put the plateau 46% wrong.

`21` is the one that needs **no catalogue entry at all**, and it is the format's first `materials` key.
A wax is what a phase-change buffer is actually made of, no library will ever ship every material
somebody has a datasheet for, and enumeration was never going to close that gap — so a scene can write
a `Substance` out in full and use it exactly as it uses `"ice"`. Nothing downstream can tell which one
it got.

Three mistakes in a declaration are refused rather than accepted: an impossible substance, a name that
shadows a catalogue material, and — the one that is not obvious — a declaration **nothing goes on to
use**. `material` on a block is optional and defaults to aluminium, so declaring a wax and forgetting to
name it gives a block of metal that runs, audits and renders while answering about the wrong substance.
That is the same failure `19`'s empty region has, one level up, and it gets the same treatment.

This scene also carries the first `liquid` block: the melt conducts 0.152 W/m/K against the solid's
0.358, so the buffer insulates itself as it melts, which is the physically interesting thing about a wax
and something a single conductivity cannot say.

`19` is the one that needs `Solid3D::fill`, and it is the format's first `regions` key. A region is
**half-open** — `to` is one past the last cell — and one that selects no cells is refused rather than
ignored, because a mistyped bound is otherwise invisible: a block of one material runs, audits and
renders, and nothing anywhere says the coating was not applied.

`15` is the first scene whose field is a volume, and the report draws it as **every z-slice at
once** rather than as one plane with the rest behind a slider — because a viewer who never
touches a slider would see a picture of a solid that was really a picture of one plane through
it. The filmstrip has no room for that, so it draws the middle slice and *says* `z-slice 5/9` in
the label.

Its claim in `tests/scene.rs` is the one only three dimensions can make: the neighbour one cell
away along **z** is exactly as warm as the one along x. A model that resolved a plane and stacked
it, or that used the wrong spacing on one axis, fails there and passes everything else.

| Scene | What it shows |
| --- | --- |
| `11-motor-thermal-network` | A copper winding, electrical steel around it and an aluminium housing: three materials, two joints, and the drop across each. **It starts under load** — 36 W for 300 s and 12 W after, three times the loss, which is 1.73 times the current — and that is the whole point of it. The winding peaks at **72.82 °C** and settles back to 59.64; the same scene at a constant 12 W climbed monotonically to 55.04 and reported that, which is **17.8 K low** for choosing an insulation class. A first-order network under a constant source cannot overshoot at all, so the peak being above the end is a claim `stages` is required for. The shape is the winding’s own time constant: 62 J/K across 0.9 W/K is 69 s, against 216 s for the stator, so the winding tracks the load while everything behind it is still climbing |

The first of three scenes with **nothing to draw** -- `11`, `12` and `13`, which is the list
`tests/scene.rs` holds as `NOTHING_TO_DRAW`. A network's nodes have capacities, not positions, so
`as_field` declines to invent a continuum and the renderer has no panel to make — the numbers
are the output. That is also the point: a `lump` would report the motor as one temperature, and
the winding runs 13 K above the housing at half a time constant in. The thing that fails is not
the thing you can measure.

The scene test carries a matching list, so an undrawable scene has to earn its place with an
explicit check rather than passing by having nothing to check.

## Electricity — `pantometry-electrical`

| Scene | What it shows |
| --- | --- |
| `17-a-busbar-with-a-notch` | A 12 × 5 × 5 mm copper busbar with a notch three cells deep, leaving a **2 mm ligament**, driven at 1 mV. The resistance is **solved** rather than stated: 12.39 µΩ, against 8.28 µΩ for the full section. It is also **4.65% above** the answer the shape has, and that is the interesting part — its sweep would not run at all until a blocked cell learned to refine into its eight children, and what it then measured is that a notch converges at `h^(4/3)` and not `h²`. Current crowds into a re-entrant corner: the conducting wedge subtends `3π/2`, so `λ = π/ω = 2/3` and a quadratic functional of the field goes as `h^(2λ)`. Measured 1.33534, 1.33296, 1.33318, 1.33367 across six grids, bracketing 4/3. The slot is one cell wide here, so the notch roots are one cell apart and the corner is not resolved at all. Refining is cheap — `R·t` is exactly thickness-independent — but the scene *refuses* 0.0156 mm on its own 1e-9 drift budget, so the finest thing it can say about itself still has 0.046% in it |
| `18-an-espresso-shot` | Two 30 mm baskets, identical but for the ring against the wall. Nothing states a flow rate — Darcy’s law is solved on the permeability that 250 µm at `ε = 0.45` gives. **It is a shot now**: 25 s, 10.26 g from a 4.29 g dose, **19.90% extraction** — the middle of the 18–22% a barista pulls to — at 8.31% TDS. It ran for **eight seconds** before, delivering 3.28 g at 4.81%, a quarter of a shot on a scene named for one; every claim it made was about the *comparison* between the two baskets, and Darcy at a fixed pressure gives a constant flow, so nothing failed. The gap makes that ring 2.9× more permeable, so it runs **79% faster** — the flow ratio is 1.786775, which columns in parallel give exactly — and pours **18.34 g at 4.82% TDS against 10.26 g at 8.31%**: nearly twice the liquid at 58% of the strength. The diagnosis is the reading that separates the two hypotheses: an even basket extracts its ring and its core alike at 1.0000, and the gapped one’s ring outruns the core it starved at 1.1539 |

The scene that shows why a field formulation is worth the solve. `ρL/A` is a statement about a
uniform bar, and a bar with a notch is not one — so the file states a *shape* and a material, and
the resistance comes out. The test asserts two bounds rather than the measured value: above
`ρL/A` for the full section, because removing conductor cannot help, and above a naive series
estimate that treats the notched slice as a shorter bar, because the current also has to spread
back out. The excess over the second **is** the spreading resistance, and it has no closed form
for this shape.

The potential is the panel, in volts. What a picture of it shows is where the current is going.

| Scene | What it shows |
| --- | --- |
| `12-winding-heats-a-motor` | The same motor as `11`, with the watts **computed** instead of stated: 62 m of 0.35 mm² copper at 1.75 A |

`11` says 12 W. This one derives 11.93 W from a length of wire, and the two settle within a
fifth of a kelvin of each other — so the guess was a good one, and this scene is what would have
caught it if it had not been. A stated number cannot be wrong, which is another way of saying it
is not a model.

Evaluated at 90 °C, which is worth 27.5%: copper gains 0.393% per kelvin, so the same coil on a
cold bench dissipates 9.35 W. The temperature is a parameter rather than an omission, and it is
one the *simulation* does not set — a domain cannot read another's state inside the step loop,
so closing that feedback is the caller's job and `pantometry-electrical` says why.

| Scene | What it shows |
| --- | --- |
| `13-winding-that-heats-itself` | The same coil with `tracks` set: its resistance follows its own temperature, so it settles 16 K hotter than one held at ambient |

**The only place in this repository where two domains are coupled by hand**, and it is worth
saying why that is allowed. Domains never read each other *inside* the step loop — they meet on
the bus, which carries amounts and not state. This runs between frames, in the code that owns
the simulation. It needed `Simulation::domain_as_mut`, which did not exist: a caller could read a
domain and not write one, so this loop was closable from nowhere at all. FRICTION 18.

The amplification is `1/(1 − g)` with `g = I²R₂₀·α·R_th`, and the test checks it as a ratio
against the same scene with `tracks` removed. Measured 1.281. Convection alone predicts 1.310;
the housing's linearised radiative conductance at 74.6 °C is 0.036 W/K against 0.294 for
convection, and including it gives 1.280 — so the 2.2% gap is radiation stiffening the heat path,
not error.

## A world — four crates at once

| Scene | What it shows |
| --- | --- |
| `14-a-world` | A laser and a lamp both heating one bar, sound crossing a room, three planets orbiting. Five domains, one clock, one bus, one audit at 1e-9 |

Every other scene here is one physics, or two meeting. This is the first with more than two, and
it exists because the crate split's claim is that domains *compose* — a claim verified in pairs
and never beyond. Building it found two things pairwise coupling cannot reach.

**A second consumer of one channel silently got nothing.** `Exchange::take` empties a channel, so
the second domain to take gets zero while every total agrees, because everything published was
consumed. Two plates under one lamp warm at the rate of one plate. That is refused now, and it
had been in every released version.

**A world's tolerance is set by its loosest domain.** The first attempt included the atom box,
which runs at 5e-2 over a few of its own time units; the rest of this scene holds 1e-9 over
0.2 seconds.
Those are not reconcilable, and that is physics rather than a defect — a Lennard-Jones fluid and
a planetary orbit do not share a clock. "Physics for simulated worlds" means *a* world, not all
of them at once.

Two consequences that shape any scene like this. Declaration order is execution order, so a
producer declared after its consumer publishes into a step that has already taken — the audit
catches it as *published but not consumed*, which is how the lamp got moved above the bar. And
the plain channel is a single global resource: many producers may write to it, exactly one domain
may consume it.

## Motion — `pantometry-mechanics`

| Scene | What it shows |
| --- | --- |
| `06-orbits` | Four satellites on circular orbits round one mass, each tilted out of the reference plane. The inner ones lap the outer ones, which is Kepler's third law drawn rather than asserted |
| `07-bouncing-ball` | A penalty contact with a dashpot, losing height each bounce — and a thermal lump taking the joules it dissipates |

`07` needs that lump. Without a consumer the contact publishes heat that arrives nowhere and
the kernel refuses the step, which is correct and is worth seeing once.

## Matter — `pantometry-molecular`

| Scene | What it shows |
| --- | --- |
| `08-atoms-crystal` | 108 atoms at `T* = 0.15`: still on their lattice sites, rattling in place |
| `09-atoms-liquid` | The same 108 at `T* = 1.4`: the lattice is gone and they wander |

Same seed, same density, same box — only the temperature differs, so the two pictures side by
side are melting. The periodic cell is drawn as a wireframe, because it is a real boundary
rather than the edge of a picture.

**Both were frozen lattices until 2026-09-04, and the sentences above were describing nothing.**
The domain is built with `LennardJones::reduced()`, where σ = ε = m = 1 in SI and the only time
scale is `τ = σ√(m/ε)` — one second. Both scenes asked for `duration_s: 6.0e-12`, which is 6e-12 τ:
ten orders of magnitude short of the 0.01 τ the domain itself suggests for a *single step*. Six
picoseconds is where argon melts (σ = 3.405 Å, ε/k = 119.8 K, m = 39.948 u give τ = 2.16 ps), and
that is the number that got written down.

Measured: the mean square displacement grew as exactly `t²` — the ratio between the last frame and
the half-way one was **4.000**, which is free flight with not one collision — and the "liquid" sat
on the crystal's own lattice sites to 1e-11 σ. Every check the scenes had passed throughout, because
the *speeds* were right: equipartition gives `⟨v²⟩ = 4.161` against `3T*` = 4.2, the 321/324 of a
box with its centre-of-mass momentum removed. The initialisation was correct and only the clock was
wrong.

At `duration_s: 6.0` the crystal saturates at **0.075** of a neighbour spacing (Lindemann melts near
0.1) with a growth ratio of 0.875, and the liquid reaches **1.174** neighbour spacings with a ratio
of 2.42, which is diffusion. `every_scene_that_ships_runs_and_says_something_true` asserts both, and
that is what the equipartition check could not see.

## Light — `pantometry-optics`

| Scene | What it shows |
| --- | --- |
| `10-lamp-on-a-mirror` | A 100 W tungsten lamp on an aluminium mirror. The mirror is worse in the blue, so the lamp's *colour* decides how much becomes heat |

`pantometry-optics` has no `Domain` in it — spectra and Fresnel coefficients answer questions
rather than march a state — so this one is written in the application, like the heater and the
beam. A flat reflectance would make the colour temperature irrelevant and the whole spectral
apparatus an expensive way to multiply by a constant, which is why the test compares 2800 K
against 6500 K rather than checking one number.

## `frames` was a physics knob and is not any more

A scene says `duration_s` and `frames`, and the run used to advance by `duration_s / frames`. Each
domain then subdivided that into whole substeps no longer than its own stability limit, so the step
was `window / ceil(window / limit)` — **a function of the frame count**. Asking for more pictures
made the answer better, and nothing said so.

Measured on `15-a-hot-spot-in-a-block`, whose peak reads

| frames | 2 | 11 | 51 | 101 | 401 | 1601 |
| --- | --- | --- | --- | --- | --- | --- |
| peak K | 5.177 | 7.027 | 7.940 | 8.065 | 8.161 | 8.186 |

which is 58% between the ends and textbook first order in the step. It shipped at 11, **14% below
its own grid's converged answer**.

**A stability limit is not an accuracy limit**, and that is the whole confusion. `Solid3D`'s limit
here is 2.414e-3 s; at a fifth of it the scene is still 14% out. `max_stable_dt` says what will not
diverge, and nothing said what would be right.

`window_s` is the separation. State it and the run takes `ceil(duration_s / window_s)` whole steps,
with `frames` choosing only which of them are photographed — so the frame count cannot move a
number, and `the_frame_count_is_pictures_and_not_physics` asserts that bit-for-bit. Leave it out and
the old behaviour stands, which is why every scene written before the key is unchanged.

What tells a scene it needs one is `pantometry verify`. Its window sweep has always measured this;
what changed is that it now divides the shift by **what the run actually did** rather than by the
reading's own magnitude — for a celsius temperature that magnitude is dominated by the 273.15 the
scale carries, so a shift worth 0.966% of everything that happened was reported as 0.212% and read
as converged. Above half a percent it is a finding and carries the exit code. Three scenes tripped
it and were corrected: `15` to 88 frames, `05` to 52, and `27` to 1600.

`27` is worth its own sentence. Its shift goes 0.608, 0.000, 0.916 and 0.467% at one, two, four and
eight times its frames — **not monotone**, because what moves in a resonant cavity is the phase the
sampling lands on rather than a convergence error. The finding still means what it says there, and
it is not an error estimate; `--deep`'s measured order is what tells the two apart.

## Every one of them is run by CI

A scene in this repository is a claim, and one that parses and then produces nonsense is worse
than none at all. `tests/scene.rs` runs all thirty on every commit and asserts one number each —
chosen to be a property of the physics rather than of the file, so it would change if the
library broke and not merely if the scene were edited. Adding a scene without a claim fails
the test rather than passing quietly. CI also runs the real binary on the real files, which is
what a reader will type.

Running a scene is not a weak check on its own: the conservation audit is live for the whole
run, so one that leaked energy, created it, or left it unclaimed on a channel fails before it
draws anything.

## Three numbers that look wrong and are not

**A room is not quite the size you asked for.** `Room::of_air` quantises the height to a whole
number of cells so they stay square and the stability limit stays isotropic. A 3.1 m room at 81
cells across is 3.08 m tall, and its (1,1) mode is 67.97 Hz rather than the 67.67 Hz the
nominal dimensions predict.

**The bar's peak is not its mean.** In `04`, six joules raise 20 mm of aluminium by 1.240 K on
average, but the peak reads 1.302 K after four seconds. Heat arriving on a plain channel has no
place, so `Bar1D` puts it in the first cell and conduction has not finished levelling it.

**A bar's panel is in kelvin and its readings are in celsius**, and both say which. The field
returns what the cells hold; celsius is a conversion a view chooses, and no view here makes it.
So `04` prints `|K| 294.4821` for a bar whose hottest cell is at 21.33 °C. It said `|C| 21.3321` while this
application sampled fields itself and applied the offset in the same expression as the label —
which meant nothing could disagree with anything. `FRICTION.md` 22 is that gap, still open.

**The bouncing ball ends at zero.** `07` prints the last frame, and by one second the ball has
stopped. The bounces are in the middle of the strip, not at its end.

## Three dimensions

The physics always had them: `NBody`, `ContactSystem` and `Fluid` all carry `DVec3`, and
flattening to a plane was the renderer's simplification rather than the simulation's. Bodies
are drawn in an axonometric projection now, sorted back to front, with radius growing toward
the viewer and colour mixed toward the plate for distance. Without all three the picture is
flat however true the coordinates are.

Not isometric, deliberately: a true isometric view puts the axes at 120 degrees and makes a
cube ambiguous, which is a bad way to read a periodic box. Rooms and bars stay as they are —
`Room` is two-dimensional by construction, and the crate says why: a third dimension costs a
factor of √3 in the stability limit as well as the obvious one in cells.
