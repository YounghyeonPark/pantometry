# Examples

What each one demonstrates and what it checks itself against. Every one runs in CI on every commit,
and every one is a closed-form check rather than a demonstration that something did not crash.

Add an output path to any of them and it draws the result:

```sh
cargo run --release --example lens_spots -- lens.svg
```

## A taste

```rust
use pantometry_optics::{Material, Spectrum, SurfaceFinish, fresnel_reflectance};
use pantometry_units::Length;

// Reflectance is not a setting. It follows from the refractive indices.
let bk7 = Material::from_catalog("N-BK7").unwrap();
let n = bk7.index(Length::nm(587.56));                 // 1.5168
let bare = fresnel_reflectance(1.0, n, 1.0);           // 0.0421 — the textbook 4%

// A coating can only scale that down, and it does so spectrally.
let coated = SurfaceFinish::broadband_ar()
    .reflectance_at(1.0, n, 1.0, Length::nm(550.0));
assert!(coated < bare / 10.0);

// A lamp has a temperature, and Planck decides what colour it is.
let tungsten = Spectrum::blackbody(3200.0);
assert!(tungsten.at(Length::nm(450.0)) < 0.45 * tungsten.at(Length::nm(650.0)));
```

The dimensions carry the coupling between domains too, and check it:

```rust
use pantometry_units::{Area, HeatCapacity, Irradiance, Length, Mass, Power, SpecificHeat, Temperature};

let absorbed: Power = Irradiance::mw_per_cm2(50.0)
    * (Length::mm(10.0) * Length::mm(10.0))   // an Area
    * 0.02;                                   // what SurfaceOptics::absorptance returns

let capacity: HeatCapacity = Mass::g(2.0) * SpecificHeat::j_per_kg_k(858.0);
let rise: Temperature = (absorbed * pantometry_units::Time::s(1.0)) / capacity;
// 0.58 mK, and every step of that chain landed on the dimension that names it.
```

## Examples

```
cargo run --example beam_hot_spot            # numbers, checked
cargo run --example beam_hot_spot out.svg    # and a picture
```

| | |
| --- | --- |
| `beam_hot_spot` | A 100 W laser on a mirror. Optics and heat meeting over a shared boundary, and the peak temperature a lumped model cannot see |
| `airy_pattern` | What a *perfect* lens does to a point: the Airy pattern, its encircled energy, and the MTF that follows from both |
| `detector_snr` | Why read noise is the first number on a datasheet, with the Poisson statistics sampled rather than asserted |
| `room_modes` | Why a small room booms at one note, and why that note is not a note. Four mode shapes and a corner trace |
| `melting` | A Lennard-Jones crystal melting, read off its own radial distribution rather than declared |
| `lens_spots` | An achromatic doublet *solved* from the Abbe numbers and then traced ray by ray in 3D — focal length, spherical aberration and colour, each against a formula the trace did not compute |
| `heat_in_three_dimensions` | A point of heat in a block of aluminium. The peak falls as `t^(-3/2)`, and **that exponent is the dimensionality** — a bar gives `-1/2`, a plate `-1`. Nothing 1D can produce it |
| `room_in_three_dimensions` | The same room as `room_modes`, with a ceiling. The floor-to-ceiling mode at 71 Hz that a floor plan does not have *at all*, and a mode count growing as `f³` rather than `f²` |
| `optical_bench` | **A 3D instrument, not a graph.** A doublet, a fold mirror turning the axis through 90°, three field angles — prescribed, traced, refocused, then *bent* until the spot falls inside the Airy disc. `optical_bench bench.html` gives a layout you rotate in a browser |
| `espresso_shot` | **A machine, and the inside of what it is doing.** An espresso basket from the pump to the cup: Darcy's law solved on the permeability the grind gives, the dissolution that rides on the flow, and a vertical cut through three baskets — even, channelled, and pulled into a cold portafilter — as the shot runs. Grind, temperature and pressure each swept on their own, against the exponent each is supposed to carry |
| `portafilter_flow` | **The machine, and the water going through it.** A shower screen, a basket, a body and a spout, with parcels of water leaving the screen and working down through the grounds — advected by the solved Darcy field at the *pore* velocity, and darkening as they pick up what they cross. Two baskets side by side, differing only in a loose ring at the wall. `flow.html` turns in a browser, `flow.gltf` opens in Blender, and `flow.json` opens in the native window with `pantometry view` |
| `busbar_rating` | **A design study rather than a demonstration.** A bolted busbar joint, from geometry to a production yield: the contact resistance solved as a field, the thermal path from a network, the electro-thermal fixed point, a rating by bisection, the margin to runaway, and 20 000 units against manufacturing tolerance. Every step against a closed form |

Two more are run by CI without being in the table, because they are checks rather than
showcases: `agents_quickstart`, the runnable form of [AGENTS.md](AGENTS.md), and
`readme_check`, which re-runs this file's own code so the snippets above cannot rot.

A fifteenth, `where_the_time_goes`, is a benchmark and is **not** run by CI. It measures rather than
asserts, and a timing threshold on a shared runner fails for reasons that have nothing to do with
the code. Run it by hand when a change should have made something faster: it is dependency-free,
takes best-of-five, and prints where a step actually spends itself.

Each one prints its numbers and **asserts** them, so CI runs all of them on every commit.
An example is a claim that the library works, which makes a quietly broken one worse than
no example at all — every value printed has been checked against a closed form or against
a calculation that did not go through the same code. Give a path and it also writes an SVG;
give none and it just checks. Nothing generated is committed — with one exception, and it is named so that this sentence stays true: `tools/presets` writes `presets.rs` and the chooser's tiles, which are build *inputs* rather than outputs and cannot be produced without a GPU.

Plotting has no dependency. SVG is text, so it is a `format!` and a file write — no
encoder, no fonts, and it opens by double-click. `examples/common/svg.rs` is about three
hundred and fifty lines and is the right size for this job; when it stops being, the answer is
a crate rather than a bigger version of it. That turned out to be a prediction rather than a
plan: `pantometry-world` had to write its own renderer, because this one lives under `examples/`
where no other crate can reach it — `FRICTION.md` finding 4.

`pantometry-view` is that crate now, and it does **not** close finding 4, which is worth being
precise about. It draws a `Frame` — a captured instant of a running simulation. An example plots
an encircled-energy curve or an MTF against spatial frequency, which is an arbitrary pair of
axes and not a frame of anything. The two want different interfaces, so the duplication is still
there and the finding stays open on its own terms rather than being closed by something adjacent
to it.

The examples also exist to keep the library honest in a way tests cannot. `ScalarField` was
written as the interface a visualiser would read a simulation through and then sat with no
implementor at all, which meant "is it the right interface" was a guess. There are two now —
a one-dimensional bar governed by diffusion and a two-dimensional room governed by a wave —
and the pair has said more about the interface than either could alone. See the next section.

An actual visualiser has since said something neither could: `ScalarField` is the right
*shape* and is unreachable through the thing you have. A renderer holds `&dyn Domain` and
there is no way to ask it for a `&dyn ScalarField`, so `pantometry-world` downcasts to concrete
types instead and knows every domain by name — which is exactly what the interface existed to
avoid. `FRICTION.md` finding 3.
