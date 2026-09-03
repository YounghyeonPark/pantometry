# The shape of pantometry

The goal is to reproduce physical law in three dimensions, and to do it in a structure that can
accept physics nobody has written yet without the parts already written having to change.

That second clause is the whole design. "All of physics" is unbounded, so no architecture can
contain it by enumeration — the only thing an architecture can do is stay open. Everything below
is in service of one property: **a new physics is a new crate, and nothing existing moves.**

This document is the map. It says what each layer owns, what is built, what is missing, and which
rules cannot be broken without losing the property above. It is not a roadmap and it does not
promise dates.

---

## The three layers

```
    ┌─────────────────────────────────────────────────────────────┐
    │  ANALYSIS      what a person or an agent asks of a run       │
    │                cameras · 3D · 2D · graphs · measurements     │
    │                                              `pantometry-view`   │
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  SCENE         where things are, and how they meet           │
    │                placement · interfaces · one clock            │
    │                                             `pantometry-scene`   │
    └───────────────────────────▲─────────────────────────────────┘
                                │  reads
    ┌───────────────────────────┴─────────────────────────────────┐
    │  PHYSICS       what evolves, and what it conserves           │
    │                the kernel, and one crate per physics         │
    │                         `pantometry-core` + eleven domain crates │
    └───────────────────────────▲─────────────────────────────────┘
                                │  fills
    ┌───────────────────────────┴─────────────────────────────────┐
    │  INPUT         what a person designed, made into cells       │
    │                a mesh · a rasterisation · what it lost       │
    │                                             `pantometry-shape`   │
    └─────────────────────────────────────────────────────────────┘

As of 0.9.0 the upper three are crates on crates.io. Before that the top two lived inside an
unpublished application, so a consumer could run a simulation and had no way to see it.
```

**The arrows point one way and that is load-bearing.** Analysis reads the scene; the scene reads
physics; physics reads neither. A domain that could see the scene could see another domain
through it, and the moment that is possible the crate split stops meaning anything.

### Layer 1 — physics

A domain owns state and advances it. It declares what it conserves, and it meets other domains
only on the `Exchange`, which carries **amounts** and not state.

The kernel — `pantometry-core` — knows conservation, integration, scheduling, sampling and boundaries.
It knows no physics. That is the rule that makes the goal reachable at all: light, heat, motion,
sound, matter and electricity have each arrived as a crate, and none of them required the kernel
to learn anything about them.

### Layer 2 — scene

Where things are, and how they meet. `pantometry-scene`.

It owns three things:

- **Placement.** A domain's state lives in its own coordinates; the scene says where those sit in
  the world. A grid gets a pose, a body set gets a pose, and a lumped model gets a *presentational*
  position — see below, because that distinction is the whole difficulty.
- **Interfaces.** Already in the kernel: `Interface` and `Flux` carry an amount *and where on a
  surface it crossed*, audited face by face rather than on the total.
- **The clock.** One `Simulation`, one schedule, one audit across everything placed in it.

`capture` turns all of that into a `Frame` at an instant, by asking each domain what it *offers*
— `as_field` for a continuum, `as_bodies` for a countable set, `readings` for scalars. It names
no domain, and `knows_no_physics.rs` demonstrates rather than asserts that: it defines a physics
inside the test file and captures it whole.

### Layer 3 — analysis

What a person asks of a run. `pantometry-view`: a filmstrip as SVG, a self-contained HTML report,
a CSV of every domain's scalars, and the frames as JSON.

The rule here is proven and should be extended rather than replaced: **views dispatch on the
shape of the data, not on the name of the domain.** Scalars over time become a chart, a 1D field
a profile, a 2D field a heatmap, points a 3D scene. A new domain gets a correct picture without
this layer learning it exists.

Its tests are driven by frames written out by hand rather than by a simulation, which is the only
way to check that rule at all: a test that ran a real scene could not tell *a heatmap because the
data is a 2D grid* apart from *a heatmap because that domain was a room*.

One more rule holds throughout, and it is the easiest to break by accident: **the scale is fixed
across a run**. A picture that renormalises per frame makes a decay look like a steady state, and
per-frame normalisation is what you get if you do not think about it.

### Layer 0 — input

`pantometry-shape`, and it is under the physics rather than beside it. Every domain
here takes a **structured grid**; every file a person designs is a **surface**. This is the bridge,
and it is the answer to a question the other three layers cannot be asked: *where does the geometry
come from?* Until now, from a closure written by hand, which is a fine answer for a test and no
answer at all for someone who has a part.

It depends on `pantometry-units` and nothing else in the workspace, and **no domain depends on it**. It
produces a predicate, `|i, j, k| voxels.contains(i, j, k)`; a domain's `fill` consumes one. That is
the entire coupling, and it is why adding geometry cost the domains nothing — `Solid3D::fill`,
`Block::fill` and `Waves::fill` already had exactly that signature for other reasons.

**The layer exists because the cell size is a physics decision disguised as a performance one.**
Picking `cell_mm` picks three things at once: the stability limit and step count as `dx²`, the
discretisation error, and — the one with no symptom — *how much of the shape survives*. A 0.5 mm rib
voxelised at 2 mm is not a thin rib, it is gone, and the simulation runs perfectly well without it,
conserves energy, audits clean and answers a question about a different object. So the crate's
output is not only cells but a `Loss`: volume error, the fraction of the volume in boundary cells,
runs one or two cells thick, and rows the rasteriser could not decide.

Two things it deliberately does not do. It does not **choose** the cell size — that would be the
library guessing, and refusing to guess is the property the rest of this document is about. And it
does not repair, refine or simplify a mesh; it reads what was exported, measures it, and says
whether the answer is trustworthy.

What is still missing is named honestly in *What is not built* below: a grid has no **void**, so the
cells a part does not occupy are some other substance rather than nothing, and two parts voxelised
separately have no way to touch.

---

## The crates

| Crate | |
| --- | --- |
| `pantometry-units` | Dimensional analysis. SI quantities and vectors whose dimension lives in the type, so `Length + Time` does not compile |
| `pantometry-core` | The kernel: conservation audits, fixed-step integrators, fields, shared boundaries, multi-domain scheduling, deterministic sampling, parallel ensembles for Monte Carlo, closed-form rigid motion |
| `pantometry-optics` | Light: spectral radiometry, surface optics, dispersion, ray geometry, diffraction |
| `pantometry-thermal` | Heat: lumped masses, conduction in one dimension and in three, networks of bodies joined by conductances, radiative and convective loss |
| `pantometry-mechanics` | Motion under force: N-body, Barnes-Hut, penalty contact, rigid rotation |
| `pantometry-acoustic` | Sound: the wave equation on a staggered grid in one, two and three dimensions, impedance boundaries |
| `pantometry-molecular` | Matter atom by atom: Lennard-Jones fluids in periodic boxes, cell lists, a Langevin bath, radial distributions |
| `pantometry-electrical` | Electricity: resistive dissipation into the heat channel, conductors whose resistance moves with temperature, and a **field** formulation where `I²R` is solved out of a shape rather than stated |
| `pantometry-fluid` | Incompressible Navier–Stokes by projection on a staggered grid. Poiseuille and Couette come out exact against the *discrete* profile, whose gap to the continuum is a closed form and not a tolerance; Taylor–Green checks the nonlinear term against `e^{−2νk²t}`; and momentum in a periodic box drifts by `1e-15` |
| `pantometry-em` | Maxwell's equations on a Yee grid, where `∇·B = 0` is an **identity of the update** rather than a tolerance: inject a divergence and 500 steps later it is unchanged to nine digits. Cavity resonances at second order, the leapfrog's energy swing against its own closed form, an absorbing boundary that leaves 0.149% of a pulse behind where a conductor leaves all of it, and a waveguide whose dispersion relation comes out of one march |
| `pantometry-elastic` | What a shape does under load: `∇·σ = 0` solved on trilinear elements, so a stiffness is a property of a geometry. Four moduli come out exactly — `E`, the constrained `M`, the bulk `K` and the shear `G` — and Clapeyron's `2U = Σf·u` says the discretisation is self-consistent |
| `pantometry-porous` | Flow through a packed bed: Darcy's law solved as a field, the heat the liquid carries, and the dissolution that rides on both. An espresso puck, and also a filter, a catalyst bed and an aquifer |
| `pantometry-quantum` | A wavefunction in a well: the time-dependent Schrödinger equation, marched with the same staggered leapfrog family the acoustic domain uses — real part on integer steps, imaginary on halves — so **probability is conserved as an identity of the update**, the way `∇·B` is on the Yee grid. Eigenvalues against the discrete operator's own closed form, Gaussian spreading and Ehrenfest's theorem at second order |
| `pantometry-shape` | Designed geometry as input: an STL read and measured, and rasterised into the cells a domain fills — with a report of what the cells could **not** hold, because a rib finer than the grid does not fail, it disappears. Depends on `pantometry-units` and nothing else |
| `pantometry-scene` | Where things are and what a run looks like: placement, capture, and the shapes a view can draw. Names no domain |
| `pantometry-view` | Drawing that: a filmstrip, a self-contained HTML report, CSV, JSON, and **glTF** so Blender, three.js and USD tools can open a result — as shaded surfaces, not a point cloud. The view is chosen by the shape of the data, never by the name of a domain. No dependencies |
| `pantometry` | A facade over the other sixteen, and where the cross-domain integration tests live — including the three that hold two domains against each other: a Yee grid against Fresnel's algebra, a field's decay in a conductor against a lumped resistance that has no frequency in it, and a diffraction pattern against the scalar theory it converges on |
| `bindings/python` | Python bindings, in their own cargo workspace and on PyPI as `pantometry`. SI floats at the boundary and the conservation audit as a catchable exception — the dimensional types are compile-time and cannot cross |
| `app/pantometry-gpu` | `Solid3D`'s stencil as a compute shader — **33–67× on a 64³ grid** and a wash at 16³, measured one grid per process by a test that prints the adapter it ran on. Single precision against the domain's double, so the CPU is the reference and the difference is measured. A scene says `"device": "gpu"` and the binary honours it |
| `app/` | Everything a person runs, as one binary: `pantometry run | check | verify | view | edit`. Its own workspace, because a GPU stack is 86 external crates and a GUI shell 371 against the library's 12. `viewer-core` inside it depends on the run **file**, not on `pantometry`, so the wire format being sufficient is demonstrated rather than claimed |
| `pantometry-world` | The first consumer, and not published. Worlds described as data: built, coupled over the bus, run and drawn, with thirty scenes across all eleven domains that CI runs. It exists to use the SDK from outside and write down where that is awkward |

The last three are the workspace's answer to the same question from three sides: what a
simulation *is* (`pantometry-scene`), what a picture of one *is* (`pantometry-view`), and what it feels
like to use both from outside (`pantometry-world`). The first two are libraries because a consumer
who can state a simulation should not have to write a plotting stack to see it.

```text
pantometry-units       no dependencies but glam and serde
pantometry-core        depends on units                     ── the kernel
pantometry-optics      depends on core     ─┐
pantometry-thermal     depends on core      │
pantometry-mechanics   depends on core      ├─  one crate per physics, and
pantometry-acoustic    depends on core      │   none of them knows another
pantometry-molecular   depends on core      │
pantometry-electrical  depends on core      │
pantometry-elastic     depends on core      │
pantometry-em          depends on core      │
pantometry-fluid       depends on core      │
pantometry-porous      depends on core      │
pantometry-quantum     depends on core     ─┘
pantometry-shape       depends on units only                ── designed geometry in
pantometry-scene       depends on core                      ── where things are
pantometry-view        depends on scene                     ── how to draw that
pantometry             depends on all of them
pantometry-world       depends on the facade, and nothing depends on it
```

**The kernel must never depend on a domain.** If a new physics needs the kernel
changed, the kernel was wrong — that rule is what makes "add sound, add fluids" a
matter of writing a crate rather than editing this one.

**And the layers above must not know one either**, which is the same rule from the other side.
`pantometry-scene` asks each domain what it *offers* — a field, a set of bodies, some readings — and
`pantometry-view` dispatches on the shape of what came back. Neither names a domain, and neither
merely claims that: `pantometry-scene`'s test defines a physics inside the test file and captures it
whole, and `pantometry-view`'s tests are driven by frames written out by hand, because a test that
ran a real scene could not tell *a heatmap because the data is a 2D grid* apart from *a heatmap
because that domain was a room*.

Six domains are now the proof rather than an assertion. Optics publishes absorbed
light as heat and thermal consumes it; mechanics publishes a dashpot's dissipation on
the same channel and the same thermal domain consumes that too, with nothing changed
on either side; acoustics publishes what an absorbing duct end radiates onto the same
channel again; molecular dynamics arrived without asking for anything; and electricity
arrived last, publishing `I²R` onto that same channel. None of the six names another and
none of them needed the kernel changed.

Electricity is the one that closes a circle. Every other producer of heat here answers a
question about *something else* that happens to warm a thing — light landing on a mirror,
a dashpot damping a bounce. A winding is the case where getting hot is the entire subject,
and until that crate existed the workspace's own examples stood a stated number of watts in
its place. A stated number cannot be wrong, which is another way of saying it is not a model.

One thing did need the kernel changed, and it is worth being precise about why that is
not a violation of the rule. The rule is that a *domain* must never force a kernel edit.
What forced this one was the coupling mechanism itself being under-specified from the
start: a bus carrying one number per channel could not say *where* on a surface a
quantity crossed, so no domain could ask for that and none of them ever did. Adding
`Interface` and `Flux` did not teach the kernel any physics — a discretised boundary is
not optics or heat — and no existing domain had to change to keep working.

They also brought the audit three conserved quantities instead of one, and the three
hold to wildly different tolerances for structural reasons rather than through
differences in effort:

| Quantity | Holds to | Because |
| --- | --- | --- |
| Linear momentum (`NBody`) | 1e-13 | Exact by construction: equal and opposite forces cancel bit for bit |
| Linear momentum (`TreeNBody`) | θ-dependent | Each body sees its own approximation of the rest, so nothing cancels |
| Angular momentum (`RigidBody`) | 1e-9 | Nothing cancels here either; it is only as good as RK4 plus a quaternion renormalisation |
| Energy across a coupling | 1e-9 | Both sides are closed-form evaluations, nothing is integrated |
| Energy through a contact | 2e-2 | A penalty contact is a non-smooth potential, so the symplectic bound does not apply |

Auditing a vector component by component also makes the *smallest* component the
binding constraint, since the absolute error is set by the whole vector while the scale
it is judged against is only that component. That is a property of the kernel's
per-quantity audit, and it is written down where it will be met.

Storage is SI base units everywhere: metres, kilograms, seconds, kelvin. Millimetres
and nanometres are entry and exit forms, and the unit-bearing constructors are the
only place a factor of a thousand can hide.

## What is not here

No scene *graph* — no parent-child transform hierarchy, no culling, no traversal order. What
there is, is flat: `pantometry-scene` gives each domain a `Pose` in world coordinates and captures
what they hold. A hierarchy is what you want when placements are relative and animated, and
nothing here has needed one.

The renderer is deliberately modest. `pantometry-view` draws a filmstrip, a heatmap, a profile, a
point scene — depth-sorted back to front, which is painter's algorithm on a 2D canvas, with no
depth buffer and no shading — and a **raycast** of a 3D field, composited front to back and
rotatable, beside a montage of every slice. The render shows shape and cannot be read for values;
the montage is the reverse, and a volume gets both rather than a choice between them. The SVG has one fixed projection; the HTML report can be dragged to
rotate and scrolled to zoom, and that is the whole camera model. It is enough to see whether a
simulation did what you expected, and it is not a visualisation package. The JSON export exists
for when it is not enough.

The **export** is where the pictures in somebody else's renderer come from, and it stopped being a
point cloud. A three-dimensional field is the surface of its present cells — one quad per face
whose neighbour is absent, so 89% of a solid block's faces are culled as unseeable and a void
inside it produces a real interior surface. Bodies are spheres. Both carry normals, so they take
light and cast shadows, which a point cloud cannot. Rendered from Blender straight off the export,
scene 23 is a hot part and a cooled lid with a real gap between them, and scene 16's room shows its
standing wave's nodal planes as dark bands across a solid.

**USD** is the other half of that, and it does what glTF cannot: `out.usda` is the *whole run*.
USD has time samples on any attribute, so the topology is written once and the colours per frame,
and usdview's timeline scrubs the physics rather than a camera move. Verified by rendering it
through Blender's OpenUSD: scene 23's part is cream at frame 0 and green at frame 20 while its lid
warms from navy, both on one scale across the run so the two frames are comparable.

It carries the numbers too. A domain with no geometry at all — a heater, a lamp, a winding — is
most of what a scene here contains, and its scalars go out as time-sampled custom attributes under
`pantometry:`, which usdview shows in its property panel and scrubs with the timeline. There is no
USD schema for a `Ledger` or a `Violation` and this invents none; a custom attribute is a number
with a name, which is what a reading is.

Still no dependency. `.usda` is USD's text serialisation and this writes it by hand, for the same
reason glTF and SVG are written by hand here — and `usdcat -o out.usdc out.usda` is one command if
a pipeline wants the binary crate.

Three things came out of doing it. The colours were being written into glTF's `COLOR_0` as sRGB
where the specification says linear, so every export this workspace has ever produced was decoded
about 2.3x too bright in the midtones. And giving each sample a full cell made an object one whole
cell larger than the extent it was sampled over — 2x on a two-sample axis — because `capture`
samples corner to corner and an end node owns *half* a cell, which is the third time that
arithmetic has been the bug here. And the USD wrote a primvar's interpolation as a separate
property where USD reads it as attribute **metadata**, so a file holding a colour per vertex was
drawn with the first one over the whole prim — found by rendering it, where a hot part and a cooled
lid came out the same colour, and invisible to any check on the text.

Modest is not the same as unreadable, and the report was the second for a while. Axes are in
metres rather than in cells, hovering reads back the sample under the cursor, the scalar chart has
an axis per unit instead of normalising every series to itself, and the colour scale is
constructed in CIE LCh with lightness linear in the value — properties `pantometry-view::ramp`
pins with tests rather than a designer's judgement. The viewer that draws all of it is executed by
`tools/report-check`, which is new for the ordinary reason: it was four hundred lines of
JavaScript that nothing had ever run.

The JSON scene *format* is still `pantometry-world`'s and not the library's, which is why that crate
is unpublished. A file format is a compatibility promise, and this one is not ready to make one:
it renamed a field once already and nothing failed, because serde discards unknown keys by
default — which is how `deny_unknown_fields` came to be on every type in it.

Rigid bodies are spheres where they collide. `RigidBody` rotates with an applied torque
and an arbitrary inertia tensor, but `Sphere` and `Rolling` are the only things that
touch anything, and a sphere is chosen because its inertia is the same about every axis
— so a contact does not have to carry the orientation through. Boxes hitting boxes is a
different module.

Acoustics is linear. A tube, a room and now a hall — `Hall` is the wave equation in three
dimensions, which is what the vertical and oblique modes need, and a floor plan does not have
them at all rather than having them inaccurately. Still no scattering geometry and no
nonlinearity, so nothing here shocks up or distorts, and no absorbing wall model beyond an
impedance boundary.

Molecular dynamics is monatomic, so the crate's name is aspirational by one step: no bonds,
angles or torsions, and therefore no molecules. And no electrostatics, which for a charged
system is the hard part rather than a missing feature — Coulomb falls off as `1/r` and cannot
be cut off at all, so it needs Ewald summation or a particle-mesh method, and that is a larger
piece of work than everything currently in the crate. No constraints, no barostat, no
free-energy machinery.

Both acoustic domains had a boundary defect until recently, and how it was found is worth
more than the fix. See [EVIDENCE.md](EVIDENCE.md#a-boundary-defect-and-the-thing-that-found-it).

**There is no fluid domain, and that is deliberate.** Sound *is* the fluid domain here:
it is what a fluid does when the variations are small enough to linearise, which is
exactly the regime where every answer has a closed form to check against. Full
Navier-Stokes has none, and a solver that could not be validated against anything would
be decoration. See the note below on turbulence.

Gravity comes both ways and the choice is a real one. `NBody` sums every pair: exact,
momentum conserved to the last bit, `O(n²)`, and awkward to parallelise precisely
because the `i < j` pairing that makes it exact has two threads writing to one body.
`TreeNBody` is Barnes-Hut: `O(n log n)`, embarrassingly parallel, and it **gives up
exact momentum** — each body sees its own approximation of the rest, so their mutual
forces no longer cancel. The drift is a knob rather than a defect, closing with the
opening angle and vanishing at `θ = 0`, and the audit tolerance has to be the one the
angle earns. The expansion carries the quadrupole as well as the monopole, which buys
back most of that accuracy at close angles — a factor of 6.5 at `θ = 0.3` — but nothing
past it, so `θ` above about 1 is still asking a centre of mass to stand in for a group
that is not far enough away.

Fields propagate between planes by the angular spectrum, but only through free space:
there is nothing to put in the beam's way except an aperture, and the grid's reach is
`NΔ²/λ`, past which the propagator refuses rather than aliasing.

Partial coherence is a transfer function rather than a simulation. The two exact limits
are there, and so is the coherence a source has at a distance, but computing an
arbitrary object's partially coherent image needs the transmission cross-coefficients —
a four-dimensional integral, and not here.

The tree stops at the quadrupole. Higher multipoles would buy another order in the
opening angle, and a proper Fast Multipole Method would change the complexity rather
than the constant.

Meshes and grids are no longer excluded. They were, on the grounds that adding them
before a second consumer would be guessing at an interface; finite elements and
finite differences need them, so that decision has been reversed deliberately rather
than drifted away from.

An `Interface` is one-dimensional: a boundary is a sequence of faces, and `Flux::resample`
remaps between two of them by overlap in cumulative area. That is enough for a mirror
face, a bar's side, a row of pixels, and any boundary whose faces have a natural order —
and it is not enough for a triangulated surface or an arbitrary mesh-to-mesh projection,
where the overlaps are not an interval intersection and there is no order to walk. The
conservation argument generalises; the implementation does not, and the doc comment says
so where a caller will meet it rather than only here.

An `Interface` also carries areas and an order and nothing else — no coordinates, no
normals, no connectivity. A domain that needs to know where a face *is* in space still
has nowhere to put that, so a beam profile has to be handed over in the boundary's own
coordinate rather than computed from geometry. That is the next thing this layer is short
of, and it is deliberately not guessed at before something needs it.

And some things stay out for reasons that will not change: general relativity and
quantum field theory are research subjects rather than simulation targets, turbulent
DNS at world scale is a question about supercomputer budgets rather than about API
design, and no single `f64` state vector spans the fifteen decades from a nucleus to
a galaxy — a simulation has to declare which regime it is in.

## Where the work is: 3D is not the default yet

The physics layer is eleven crates deep and dimensionally uneven. This is the honest state.

| crate | space it lives in | for the goal |
| --- | --- | --- |
| `pantometry-mechanics` | **3D** — bodies, contacts, rigid rotation | done |
| `pantometry-molecular` | **3D** — atoms in a periodic box | done |
| `pantometry-optics` | **3D rays**; no volumetric field | rays done, fields missing |
| `pantometry-acoustic` | **3D** `Hall`; **2D** `Room`; **1D** `Tube` | done |
| `pantometry-thermal` | **3D** `Solid3D`; **1D** `Bar1D`; `ThermalNetwork` is a graph with no space | conduction done |
| `pantometry-electrical` | **3D** `Conductor`; `Winding` is a lumped `I²R` | done |
| `pantometry-porous` | **3D** — Darcy flow, advected heat and dissolution in a packed bed | done |
| `pantometry-elastic` | **3D** — linear elasticity on trilinear elements, statics and waves | done |
| `pantometry-em` | **3D** — Maxwell on a Yee grid, conducting walls | done |
| `pantometry-fluid` | **3D** — incompressible Navier–Stokes by projection | done |
| `pantometry-quantum` | **1D** `Well` — a wavefunction between walls | arrived one-dimensional, as every wave here did; higher dimensions are cells and cost, not new physics |

Two observations follow, and they point in opposite directions.

**The reductions were deliberate.** `Room` is 2D because a third dimension costs √3 in the
stability limit on top of the obvious factor in cells, and it said so before anyone asked.
`ThermalNetwork` refuses to have positions because a conductance is not a distance. These are not
unfinished work. They are the cheap models that answer engineering questions, and a platform that
deleted them would be worse at its job.

**And they are not the goal.** A 3D platform needs 3D conduction, a 3D wave and a field
formulation of current, and those are three new crates or three additions, none of which requires
the kernel to change.

Both are true. The resolution is in the next section.

---

## Lumped models are reductions, not exceptions

A lumped thermal mass is a 3D conduction problem with the interior collapsed to one temperature.
That is valid when the Biot number is below about 0.1, and `LumpedMass::biot_number` exists so a
caller can find out rather than assume.

So in a platform whose goal is 3D, a lumped model is **a reduction with a stated validity**, and
it earns its place by being fast where the reduction holds. The alternative is not a purer
platform; it is an unusable one. Heat conduction through a motor housing, resolved as particles,
is picosecond steps against a two-thousand-second time constant — about 10¹⁵ steps for one answer
a graph of four nodes gives immediately.

What the platform owes such a model is not deletion but **honesty about what it is**: which 3D
problem it reduces, under what condition, and what it cannot show. `ThermalNetwork` cannot show a
hot spot, and its documentation says so.

---

## Placement: built, in two halves that must not touch

The kernel has `ScalarField` and `VectorField` — functions of position — `Interface` and `Flux`
for discretised boundaries, and now `Pose`. Until `Pose` a domain's coordinates *were* world
coordinates, implicitly, and two grids could not be placed against each other at all.

`Pose` is a rigid motion and deliberately nothing more: a rotation and a translation, no scale,
no shear, no projection. An isometry preserves every distance and angle exactly, and that is the
only class of placement a physics can be moved by without its physics changing — a conservation
law stated over a sheared volume is a different law, and a scaled metre is not a metre.

Placement has two uses, and they must not share a type:

- **Physical placement** — changes what the physics computes. Two solids in contact, a lens at a
  distance, a grid rotated against another. This is `Pose`, and it is **built**. The scene layer
  assigns it; a domain reads only its own coordinates.
- **Presentational placement** — a position given to something that has none, purely so a viewer
  can draw it. A thermal network node on a diagram. This is `Placement::marker`, and it is
  **built**, in `pantometry-scene`: above the kernel, above every domain, where no physics can reach
  it.

Keeping them apart is structural rather than a naming convention. They are separated by *which
crate they live in* — if they shared a type, someone would eventually feed a drawing coordinate
into a conductance and nothing would fail loudly.

`Placement` also carries an `Extent`, and that third field is the one nobody predicted. A
`ScalarField` is a function of position and does not stop anywhere; a field that knew its own
bounds would be a mesh. So the region to sample has to come from above, and the scene is where
the size was written down in the first place.

**And a third half, which nobody predicted either: a placement has to leave.** Two halves that do
not touch is a statement about the *inside* of a run; the outside is **six** readers — the glTF
and USD writers, the standalone viewer's path projection, the editor's shaded viewport and its flat
painter, and the HTML report — and `capture` was giving them half of each `Placement`. A field was sampled in the domain's own coordinates and the
pose dropped; a body had it multiplied into its position. One file, two conventions, and nothing
in it saying which a panel used. Two blocks half a metre apart came out of glTF at identical
coordinates.

`Panel::place` and run **format 2** are the fix: every panel is in its own frame and says where
that frame is. Baking is not wrong for a *vertex list* — a rigid motion loses nothing there — but
it is for a **box**, where the axis-aligned box around a rotated cell is bigger than the cell. That
is the asymmetry the whole thing turns on, and it is why a designed mesh's triangles may be placed
where a field's extent may not.

The glTF writer was fixed in the same commit; the other five took **three more**, found one or
two at a time, because **under the identity a reader that drops a placement gives exactly the right
answer** — and every shipped scene was the identity. Scene 30 is the one that is not, and it exists
for that reason alone.

And the last thing the episode taught is why there were six. Seven *places* were applying the
placement — three in the editor's shaded pass, three in its flat painter, one in the standalone
viewer — each in its own four lines, so a seventh careful reader was never going to be the last
one. `Panel::placed_positions`, `placed_vertices` and `placed_corners` are the one place now, and
every painter goes through them. The writers keep their own, because a placement leaves glTF and
USD as a **node transform** rather than as moved numbers.

---

## Rules that hold the whole thing up

These are not style. Each one is what makes some part of the goal reachable.

1. **The kernel must never depend on a domain.** Without this, "add a physics" means "edit the
   kernel", and the goal is a rewrite each time.
2. **No domain may depend on another.** They meet on the bus. Eleven domains have now been added
   without this breaking, which is the evidence that the split is real.
3. **The arrows point one way.** Analysis → scene → physics. A domain that can see the scene can
   see another domain through it. This is enforced by cargo rather than by discipline now that
   each layer is a crate: `pantometry-scene` does not appear in any domain's manifest, so a domain
   reaching upward does not compile.
4. **Conservation is audited, not assumed.** The audit is what makes an unfamiliar coupling
   trustworthy, which matters more as the number of domains grows, not less.
5. **Results are bit-for-bit across platforms.** A simulation that gives a different answer on a
   different machine is not a measurement.
6. **A number in prose is a number under test.** This repository has shipped stale figures
   repeatedly; the ones that stopped recurring are the ones a test now checks.

Rule 4 had a known limit and it is closed. The tolerance used to be one number for the whole
simulation; it is now **one per quantity**, and a domain can additionally claim
`books_balance` and be checked **on its own scale** rather than against the sum of every ledger.
Both halves matter and neither substitutes for the other: the first separates schemes carrying
different quantities, the second separates domains carrying the same one.

---

## What is missing, in the order it blocks things

1. ~~**Placement in the kernel.**~~ Done: `Pose` is a rigid motion — rotation and translation,
   no scale, no shear — so two things can be positioned relative to each other. The
   *presentational* half is still missing on purpose and waits for the scene crate, where the
   physics cannot reach it.
2. ~~**Scene and analysis as libraries.**~~ Done. Both were inside `pantometry-world`, which is
   `publish = false`, so a consumer who could state a simulation and run it could reach neither
   the shape of the answer nor any view of it.

   `pantometry-scene` is layer 2: `Placement`, `Extent`, `Frame`, `Panel`, `capture`,
   `settle_framing`. `pantometry-view` is layer 3: a filmstrip, a self-contained HTML report, CSV and
   JSON, with the view chosen by the shape of the data.

   **Neither names a domain, and both prove it by construction.** `pantometry-scene`'s test defines a
   physics inside the test file and captures it whole; `pantometry-view`'s tests are driven by frames
   written out by hand, which is the only way to tell "a heatmap because the data is a 2D grid"
   apart from "a heatmap because that domain was a room".

   Getting there needed **five** things the layer had been doing by matching on domain types:
   `Domain::readings`, `Domain::as_bodies`, `Simulation::domains`, `ScalarField::unit`, and
   `Placement::extent` for the region a field occupies. Each was invisible while one crate did
   everything, and each is a thing the *kernel* was missing rather than an invention of the
   layers above it.

   What is left in `pantometry-world` is what an application actually is: a file format, the domain
   types that format names, and one place saying how far each field extends.
3. ~~**3D field domains.**~~ Done. `Solid3D` is conduction through a block — a seven-point
   stencil, insulated faces, `dx²/6α`, checked against the exact eigenvalue of its own discrete
   operator. `Hall` is the wave equation with a ceiling — a staggered grid, rigid surfaces,
   `dx/(c√3)`, checked against the rigid-wall mode frequencies and a second-order convergence
   rate measured across three doublings.

   Nine of the eleven domains are three-dimensional now. What is left is `pantometry-optics`, whose rays
   are already 3D and whose *fields* are not — and gap 4 below, which closed `pantometry-electrical`,
   is why this sentence used to say eight.

   Building the first one **found a gap in the layer above it**, which is what a first
   three-dimensional anything is for. The second one found nothing, which is the evidence that
   the first fix was the right shape: `Hall` needed no change anywhere outside its own crate and
   the scene format. `Extent::samples` was a pair and the sampler built its
   position as `(u, v, 0)`, so a solid would have been captured as its `z = 0` face — silently,
   because a slice of a block is a perfectly plausible picture of a block. Nothing in
   `pantometry-scene` could have noticed on its own; every field it had ever been handed was flat.
   `samples` is a triple now, `PanelData::Field` carries `nz`, and the type system made all three
   view sites decide what to do about it.
4. ~~**A field formulation of electricity.**~~ Done. `Conductor` solves `∇·(σ∇φ) = 0` by
   conjugate gradients and reads `J = −σ∇φ` off it, so a resistance is a property of a shape.
   `ρL/A` comes out exactly for a uniform block; a notch comes out as whatever the notch gives,
   which is the point — spreading resistance has no closed form for an arbitrary geometry.

   It is the first **elliptic** domain here, and the first whose failure mode is a solver rather
   than a stability limit. An iterative solve stopped at its iteration cap returns a field that is
   smooth, bounded and shaped exactly like an answer, so `step` refuses one that did not converge
   and the residual is a *reading* rather than an internal number.
5. **Per-quantity tolerances.** Done. `Tolerances` and `Simulation::conservation_tolerance_for`
   give each conserved quantity its own number, with a default for the rest. A Barnes-Hut tree
   gives up exact momentum by construction while energy in a rigid room is exact to `1e-15`, and
   under one number either the momentum check refuses a correct run or the energy check stops
   being able to see anything. Both failures are demonstrated in
   `per_quantity_tolerances.rs` rather than asserted.

   **The other half is done too**, and it turned out easier than this document predicted.
   `Domain::books_balance` is an opt-in claim that a domain's ledger changes by exactly what it
   took from the bus minus what it published, and a domain that makes it is checked **on its own
   scale** every step. A domain holding a microjoule beside one holding a kilojoule can lose a
   fifth of itself without moving the sum by more than `2e-10`; on its own scale a fifth is a
   fifth.

   This document predicted the blocker would be ledgers that are honest approximations, and named
   `Room::startup_adjustment` as the example. **That was wrong.** `Room::energy` reports the
   released state's datum plus an offset chosen so the two agree, so its books balance exactly
   from the first step — the `O(h²)` correction is *inside* the number it reports and does not
   move it. Every domain in the workspace except `LumpedMass` takes the claim and passes.

   `LumpedMass` correctly declines, and it is the reason the check is opt-in rather than
   automatic: it loses heat to an environment that is not on the bus, so its ledger does not
   balance against bus traffic alone. That is a boundary being modelled, not a leak, and a check
   that accused it would be the wrong check.
6. **A renderer with depth.** Partly answered, and the answer turned out not to be a depth
   buffer. Content came first, as this list said it should: three domains produce volumes now, and
   what a volume wants is not occlusion but **integration along a ray**. `pantometry-view` raycasts a
   3D field — trilinear sampling, front-to-back compositing, rotatable — beside the slice montage,
   because a render shows shape and cannot be read for values while a montage is the reverse.

   **And now something has surfaces.** `pantometry_view::mesh::isosurface` draws where a field
   reaches a value — marching tetrahedra, so there is no 256-entry table to be wrong and no
   ambiguous case to leave a hole — and the editor's viewport, which has had a depth buffer since
   it was written, is where it is drawn. That closes the condition this entry set for itself: the
   depth buffer was waiting for something with a surface, and an isosurface is one.

   **And the other half is closed too: a designed mesh reaches the screen.**
   `pantometry_view::mesh::mesh_surface` takes triangles and returns a `Surface`, and the editor
   draws every `parts` entry of a scene where its voxels will be — so an assembly is visible while
   somebody is still authoring it, before there is a run and therefore before there is any field
   to draw a surface of. Until now the only picture of a designed part was the staircase it
   rasterised into, which is a picture of the grid as much as of the object.

   The layering held without an exception. `pantometry-view` still does not depend on
   `pantometry-shape` — `mesh_surface` takes the same plain arrays the other three producers take,
   and `editor-core`, which links both, does the conversion and applies the pose. It is the *same*
   pose `Voxels::onto` rasterises against, and that is not a detail: a mesh drawn a millimetre
   from its own voxels would be two pictures of one part that disagree, which this module's own
   header calls worse than either being wrong.

   Bodies are still points and a field is still composited rather than occluded, so the *flat*
   renderers in `pantometry-view` need no depth. That is still true and still fine.

   **The last thing this entry wanted is done, and not the way it asked for.** It said a designed
   mesh had to travel in a `Run` — a fourth `PanelData` beside fields, bodies and paths — and
   called that the next decision. Sizing it turned the question round: **a run is the simulation's
   output and an STL is its input**, and copying an input into an output is what makes the other
   answers bad. Repeated per frame it is unbounded — a bracket's 24 triangles are free and a
   million over a hundred frames are not. A static section is a new *kind* of thing in a format
   that has never had one. A path resolved on read makes a run depend on a filesystem, whose
   failure is an empty picture rather than an error.

   So `mesh::Drawing` is an argument to the writer instead. The **three** writers that draw
   geometry take one — `gltf_with`, `usda_with`, `html_with` — and the editor's viewport had the
   design before any of them, from the scene rather than from a run. The filmstrip and the CSV
   draw no geometry and are unchanged. The wire format did not change either, and the
   `deny_unknown_fields` reader this entry was worried about was never asked to.

   What the pair is *for* is the rasterisation loss. `Rasterised` has reported it as a volume
   error since designed parts existed, and a number in a terminal is not what a designer looks at.
   Scene 29 through OpenUSD: the cells reach `0.05096` where the design ends at `0.05`, and the
   overshoot is the grid.

   It is drawn **uncoloured** in all four, from one constant. The solver ran on cells; tinting a
   smooth surface from the field would claim a resolution the computation never had, which is the
   per-frame colour scale's mistake wearing different clothes.

   **What did change the run format was something else entirely**, and finding it is why this
   entry took five commits rather than one. A run had no coordinate frame: `capture` held each
   domain's `Placement` and used half of it, sampling a field in the domain's own coordinates
   while multiplying the pose into a body's. Two blocks half a metre apart exported to glTF at
   identical coordinates. `Panel::place` and run **format 2** are that fix, and then three more
   readers of it had to be found one at a time — USD's writer, the editor's viewport, the HTML
   report — because under the identity a reader that drops a placement gives exactly the right
   answer and every shipped scene was the identity. Scene 30 is the one that is not.
7. **Geometry from a designed file.** Done as far as one object goes: `pantometry-shape` reads an STL,
   measures it, and rasterises it into the predicate a domain's `fill` already took. Nothing in the
   physics changed.

   The three that blocked a real assembly are done, and each closed differently:

   - ~~**A grid has no void.**~~ Done. `Solid3D::empty` marks cells as nothing: no capacity, no
     conduction across any face they touch, no share of what arrives on the bus, no vote in any
     average, and **no temperature** — `temperature_at` returns `NaN` there, because a zero or an
     ambient is a value somebody would plot and believe.

     It cost almost no new arithmetic, which is worth recording: the face mean is already the
     harmonic one and already guarded at zero, so a face touching a zero conductivity carries
     zero without the face knowing what void is. The one special case is that a cell with no
     capacity has no `dx/C` either.

     **What it did cost was every reader downstream, and that was not noticed for two commits.**
     `temperature_at` returned `NaN` from the first day and `ScalarField::at` — which is what
     every exporter actually calls — read the raw cell array, where an emptied cell still holds
     whatever it held when it was emptied. So a clearance left this workspace as a piece of the
     block sitting at its start temperature forever. Measured on scene `23`: the glTF carried
     **252** points for a grid with **120** solid cells. The lesson generalises past void — an
     invariant that only the domain's own accessor enforces is an invariant the platform does not
     have, because the platform reads through the trait.

     The measurement that says a low conductivity is not the same answer: three copper bars
     differing only in the middle cell — copper, the catalogue's poorest insulator, and nothing.
     The far end warms by 50 K through the first, **still warms** through the second, and sits
     exactly where it started through the third.

     A gap now **radiates**: two solid cells facing each other along a grid line across void
     exchange `σA(T₁⁴ − T₂⁴)/(1/ε₁ + 1/ε₂ − 1)`, the parallel-plate series, applied inside the
     sweep and antisymmetrically so the pair conserves to the bit. For a vacuum clearance that is
     not a correction, it is the whole answer.

     Two limits are pinned by tests rather than left to be found. The view factor is **one at
     every width**, so a gap wide compared with its faces is coupled harder here than it really
     is; the narrow gap a joint or a clearance actually is, is the case this gets right. And
     nothing **convects**, because that needs a Rayleigh number and a correlation rather than a
     closed form — so a gap in air is coupled less than it really is, and the answer there is a
     lower bound rather than an estimate. View factors would be a radiative-exchange solver, and
     that is a domain rather than a boundary condition.
   - ~~**Two parts have no way to touch.**~~ Done, and it cost one constructor plus nothing.
     `Voxels::onto` rasterises a mesh onto a grid the *caller* states, so two parts land in one
     array at the cells their own coordinates put them in — an STL carries absolute positions, so
     where they sit relative to each other is what their files already say. A mesh that would not
     fit is refused with both boxes named rather than cropped.

     The other two halves needed no new mechanism, which is the part worth recording. **Contact**
     is already there: `Solid3D`'s stencil crosses a face between cells of different materials
     with the harmonic mean of their conductivities, exact for a layered wall at every resolution.
     And the **rule for a contested cell** is `regions`' own — last writer wins, "how a coating
     *on* a layer is written and there is no other way to mean it". Inventing one here would have
     been a second answer to a question this format had already answered.
   - ~~**Nothing chooses the cell size.**~~ Done, and without the guess this entry was worried
     about. `pantometry_world::fit` **rasterises every candidate** and prints what each one cost —
     filled cells, volume error, boundary fraction, thin runs, features below the cell, undecidable
     rows — so the numbers are measurements and the only judgement left is which row to recommend.
     That rule is one sentence: the coarsest grid on which every part is present, nothing was
     undecidable, and the worst boundary fraction is under the bar. It can be disagreed with by
     reading the table it came from, which is what "guess visibly" was asking for.

     Predicting would have been wrong in a way `Loss` documents about itself: `volume_error` is a
     lattice-point count after the bulges and the cuts cancel, and a sphere at 2.5, 2.0 and 1.5 mm
     gives `+4.9%, +5.8%, −2.3%` — the first refinement makes it worse and the second changes its
     sign. A tool that extrapolated one row to the next would be confidently wrong on the commonest
     shape there is. So the rule steers by **boundary fraction**, which is the rasterisation's
     uncertainty rather than its error, and does fall monotonically.

     The ladder is a statement about the assembly rather than about millimetres: the thinnest
     dimension of any part gets 1, 2, 4, 8 … cells. Measured on a 45 x 30 x 4 mm bracket with a
     9 x 6 x 8 mm boss, the 4 mm row leaves the boss **four cells and 41% short by volume** and the
     0.5 mm row is exact with 35% of its volume in boundary cells. That first row is the failure the
     whole step exists to prevent, and it is not an error: a part finer than the grid rasterises to
     nothing and the run is perfectly well behaved about a different object.

     All three of this section's gaps are now closed.

Beyond these, the physics itself is open-ended, and two of the three that list named are done.

**Elasticity** is `pantometry-elastic`: the third elliptic domain and the first whose unknown is a
vector. **Electromagnetism** is `pantometry-em`: Maxwell on a Yee grid, which is the second hyperbolic
domain and the first to carry a constraint — `∇·B = 0` — that the update preserves as an *identity*
rather than to a tolerance. Neither needed the kernel or either layer to change at all, which is
the claim rule 4 makes and the ninth time it has held.

**Fluid dynamics** is `pantometry-fluid`, and it was the hardest of the three to make checkable for
the reason predicted: few exact solutions, schemes that trade stability against diffusion, and "it
looks like a fluid" as the easiest wrong answer in computational physics to accept. It is built
around the three that exist — Poiseuille, Couette and Taylor–Green — each chosen to be blind to a
different mistake, with two machine-precision statements beside them that a decay rate is too
coarse to see.

That closes the list this section opened with. **The physics layer is where it was aimed**: ten
crates, each one added without the kernel or either layer above it changing, which is this
document's rule 1 held ten times. `CLAUDE.md` numbers the same rule 4, which is why a reader who has
both open should trust the wording over the number.

What was open was not a list of missing physics but a list of missing *depth* in what is here, and
that list is closed too. All three entries were the second half of a crate that already existed, and
none of them cost a new one.

**Small-strain** was the last, and the answer turned out not to be finite strain. `pantometry-elastic` was
static — `∇·σ = 0`, with a `density` field its own documentation called unused — so the missing half
was **inertia**, not large deformation. `Waves` is `ρü = ∇·σ` on the same trilinear element, marched
with central differences and a lumped mass, and it is checked against two exact speeds where the
static solver is checked against four exact moduli.

The sharp statement is `c_p/c_s = √(2(1−ν)/(1−2ν))`, which holds to about `1e-8` — four orders tighter
than either speed alone, because `E` and `ρ` cancel out of it algebraically *and* the two modes share
a shape so the mesh error cancels out of it numerically. Each speed on its own is second order:
0.161%, 0.040%, 0.010% over 16, 32 and 64 elements.

Finite strain remains absent and is now a genuine *new* thing rather than a missing half. It needs
Newton iteration, and more importantly its closed forms are weak — a neo-Hookean bar's uniaxial
extension is an answer that depends on which material model was chosen, which is the closest this
workspace gets to the thing it wrote down as verifying nothing.

**Single-phase** was the second of the three, and `Solid3D` accounts for latent heat now — in **two
phases** as of the release after, where the liquid has its own conductivity and specific heat and the
check is the two-phase Neumann solution. The one-phase model is not a worse version of it: it is exact
while the liquid sits at the melting point, and switching two phases on costs first-order accuracy at
the interface, so `Substance::ice` keeps `liquid: None` and two phases are opt-in. The reduction between
them is the check on the algebra — zero superheat gives back the one-phase condition, proportionally, at
`9.4465e-3` per kelvin. `Substance`
carries `FusionProps`, `Substance::ice` is the catalogue's only entry that melts, and the check is
**Neumann's exact solution** of Stefan's problem — which is why this problem and not another: almost
nothing with a moving boundary has one.

```text
  X(t) = 2λ√(αt)      where  λ e^{λ²} erf(λ) = St/√π,   St = c(T_m − T_s)/L
```

A *position*, not a rate and not a limit. Measured at 1 mm for ice under 20 K of undercooling: 5.311
against 5.288 mm at 100 s, 10.570 against 10.567 at 400, 15.843 against 15.848 at 900. And nine times
the time gives 2.9984× the depth against `√9 = 3`, which is the `√t` law with no closed form in it at
all.

The scheme is **enthalpy**, bookkept as a temperature and a melted fraction: a cell's state is one
monotone number — `T − T_m` below, `φ·ℓ` inside, `ℓ + T − T_m` above — so the sweep adds energy and
inverts. Nothing tracks the front, nothing iterates, and energy is conserved because energy *is* the
state. It is not **apparent heat capacity**, which smears `L` over a temperature interval and lets a
step big enough to cross the interval skip the latent heat entirely; here an overshoot's remainder
lands on the far side because the inverse map says where that much energy goes, verified with ten
times the latent heat in one delivery.

Two things the tests found that are worth carrying: the profile's convergence order is **not** clean
between adjacent resolutions (0.78 then 1.84) because a fixed grid makes the front advance in a
staircase and the field depends on where it sits between cell centres — over the full fourfold
refinement it is 6.2×, between first and second order. And a checkpoint of temperatures alone is not
a checkpoint: 0 °C is ice, water or any mixture, so `checkpoint` carries the phase, which
`Schedule::Iterative` and the audit's retry both depend on.

**Single-material** was the third entry and `Solid3D::fill` is the first answer to it — and it
demonstrates the shape the rest of this list wants, because it changed nothing outside the crate
and it turned on a closed form rather than on a feature.

The number is the conductivity **on a face**, which is the harmonic mean of the two cells' and not
the arithmetic one. That is not a coarser convention: for aluminium against borosilicate the two
are 2.21 and 84.1 W/m/K, a factor of 38, and the arithmetic one short-circuits the interface. The
harmonic mean earns an *equality* — with the material interface on a cell face, the discrete chain
of face resistances is exactly `Σ Lᵢ/(kᵢA)` at every resolution, so a layered wall's resistance has
no discretisation error at all. The arithmetic mean's is first order, which is the dangerous kind:
8.9% at twelve cells, 1.0% at ninety-six, and about a thousand cells to reach in 0.1% what harmonic
has at three. A single-resolution check would have read that as a tolerance.

It also moved the stability limit off being a material property, and in the direction nobody
expects. `max_stable_dt` sums the actual face conductances now, so it is `minᵢ Cᵢ/(dx·Σ_f k_f)` —
which is `dx²/(2α)`, `dx²/(4α)` or `dx²/(6α)` according to how many axes have more than one cell
(this domain charged every shape the 3D rate, and a bar-shaped block three times the steps it
needed), and for a filled block is usually far *looser* than `dx²/(6·α_max)`. One aluminium cell
inside borosilicate is stable at exactly `k/k_face` = **75×** aluminium's own limit, because heat
cannot reach a cell faster than its worst face delivers it.

The consumer half is `regions` in the scene format and `19-a-coating-stops-the-heat`, where the
whole temperature drop sits on one face. A region that selects **no cells** is refused rather than
ignored: a mistyped bound otherwise gives a block of one material that runs, audits, renders and
answers the wrong question with nothing saying the coating was never applied.

### Every material, and the honest answer about a mixture of two

Resolving a material *per cell* answers "this part is two things". It does not answer "this part is
made of a thing that is itself two things", which is what a motor, a populated board, a printed part
and a phase-change buffer all are — and answering that turned out to need something the kernel did not
have.

**A catalogue is not the answer to "any material".** Nine entries against hundreds of thousands, and
adding a tenth does not change the shape of that. The answer is that a `Substance` is *data* —
`Deserialize` and a validator — so anything with a datasheet can be declared without the library
learning it exists. `Substance::from_name` puts the catalogue's spelling in one place so a consumer
writing a data-driven front end reads it rather than retyping it; that is the same defect three times
now, and the third copy was in the Python bindings holding five of the nine.

**`Mix` is the kernel's answer to a composite, and its shape is a refusal.** A mixture's properties
divide into three kinds and conflating them is the whole failure mode:

| | what a mixture rule can say |
| --- | --- |
| density, volumetric heat capacity, latent heat | **exact**, from conservation alone |
| conductivity, stiffness | **bounded**. No single value exists without the microstructure |
| emissivity | **nothing**. It is a property of the surface, and a mixture has no surface |

The middle row is why this belongs in the kernel rather than in a helper somebody writes once. Half
aluminium and half borosilicate conducts anywhere between 2.21 and 84.06 W/m·K — a factor of **38** —
and *both ends are attained by real arrangements of those two materials*. A library that returned
`0.5·167 + 0.5·1.114` would hand back the top of a 38-fold range as though it were a measurement. So
`Mix` returns the pair, and `as_substance` makes the caller choose inside it and refuses a choice
outside.

That this belongs to the kernel and not to a domain is worth stating as an instance of the rule: `Mix`
knows conservation and algebra and no physics. Every check on it is a domain's — a laminate in
`Solid3D` attaining Reuss across its layers and Voigt along them, a laminate in `Waves` and `Block`
attaining the shear pair, a checkerboard landing inside Hashin–Shtrikman — and none of those required
the kernel to know they existed.

**The bounds are checked against closed forms that were derived elsewhere**, which is rule 1 applied to
algebra rather than to a solver: Hashin–Shtrikman equals Maxwell–Garnett for conductivity to `6.7e-16`
and Mori–Tanaka for stiffness to `2.2e-16`, and both of those are separately derived rational functions
that agree with it at every fraction rather than in a limit.

The consumer half is `materials` and `composites` in the scene format, and scenes `21` and `22`. The
pair is the demonstration: the same wax buffer, then four fifths of it with an aluminium matrix, melting
at exactly `1/0.8` the rate — a ratio with the density and the latent heat cancelled out of it.

Two entries that used to be on this list are answered, and both the same way — not by a domain
absorbing another's job, but by a test that fails if two of them ever stop being limits of one
physics. There are three of those now:

- `fields_and_rays.rs` — a Yee grid against Fresnel's algebra, to 0.5% at eighty cells per
  wavelength, converging at second order, with a quarter-wave coating taking 14.4% to 0.23%.
- `loss_and_lumps.rs` — a field's decay in a conductor against the exact `α(ω)`, then against its
  own famous `√(2/ωμσ)` limit, beside a `Winding` answering `0.1724 Ω` at every frequency.
- `a_slit.rs` — a diffraction pattern against `sinc²`, which is **not** a comparison the field
  solution should win or lose. Scalar theory assumes Kirchhoff's boundary condition and Maxwell's
  equations do not, so the two agree in a limit: measured, the difference falls from 0.277 at one
  wavelength to 0.0057 at twelve, a factor of 48. The value is knowing *where* a closed form is
  valid, which is worth as much as the closed form.

**`pantometry-optics` having rays and no fields is the one already answered**, and not by adding fields
to it. `pantometry-em` *is* the field formulation, and `crates/pantometry/tests/fields_and_rays.rs` is where
the two have to agree: a Yee grid and Fresnel's algebra, sharing no code and not even depending on
each other, land on the same reflectance to 0.5% at eighty cells per wavelength, converging at
second order. That is the shape the rest of this list wants — not a domain absorbing another's job,
but a test that fails if they ever stop being two limits of one physics.

8. ~~**A shell that can draw a scene.**~~ Done, and it is here because it was invisible for as long
   as it existed. `pantometry view` drew `Panel::Paths` and refused everything else — *"fields and
   point clouds are different pipelines and are not built yet"* — which is a clear thing to say and
   hid how far it went: `PanelData::paths` is built by two examples and by tests and by **no
   `Domain` at all**. Not one of the thirty shipped scenes produced a panel it could draw. The
   refusal was the only sentence it ever said about a scene, and it read as a limitation rather
   than as *nothing works*.

   They are not different pipelines. A body and a field sample are points, and a point is two short
   segments in screen space, which is the pipeline that was already there — so `viewer_core::
   segments` grew two arms and the shell, the depth sort and the snapshot were not touched.
   **Twenty-seven of the thirty draw now**; the other three are a `network` and two `winding`s,
   whose domains report readings rather than places, and that set is pinned in both directions.

   The part worth keeping is what it nearly cost. *Where is sample `i` in the world* lived in
   `editor_core::field_splats`, and the viewer needed the same answer: two copies of that
   arithmetic are two pictures of one run that can disagree about where a hot spot is. It is
   `viewer_core::field_points` now, in the crate that owns the panel, with both shells walking
   through it — and a closed-form test over the mapping itself, because the pixel count that says
   a scene drew *something* passed a sabotage that put every sample at the origin.

---

## Rendering, and the decision not to compete on it

RTX rendering, GPU physics, USD and an interactive editor were asked about directly. Three of the
four are additive and only heavy; the fourth is not.

**GPU physics conflicts with rule 5**, and `app/pantometry-gpu` is the first answer to that. Atomics
and warp scheduling reorder floating-point reductions, and addition is not associative — so
bit-for-bit, which `rng::tests::the_stream_is_pinned` holds to a digest, would be given up.

The arrangement that survives it: **the CPU domain is the reference and the GPU is an
accelerator**, with a test that measures how far apart they land. The stencil itself is safe —
each cell reads six neighbours and writes itself, with no reduction and so no order to depend on
— and reductions stay on the CPU, summed in index order after a readback.

The measurement is **33–67× on a 64³ grid**, a wash at 16³, and 41–43× at 128³ — it peaks near 64³
and comes back down, because past there the device pays for cells too and the CPU sweep picks up
threads. The band is the CPU column's: the device holds `2.9e-5` to `4.0e-5` s a step over eight runs
while this laptop's CPU stencil ranges over a factor of two.

The `191×` that stood in this paragraph for months had **nothing measuring it** — the test ran at one
hard-coded grid size and printed one row. Replacing it needed the instrument fixed before the number
could be: seven sizes in one process measured the *order* of the sweep, by a factor of two to three,
because the machine drifts under load and whichever went last paid for the rest. One process per
grid, ascending and descending agreeing row for row, is the check. `app/pantometry-gpu/README.md` has it.

It comes with a cost that is not about speed: WGSL has
no `f64`, so an accelerated domain is single precision against the library's double. It conserves
to `1.45e-11` where the CPU holds `9.1e-15`, which is below `Simulation`'s default `1e-9` audit — so a
scene using one has to loosen `conservation_tolerance_for(ENERGY, ..)`, and choosing that number
is choosing what the run may lose.

That crate also found the thing worth knowing about single precision here: it was never the
problem. Storing absolute kelvin made `sum − 6·centre` subtract two numbers agreeing to five
digits, keeping less than one digit of the answer. Storing the **deviation** from the reference —
which the linear stencil commutes with exactly — improved the divergence 1660× at no cost.

The other three belong **above the layers**, in workspaces of their own, for the reason
`bindings/python` already established. Measured: the library resolves 12 external crates, the
python bindings 15, and a wgpu stack **86**. They are one workspace, `app/` — the argument is about
the boundary between the library and everything above it, and one is all it supports. See
`app/README.md`.

**A shaded viewport is not a renderer, and the distinction is the whole of this section.** The
editor draws surfaces with a depth buffer and two lights — which is what every DCC viewport does,
and what tells a reader the shape of the thing they are looking at. It is a few hundred lines of
GLSL and glow inside a workspace that already carries a GUI stack. What is *not* being built is
below: no path tracing, no global illumination, no material model. The picture that leaves this
workspace for a photoreal render leaves as USD.

And on the rest, the recommendation is not to build them. Omniverse and Isaac Sim exist and are
enormous; what this workspace has that they do not is physics that is audited, deterministic and
checked against closed forms. So the move is to be *reachable from* them — `pantometry_view::gltf`
does that in a few hundred lines and no dependency, and reaches Blender, three.js and every USD
pipeline. A USD writer is the next rung, and only worth it if the readers really are Omniverse.

**The rung was taken.** `pantometry_view::usda` writes the whole run as a USD stage, in text, still
with no dependency — the condition above is what made it wait and the condition was met. What it
buys is the thing glTF structurally cannot do: USD has time samples on any attribute, so one file
carries the topology once and the colours per frame, and usdview's timeline scrubs the simulation.
It also carries the **numbers**: a domain with no geometry at all is most of what a scene here
contains, and its scalars go out as time-sampled custom attributes under `pantometry:`.

### Where each of the four stands

| | |
| --- | --- |
| a viewer | `pantometry view` — a wgpu window that reads the run **file**. `viewer-core` does not link `pantometry`, and `the_wire_format_is_enough` holds that now the workspace boundary does not |
| export | `pantometry_view::gltf` — one frame as **surfaces**, with normals and linear colour: a field becomes the boundary of its present cells, a body a sphere. `pantometry_view::usda` — the **whole run**, animated, with every domain's scalars as time-sampled attributes. Both no dependency; the geometry is `pantometry_view::mesh`, shared, so the two cannot disagree about a solid's size. Both take a `mesh::Drawing`: which surface a field becomes, and the shapes the scene **designed**, drawn uncoloured beside the cells they became so the rasterisation loss is a picture rather than a printed number. Both carry `Panel::place` as the format spells it — a glTF node's `[x, y, z, w]`, a USD `quatd`'s `(w, x, y, z)`, measured against `pxr` rather than remembered |
| GPU physics | `app/pantometry-gpu` — 33–67× at 64³, a wash at 16³, single precision, CPU as the reference. A scene states `"device": "gpu"` and an **application** honours it through `Accelerator`, because the library's workspace cannot carry the stack |
| an editor | `pantometry edit`: the scene's JSON checked as you type, an outliner and an inspector that **writes** — every number a selected domain states, draggable, spliced back into the text byte for byte so the file keeps its own formatting — run and verify as buttons, and a **shaded viewport** — surfaces from `pantometry_view::mesh` on the GPU with a depth buffer, lit, with the extents as wire boxes occluded by what is in front of them. `viewer-core`'s camera, now also as a matrix, and `one_camera_two_paths` holds the two expressions of it to `2.4e-7` |

**The editor was last for a reason that is now gone.** An editor writes files, and the scene
format is `pantometry-world`'s, which is `publish = false` precisely because a file format is a
compatibility promise this one was not ready to make — it had already changed under its own users
once, silently, when `mode` became `release`.

It carries a `format` number now, with the narrow promise attached: **within one version, a file
that loads today loads tomorrow.** Absence means 1, which is what every scene written before the
field existed is, so old files are readable by construction rather than by special case. A version
this build does not know is refused rather than half-run — `deny_unknown_fields` catches a key that
was *added* but cannot catch one whose meaning changed, and that is exactly the gap a version
number closes.

`--check` is the other half: parse and build without running, reporting a parse failure as
`file:line:column` with the keys that were expected. That is what an editor needs while somebody is
typing, and CI runs it over every scene so it is not the one entry point nothing exercises.

What is left for an editor is a GUI, and that is a product decision rather than an architectural
one. The format is ready to be edited — and `pantometry edit` is the first GUI over it: a
workspace, linking `pantometry` where the viewer deliberately does not, reusing `viewer-core`'s
camera rather than writing that arithmetic a third time, and split into a GUI-free
`editor-core` (checked, placed, run, verified — tested headlessly) and a shell that paints only
shapes. The platform section below is the set of rules it was built under.

---

## A platform of its own, and the rules that keep it open

The section above says what not to build and where to be reachable instead, and that answer is
complete for everything an external tool can express. **Most of what this workspace holds, no
external tool can express.** There is no USD schema for a `Violation` — nor for a `Ledger`, an
impedance boundary, a multirate schedule, a radial distribution function, a mode shape, a
spectrum, a Biot number or a `Loss` report. Subtract geometry and appearance from what pantometry
models and nearly everything left has no authoring or inspection surface anywhere. Interop cannot
reach that; only a native platform can — authoring the scene format, running it, inspecting the
result, and verifying it.

So the direction is both, with the boundary stated once: **what an external tool can express is
reached through files** — glTF today, a usda writer as the next rung — **and what none of them
can express is the native platform's reason to exist.**

### Verification is the platform's identity

The most dangerous errors this workspace has found were all passed by the audit, and that list is
the platform's specification rather than a feature backlog:

- **The coupling window.** A 10 ms window delivered exactly the right joules and read the peak
  temperature 12% low — the one error class no conservation check catches.
- **A collapsed convergence order.** Refining the grid halved the error instead of quartering it,
  which is what found both acoustic boundary defects; no single run could have.
- **Rasterisation loss.** A rib finer than the grid does not fail, it disappears. `Loss` measures
  that, and `verify` now reads the measurement: every designed part gets a row in the report and
  a part the grid could not hold is a **finding**, with the exit code that goes with one.

  The verdict is `Loss::is_clean`'s own rather than a second opinion — a second threshold here
  would be a second number to keep in step, so the bar lives on the library's constant and the
  finding quotes it. What is added is the case a `Loss` structurally cannot see, because it is
  about a *part* and not about a rasterisation: **zero cells filled**. `World::build` refuses the
  assembly where every part vanished and nothing looked at one of two; that scene built, ran,
  conserved, drew and answered about an assembly with a piece missing.
- **Margins, not verdicts.** How far each dt sits from `max_stable_dt`, and how many digits each
  audit residual has left against its tolerance, per quantity — a pass with no margin is a
  different fact from a pass.

  Done, and it took a kernel change to finish. `advance` refuses in **three** places and only one
  of them can be redone from outside it: the transfer audit reads what is left on the bus between
  domains, and the per-domain books check snapshots one domain's ledger across its own turn, and
  both are gone by the time `advance` returns. So `Report` carries them out — which is the shape
  this list predicts, a platform question that turns out to need something of the kernel rather
  than something of the tool.

Every one of these needs more than the run file. Three need *reruns* — a second resolution, a
halved window, the same scene twice — and the fourth needs the **build**, because what a
rasterisation lost is not in the frames it produced: the joules the missing rib would have held
are not small in the output, they are absent from the problem. Either way a tool that only reads
a run file cannot do it. The platform links `pantometry` where the viewer deliberately does
not, and that is the structural reason the panel is native: an incumbent without the audit and
the determinism underneath cannot build it. Determinism adds one thing free of charge: a run
digest, so "this result reproduces bit for bit on any machine" is a checkable fact on a report
rather than a claim. And where a run stands beside a measured system — the digital-twin case —
the measurement takes the closed form's seat, and the same battery is the comparison.

Of the platform's verbs, all three exist inside `pantometry-world`: `--check` parses and builds
without running, reporting `file:line:column`; a scene runs to a file; and `verify` runs the
battery above — built before any GUI, because a CLI that earns trust is the platform with the
smallest possible surface. Its maiden run did what the battery is for: the default scene's peak
converged at first order where the scheme is second, and the cause measured out to be the room's
**height**, quantised to whole cells and converging on the stated 3.1 m at first order — a
contaminant the battery's own documentation names, found in the first file pointed at it.

### The role boundary: the format graduates, the stranger remains

`pantometry-world` holds three roles today, and they pull apart. It is the stranger — the consumer
whose naivety produces `FRICTION.md`, which requires it to stay unpublished. It owns the scene
format. And it is the runtime that builds, runs and draws what the format states. The moment a
platform exists, the second and third are product, and a product cannot be played by an actor
whose job is not knowing things.

The resolution is the move already made at 0.9.0, when scene and view left this same crate:
**the format and its builder graduate into a published crate, and the platform and the CLI both
consume it.** One schema, two front ends, no forks. That crate is the one place above the domains
allowed to name them — it is the composition root, which is already true today inside
`pantometry-world` rather than a permission being invented.

### Three rules, so the platform does not undo the library

1. **Its own workspace.** Measured three times now: the library resolves 12 external crates, the
   python bindings 15, the viewer's wgpu stack 86 — and a GUI stack is heavier than a viewer's.
   The library's lockfile, licence gate, WebAssembly and MSRV promises cannot carry it. Unlike
   the viewer it links `pantometry`, because it has to run things; the rendering half is
   `viewer-core`, reused rather than rewritten.

2. **The inspection half must not enumerate domains.** A panel hard-coded per domain means the
   eleventh physics costs a platform edit, and the property this document exists to protect is
   lost one layer up — the bottleneck has moved, not gone. The mechanism is the one `pantometry-view`
   already proves: dispatch on the shape of the data. A new domain offering an existing shape is
   drawn for free; a genuinely new shape costs one panel *kind*. The authoring half is different —
   it is the composition root and names domains — and the difference stays structural for the
   same reason `Pose` and `Placement::marker` do not share a type: the two halves do not share a
   crate.

3. **The stranger's method is inherited, not retired.** `FRICTION.md`'s findings have come from
   five sources; building the platform is the sixth — direct manipulation. Dragging two parts
   together arrives at "two parts have no way to touch"; a cell-size control arrives at "nothing
   chooses the cell size"; drawing an enclosure arrives at "a grid has no void" — the gaps this
   document already names, found again from outside, in the order a user meets them. The
   obligation that transfers is the writing-down: every place the format or the library cannot
   say what a user just tried to do becomes a finding before it becomes a feature.

   The sixth source has now produced one of its own, and it is a *platform* finding rather than a
   physics one: **a scene said where a part's bytes were, when it should have said which bytes.**
   `parts` named a path and the builder called `std::fs::read`, which is a sentence with no
   meaning in a browser — there is no filesystem in a tab, and the page already holds the bytes
   because somebody dropped a file on the window. Dragging a CAD export onto the editor is the
   thing a user tries first, and it could not be done at all.

   The seam is [`Parts`]: `World::build` reads from a disk and `World::build_with` reads from
   wherever it is given. What matters is not the trait but the assertion beside it — the same STL
   voxelises to the same block, cell by cell, from either source. Without that the web editor is a
   demo, which is a thing that looks like the product and answers a slightly different question.
   The general form is worth keeping: **an interface that names a location has assumed a machine.**

   [`Parts`]: app/pantometry-world/src/lib.rs

---

## How to judge a proposed change

Against the property at the top. A change is good if a new physics still costs one crate and
nothing existing moves. A change is suspect if it makes the kernel know a domain, lets a domain
see another, points an arrow upward, or replaces a checked number with a claimed one.
