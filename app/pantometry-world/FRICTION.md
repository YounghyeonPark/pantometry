# What the first consumer found

`pantometry-world` exists to use the SDK from outside and write down where that is awkward. A
library with no consumers is a library whose ergonomics nobody has measured, and not one of the
tests inside this repository can answer the question about itself — they are written by someone
who already knows the shape. The count used to be here and it is not a number this argument needs:
"none of them" is the claim, and a figure beside it is a figure that goes stale. It did — it said
**795** across two releases and a consolidation, while the library's suite measured 733.

Everything below was hit while building the smallest thing that loads a scene, runs it, couples
two domains over a plain channel and two more over a shared boundary, and draws the result. None of it is a bug in the physics except finding 6, which is — and which no test inside the
library could have found, because none of them was checking a rate.

**Thirty-eight of the forty-five are fixed**, and seven are recorded rather than actioned. The reasons
differ and are given in each: one because the kernel already refuses the mistake it describes,
one because it is documented rather than changed, one because the flag it wants is a breaking
change to a published crate for five readings in forty-seven, and the rest on scope. The entries are
kept rather than deleted, because what the API used to be is the argument for what it is — and because the next consumer should be able
to see that the answer to "this is awkward" was to change the library rather than to work
around it. Each fixed entry says what was done.

---

## 1. A domain cannot be built behind a `dyn`

`Simulation::with` takes `impl Domain + 'static` by value, and there is no
`impl Domain for Box<dyn Domain>`. Internally the simulation already stores
`Vec<Box<dyn Domain>>`, so the boxing happens either way — it just cannot happen on the
caller's side.

The consequence is that a data-driven builder must be a `match` with one arm per domain type,
and each arm has to call `with` separately:

```rust
sim = match spec {
    DomainSpec::Room { .. } => sim.with(AcousticRoom::of_air(..)),
    DomainSpec::Bar  { .. } => sim.with(Bar1D::new(..)),
};
```

That works and is what `World::build` does. What it forecloses is a *registry*: a third party
cannot add a domain type to the scene format without editing this match. For a physics SDK
whose central claim is that domains are pluggable, the plug is only available at compile time.

**Fixed.** The kernel gained `Simulation::with_boxed(Box<dyn Domain>)` and an
`impl Domain for Box<dyn Domain>` that delegates every method. `DomainSpec::build` now returns
a box and `World::build` is a three-line loop. The `match` over domain types still exists, but
it is confined to one function and is the *scene format's* business rather than the kernel's —
an out-of-tree domain can be boxed and added without this crate knowing.

## 2. Domain names are `&'static str`, so they cannot come from data

Every constructor takes `name: &'static str`, and `Domain::name` returns one. A name read out
of a JSON file is a `String`. `World::build` therefore calls `Box::leak`.

The leak is bounded by the number of domains in a scene, so it is survivable rather than
dangerous. But it is the API stating that names are compile-time things, and for an
application they are exactly the opposite: they are what the user typed.

This is the friction that felt worst in practice, because it is unavoidable and it appears at
the very first thing an application does.

**Fixed, and it cost less than expected.** `Domain::name` returns `&str`, every domain stores
a `String`, and every constructor takes `impl Into<String>`. Because `&str: Into<String>`,
**not one existing call site changed** — `Bar1D::new("bar", ..)` still compiles. The only
breakage in 349 tests was five comparisons against `Report::substeps`, which had to become
owned for the same reason. `Interface` followed, and `Exchange`'s spatial map is keyed by an
owned interface name now.

The `Copy` that was lost was never in a hot path: a name is read to report a violation and to
look a domain up, a handful of times per step.

## 3. Reading state back needs the concrete type

`Simulation::domain_as::<T>` needs `T`, so the renderer knows every domain type just as the
builder does. `ScalarField` is exactly the abstraction that would avoid this — sample a field,
draw it, never ask what it is — but there is no way to get a `&dyn ScalarField` from a
`&dyn Domain`.

The result is a second `match` over the same enum, in `World::capture`, for no reason other
than downcasting.

**Fixed.** `Domain::as_field` returns `Option<&dyn ScalarField>` and defaults to `None`;
`Bar1D` and `Room` implement it in one line each; `Simulation::field(name)` returns one.
`World::capture` no longer mentions `Room` or `Bar1D` at all — it asks each domain for a field
and samples it. That is what `ScalarField` was written for.

One thing the fix does not give away: a `ScalarField` is a function of position and does not
know where it stops, so the *extent* to sample over still comes from the caller. That is the
right division — a field that knew its own bounds would be a mesh — and the scene has the
bounds already.

## 4. The examples' plotting is not reachable

`crates/pantometry/examples/common/svg.rs` is about three hundred and fifty lines of dependency-
free SVG plotting, and it lives under `examples/`, so no other crate can use it. This crate
has its own smaller renderer that overlaps with it substantially.

Not obviously wrong — the examples are meant to be self-contained, and a plotting API is a
commitment. But it means the first thing a consumer wants to do after running a simulation is
something the workspace already solved and cannot share.

**Not fixed, and declined on scope rather than deferred.** Sharing it means either a `pantometry-plot`
crate or a feature-gated module in the facade, and either way it is a public API for drawing
that would have to be supported, versioned and documented — for a workspace whose stated
scope excludes rendering. Two overlapping private renderers is the cheaper mistake for now.

Revisit when there is a second consumer. One application writing its own hundred lines of SVG
is not evidence; two would be.

**There is one now, and it declined to be evidence.** The sizing tool behind findings 19–21 draws
nothing at all — it prints five numbers and exits, because a sizing question wants a settled
answer rather than a picture. So the second application did not want the plotting, which leaves
this finding exactly where it was rather than settling it. That is a real answer and not a
dodge: the case for a `pantometry-plot` needs two consumers that want to plot, and there is still
one.

## 5. `Room` is not in the prelude

`pantometry::prelude` re-exports `Tube` but not `Room`, though they are the two headline types of
the same crate. Reached through `pantometry::acoustic::Room` instead.

**Fixed.** One line. It was an oversight, as suspected.

## 6. `Room` has a first-order startup error — a real defect

This one is physics, not ergonomics, and it was found by the app checking itself against a
closed form rather than by any test in the library.

A room released in its `(1, 1)` mode should follow `|cos(2 pi f t)|` at every point. It does,
but the gap converges at **first** order against grid resolution, where the scheme's interior
is second:

```text
  31 cells   0.0528       241 cells  0.0076
  61 cells   0.0265       481 cells  0.0039
 121 cells   0.0151
```

Measured through this crate's own harness — `World::run` over 0.02 s in forty advances, each
multirate-subcycled. Stepping the room directly at its stability limit gives the same second
order and smaller absolute numbers, because the substep pattern differs. The *rate* is the
claim; the column is only reproducible where it was taken.

Halving on refinement, not quartering. The cause looks like the leapfrog's startup:
`Room::released_from` sets the velocity array to zero at `t = 0`, but a staggered scheme
carries velocity at half steps, so what is wanted is `v(-dt/2)`. For a mode released from rest
that is `-sin(pi f dt)`, not zero — an `O(dt)` error, and `dt` follows `dx` through the CFL
condition, giving exactly the first order observed.

This is the same shape as the wall-weighting defect the workspace already found and fixed: a
second-order interior dragged to first order by how the boundary — here the boundary in
*time* — is handled. It was found the same way too, by the rate rather than the size.

**Fixed, and `Tube` had it too.** The first velocity update now travels half a step; every
one after it travels a whole one. Second order, and the error at 31 cells fell by a factor of
22 — 0.0528 to 0.00238. `tests/scene.rs` and a new pair in `pantometry-acoustic` pin the rate.

Two things the fix turned up that were not visible from the outside:

- **A test had turned the bug into the specification.** `one_step_from_rest_is_the_laplacian_the_field_reports`
  asserted that one step from rest moves the pressure by `h²c²∇²p`. From rest `ṗ(0) = 0`, so
  Taylor gives `½h²c²∇²p` — the test was missing the half, and it passed because the scheme
  was missing it too. `Tube` had the matching test with the matching error. Both were written
  by reading the implementation, which is the failure mode a test written from the closed form
  does not have.
- **The old startup conserved energy *exactly*, and the fix does not.** Not a regression: with
  `v = 0` treated as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first step's
  energy change cancel to the last bit. Starting correctly breaks that cancellation by
  `−h²Σ(∇p)²/8ρ` — 0.42% of the total at 31 cells, quartering on refinement, and *only at the
  first step*; from there the invariant holds to 1e-15.

  So the old code bought exact bookkeeping by making the scheme first order. That is the
  workspace's own documented trap — "the energy functional and the update were consistent with
  each other and both wrong" — appearing a second time in the same crate.

  The energy is now reported against the released state as its datum, with the one-off
  difference kept in `Room::startup_adjustment` where it can be asked for, and bounded at 25%
  so a real first-step bug cannot hide in it.

---

## 7. The name change had missed `Bar1D::exposing`

Found by going looking, after the coupling scene needed a boundary. `Domain::name` and every
constructor took `impl Into<String>` after finding 2, but `exposing(boundary: &'static str, ..)`
did not — the sweep had matched on the parameter being called `name`. A boundary name is data
in exactly the same way and for exactly the same reason: two domains agree on it, and what they
agree on can come from a file.

**Fixed.** One signature. Worth its own entry because it is what an incomplete refactor looks
like from outside: the API is *mostly* consistent, and the one place it is not is the place
nobody had reached yet.

The crate-level documentation in `pantometry-core` was also still teaching
`fn name(&self) -> &'static str` in both of its worked examples. They compiled, so nothing
failed; they were simply showing the reader the idiom that had just been removed.

## 8. Nothing checks a schedule against the domains until the first step

A scene picks its schedule by name. `staggered` with a half-second frame is thirty-eight times
the bar's explicit-diffusion limit, and the run is refused — correctly, by name, with the
limit and the value, which is the whole argument for this library and it works.

But it is refused *when the step is taken*, not when the scene is built. `Domain::max_stable_dt`
is public, so an application can ask every domain what it can survive and refuse at build time
where the message can name the file and the line. This one does not yet.

Not a library defect. A note about where the natural seam is, and the sort of thing only
somebody loading scenes from disk would think to want.

**Half of it was already done and this file did not say so.** `verify::stability_hazard` asked
the same question of the same built, unrun world at the same `t = 0` and appended the answer to a
failed run as "likely why". So the finding was true of `--check` and of the editor and false of
`verify` — which nobody could tell from reading it here, and which is the thing a findings file is
for.

**Fixed, by moving that check rather than writing a second one.** `World::build` asks every
evolving domain for `Domain::max_stable_dt` and refuses any non-subcycling schedule whose frame
window is longer, naming the domain, both numbers, the ratio, and the frame count that would fit.
The editor checks as you type through the same `build_with`, so a scene that could never have run
says so while somebody is still writing it, and `verify` passes the same refusal through.

The battery's wording is kept — "does not subcycle", "silently unstable" — because a reader who
knows those words should meet them in the new place. And merging the two surfaced a gap in the
one being written: the first draft exempted everything but `staggered`, while the battery's
exempted only `multirate`, so a **`one-way`** scene would have been let through. It takes the
whole window in one step exactly the same way. `one_way_does_not_subcycle_either_and_is_refused_too`
is that case.

Two things the fix had to get right and one it deliberately does not.

**The suggestion is a suggestion.** "Raise `frames` to at least N" is a number, and a message
naming a number nobody tried is a plausible number. `the_frame_count_it_suggests_actually_builds`
builds at N and asserts N-1 is still refused, so N is the threshold rather than a comfortable
round-up.

**It is necessary and not sufficient, and the message says so.** `max_stable_dt` takes a `now`;
the build sees the initial state, and a domain whose limit tightens as it runs is still refused
inside `step`. Both checks exist and the one in `step` is the complete one. Putting that in the
refusal rather than only in a doc comment is the difference between a reader who knows what the
check covers and one who thinks a green `--check` is a guarantee.

**Only `staggered`.** Under `multirate` a domain subcycles the window to its own limit, so a tight
limit is a *cost* and not a refusal. Refusing there would be answering a question about affordability
with an error, which is not this check's business.

**The finding's "thirty-eight times" is in Fourier units.** Measured while writing the refusal,
which reports a ratio of *times*: the 21-cell bar in `scene.rs` is 38 in Fourier and **76.1x** in
time, and `0.5 x 76.1 = 38.05` is the relation. The 61-cell bar the new tests use is **642x**. Both
numbers are right for their own bar and neither is right for the other, which is why the refusal
prints the two times rather than a ratio alone.

All thirty shipped scenes still build. That is asserted rather than assumed: `max_stable_dt` is
state-dependent, so a domain whose limit is tightest at `t = 0` and loosens as it runs would have
been refused at build for a run it survives, and this would have been a regression rather than a
fix.

Five sabotages, four caught. The fifth — deleting the `Kind::Evolving` guard — **changed nothing**,
because every quasi-static domain in this tree returns an infinite limit and the next guard drops
it anyway. The guard stays: the trait permits a finite one, and `Simulation::sweep` gives a
quasi-static domain one step whatever the window is, so its limit could never be a reason a scene
cannot run. Recorded because a guard nothing exercises is a guard somebody deletes, and "the test
did not notice" is worth more written down than quietly left out of the tally.

## 9. A spatial coupling makes both sides state the discretisation, twice

`Bar1D::exposing(name, face_area)` builds its own `Interface` with one face per cell. A
publisher has to build a matching one, and there is no way to derive it from the bar: by the
time a `Box<dyn Domain>` exists, its interface is behind the trait, and a per-spec builder
cannot see another spec's product anyway.

So a scene says the face count twice — once as the bar's `cells` and once as the beam's
`faces` — and can say it inconsistently.

**Not fixed, and the kernel is the reason it does not need to be.** `publish_on` refuses a
flux whose face count differs from the interface's and reports both numbers. That is the right
place for the check: silently padding or truncating would put energy on the wrong part of the
boundary while keeping the total exactly right, which is the one failure a conservation audit
cannot see and the whole reason the spatial channel exists. `a_boundary_the_two_sides_cut_differently_is_refused`
asserts it.

What would remove the duplication is `exposing` taking an `Interface` rather than building
one — then a scene constructs a single boundary and hands clones to both sides. Worth doing
when a second spatial consumer exists; with one, the duplication is two integers in a file
and the kernel already refuses the mistake.

## 10. A sampled field is not the state, and averaging it is not averaging the state

The renderer samples `ScalarField` at evenly spaced points including both ends. `Bar1D`'s grid
is cell-centred, so the two end samples sit half a cell outside the outermost cell centres.
Averaging the samples therefore comes out about `1/2n` low against averaging the cells — 1.2%
at 41 cells.

Found by an assertion failing, not by reasoning: a test checked that the bar held every joule
the beam paid, computed the mean from the render panel, and missed by 1.2%.

**This said "not a defect anywhere", and that was wrong.** Fixed 2026-09-07.

The reasoning was that `ScalarField` is a function of position behaving as documented and the
renderer was sampling it as it should. Both halves are true and the conclusion does not follow,
because the *array* was never labelled as samples of a function — it is `nx * ny * nz` values
beside an `nx, ny, nz`, and every consumer in the workspace read it as the cells: the viewer, the
report, glTF, USD and the CSV. What made it visible was a six-cell block holding
`0 0 100 0 0 0 °C` written into the run file as `0 0 90 0 0 0`, with `readings` in the same frame
saying `peak 100`. One frame, two answers, and the file said nothing about which.

The fix is not a better place to sample. A field now says where its own values are —
`pantometry_core::Lattice`, `Nodal` or `Centred` — and the sampler asks. Both conventions are
real here and a single one would have been wrong for half of them: `Room` and `Hall` take
`dx = width/(nx-1)` so their outermost values sit *on* the wall where a pressure antinode is,
while `Bar1D` and `Solid3D` hold cell averages with nothing on the boundary. `PanelData::Field`
carries the answer, absent meaning nodal, so every run written before this reads back as what it
was.

It also found a second half of the same shape. A `Room` quantises its height to whole cells, so a
4.4 x 3.1 m room asked for 81 across is **3.080 m** tall — and the extent the world declared said
3.1, putting every sample in y between nodes. `03-room-pulse`'s panel peak was 0.149% below the
value its own cells held; it is exact now.

The general lesson is the one this repository keeps meeting from a new direction: **the thing that
holds a value should be the thing that says where it is.** Both halves here were something else
guessing — the sampler about the field, the world about the room.

## 11. `as_field` covers half the domains, and there is no counterpart for the other half

Finding 3 gave `Domain` an `as_field`, and a renderer stopped needing to know what a room or a
bar was. Then the scene format grew orbits, a bouncing ball and a box of atoms, and the
renderer went straight back to `domain_as::<NBody>`, `domain_as::<ContactSystem>`,
`domain_as::<Fluid>`.

Not because the fix was wrong. Those three genuinely are not fields: they are a countable
number of bodies at places, and rasterising them would invent a continuum they do not have.
`as_field` returning `None` for them is the honest answer.

**Fixed.** `Domain::as_bodies` returns `Option<&dyn Bodies>`: count, position, a value to colour
by, and a *real* wall or `None`. Sixty lines of downcasting in this crate became one call.

It sat here for months and was paid the moment the layers were split apart -- a scene layer that
must name three physics to find out where anything *is* needs editing every time a fourth
arrives, and that is the one thing the structure exists to prevent. The trait draws a line the
old code could not: a periodic cell is a boundary condition and the domain reports it, while an
orbit's box is a property of the picture and nothing physical sits at its edge.

## 12. Four mechanics domains had never opted into `as_any`

`NBody`, `TreeNBody`, `RigidBody` and `Rolling` all returned the default `None`, so
`Simulation::domain_as` could not reach any of them. The orbit scene ran, conserved, and drew
nothing at all — `bodies()` returned `None` and the panel was silently dropped.

**Fixed.** Four one-line impls. The same shape as finding 7: an opt-in that everything with
state to show is expected to take, not taken in the places nobody had needed yet. Optics,
thermal, acoustic and molecular had all taken it, because tests inside the workspace had
reached for them; mechanics had not, because none had.

Worth noticing how it failed. Not a compile error and not a violation — a picture with nothing
in it. A renderer that skips what it cannot read is reasonable on its own and produces the
least debuggable outcome there is.

## 13. `Schedule::Multirate` front-loads a coupled quantity, and the audit cannot see it

**The most serious finding in this file, and it is in the kernel.** Found by building a lumped
plate under a lamp against the *published* 0.1.0 and comparing it to the closed form of its own
scheme — not by reading the code.

`Simulation::sweep` steps one domain to completion before the next. A quasi-static publisher is
never subcycled, so it puts a whole outer step's joules on the bus once; a subcycling consumer
then calls `Exchange::take` on its **first** substep and takes all of them. Every joule of the
interval is deposited at its beginning and decays for the rest of it.

So subcycling does not refine the answer. The limit of `u ← u·gⁿ + (P·dt/C)·g^(n−1)` with
`g = 1 − h/τ` as `n → ∞` is `u·e^(−dt/τ) + (P·dt/C)·e^(−dt/τ)`, which is not the solution: the
error is first order in the **outer** step and independent of the substep entirely.

```text
  outer dt   staggered   multirate    analytic    stag err   multi err
    300 s    303.670     300.033      301.920      1.75       1.89
    150 s    302.678     301.257      301.920      0.758      0.663
     75 s    302.276     301.758      301.920      0.356      0.163
```

At 300 s the schedule chosen *for* accuracy is the worse of the two, with the errors on opposite
sides. And every one of those runs passed the conservation audit at around 1e-12: the total that
crossed is exactly right and only its distribution in time is wrong, and a `Ledger` has no
representation for *when*. This is the time-domain twin of the reason `audit_transfers` had to
become a per-face check in space.

**Fixed.** The recommendation was to document it and wait for a second consumer, and the
argument changed on inspection: a *shipped* scene already had it. `04-heater-and-bar` runs a
quasi-static heater beside a bar that subcycles hard, so this was not a speculative API.

`Exchange::take_share(channel, dt)` is the fix. `Simulation::advance` tells the bus what interval
the sweep covers, and a subcycling consumer asks for its substep's share instead of the lot.
`Bar1D` and `LumpedMass` use it. At a 300 s outer step multirate went from 1.89 K of error to
0.304 K — from the worse of the two schedules to fourteen times better than the alternative.

The share is apportioned against the time **remaining**, not against the whole interval, and that
is what makes it exact: handing out `A·dt/T` and reducing both leaves `A/T` unchanged, so the last
substep receives the remainder and the channel ends empty. Against the whole interval instead,
`n` shares leave `O(n·ε·A)` stranded, and `audit_transfers` uses an *absolute* tolerance that
would eventually refuse a run that was arithmetically fine. Even so the comparison needs a slack
of `1e-12` of the interval, because three substeps of a third do not sum to one in binary and an
exact test misses the final share.

Two things the fix exposed, both mine:

- **My first reference was wrong, and refining the step made the disagreement worse.** I compared
  against `T_a + (P/hA)(1 − e^(−t/τ))`, which is the closed form of *linear* loss, on a plate whose
  `Environment` also radiates. So the real equilibrium was lower, the run sat below the reference,
  and finer steps moved *away* from it. Explicit Euler must overshoot, so a scheme sitting below
  a reference and diverging from it on refinement is a reference that is wrong. Fixed by setting
  the emissivity to zero — which is buildable from the prelude only because of finding 15.
- **A test of mine was averaging the render panel again.** `a_heater_and_a_bar_meet_on_the_bus`
  measured the bar's mean from the sampled field rather than from the cells, and passed at 1e-6
  only because the field was nearly uniform by the end. Changing *when* the heat arrives changed
  the profile enough to expose it at 4.1e-6. That is finding 10, in a test written after finding
  10 was written down.

## 14. `Report` cannot be named without going through a module

`Report` is the return type of `Simulation::advance`, the most-called method in the library, and
it was reachable only as `pantometry::core::sim::Report`. You could use it inferred; you could not
write a function signature over it. Found by a consumer wanting a helper that takes one.

**Fixed.** Added to `pantometry-core`'s root re-export and to the prelude. Two lines.

## 15. `Substance` was in the prelude and could not be built from it

`Substance::bulk` leaves `thermal: None`, which `LumpedMass` rightly refuses to step. Supplying
one needs `ThermalProps` and two unit types, and the prelude exported twenty-five unit types
without those two. So the material set reachable from the prelude was the three catalogue entries
and one that cannot be used.

**Fixed.** `ThermalProps`, `MechanicalProps`, `AcousticProps`, `ThermalConductivity` and
`ThermalExpansion` are in the prelude.

## 16. The first error a consumer ever saw was ungrammatical

`Violation::at` builds the cases that are not a before/after comparison — a substance with no
heat capacity, an iteration that never converged — and carries a *message* in `quantity`. `Display`
had no branch for them, so it read the message as a quantity name:

```text
substance has no heat capacity is not conserved at plate: inf
```

Correct use of the constructor, correct field, unreadable sentence.

**Fixed.** A third `Display` branch for the `tolerance == 0.0 && before == after` case:
`at {site}: {quantity} ({before})`.

---

## 17. `ThermalNetwork::nodes()` returned a count that could not be turned into anything

A `Node` handle can only come from `node` or `node_losing_to`. That is deliberate and it is the
reason a link naming a node that does not exist is unrepresentable — the case the conservation
audit is structurally blind to, since a link's `+q` and `−q` cancel identically.

But it also meant a caller holding a network it did not build had no way in. `nodes()` gave a
count, `label` and `temperature` both needed a handle, and there was no way to obtain one except
`node_named`, which needs a name you do not have. The count was information you could not act on.

Found the first time this application tried to print a network's node temperatures — the one
thing the domain exists to produce.

**Fixed.** `handles() -> impl Iterator<Item = (Node, &str)>`. A dozen lines, and it makes `nodes()`
mean something.

The same shape as findings 4, 8, 9 and 10: **the API is comfortable when the parts are known at
compile time and awkward the moment they are not.** A caller writing `let winding = net.node(…)`
never noticed, because they were holding the handles already. That is now five of seventeen
findings with one underlying cause, and the count is the argument.

---

## 18. A caller can read a domain and not write one, so a feedback loop is unclosable

`Simulation::domain_as::<T>` hands back a `&T`. There was no `&mut T`.

That is fine for the reader it was built for — a test asserting a profile, a renderer sampling a
field. It is not fine for the one thing an application genuinely has to do that the library
cannot: **close a feedback loop the bus cannot carry.**

Copper's resistance rises 0.393% per kelvin, so a winding that heats up dissipates more. The
temperature lives in `pantometry-thermal`, the resistance in `pantometry-electrical`, and neither can see
the other's state — correctly, since domains meeting only on `Exchange` is the property the crate
split exists to hold. The caller between frames can see both. It could read the temperature and
had no way to write the resistance, so the loop was closable from *nowhere at all*.

**Fixed.** `Domain::as_any_mut` and `Simulation::domain_as_mut`, mirroring the existing pair.
This does not weaken the rule: it is about what happens inside `step`, and this runs between
frames in code holding `&mut Simulation` that could drop the domain and rebuild it — so denying
it a write was never protecting anything.

The same opt-in hazard as findings 7 and 12, and this time it was handled in the same change:
`as_any_mut` defaults to `None`, so a domain that forgets it is silently unwritable rather than
broken. Every domain that implements `as_any` got the counterpart beside it — twelve at the time, and
fifteen sites now including the two the application defines — along with the `Box<dyn Domain>`
forwarding impl and the kernel's own front-page example.

Scene 13 is the loop closed, and it measures what the application-level version costs.

## 19. Dimensioned constructors have no unit a person types

Building a three-node network was six lines of `Volume::from_si(x * 1e-6)`. `Length` has `mm`,
`m` and `cm`; `Volume` had only `from_si`, and `Area` had none at all.

The type system's whole promise is that a factor of a thousand appears in exactly one place — a
unit-bearing constructor. Where there is no such constructor the factor moves to the call site,
which is precisely where the promise said it would not be.

**Fixed.** `Volume::cm3`/`mm3`/`m3`/`litres`, `Area::cm2`/`mm2`/`m2`.

## 20. `runaway_current` wanted a number the network held, and the hand version was wrong

`Winding::runaway_current(g)` takes the conductance of the whole path to ambient. A caller with a
`ThermalNetwork` had to assemble it — `1/(1/K₁ + 1/K₂ + 1/(h·A))` — out of the links and
environment the network is already holding.

Tedious, and worse than tedious. That formula is convection-only: the housing also radiates at
its operating temperature, so it gave 0.203 W/K where the truth is 0.220, and a threshold of
4.11 A where the truth is 4.28. **The library's own documentation was quoting the wrong number**,
and a network with one more joint, or an environment on an interior node, would have been wrong
by more with nothing to say so.

**Fixed.** `ThermalNetwork::path_conductance(node, at)` takes the slope of its own solved
balance, so every path out is in it by construction rather than by the caller remembering.
Shipped as 0.7.0 on its own rather than batched, because a wrong number in published
documentation is something a reader copies.

## 21. The electro-thermal fixed point is eight lines every consumer writes

Dissipation depends on winding temperature, which depends on dissipation. Closing it is a loop
over `steady_state` and `dissipation_at` — eight lines, and the same eight lines for anyone who
wants a settled coupled answer.

**Not fixed, and declined rather than deferred.** The library cannot close it: a domain has no
way to read another's state, and that is the property the crate split defends. A helper would
have to depend on both crates, which only the facade may do, and a physics-specific solver in
the facade is a worse precedent than eight lines in a consumer.

The general form is a state channel on `Exchange`, which stays undecided — and this entry is
evidence *against* it. Twice the hand-written loop has agreed with the stepped answer: 99.0 °C
from the sizing tool's fixed point, 99.02 °C from scene 13's marching. Nothing has yet needed
the kernel to carry state.

---

## 22. A view wants a unit the domain does not store, and there is nowhere to say so

`Bar1D` holds kelvin. Every picture of a bar in this application is in celsius. While this crate
sampled fields itself that was four characters — an offset carried beside the extent — and the
panel came back labelled `"C"`.

Splitting the scene layer out took the sampling with it, and the offset had nowhere to go.
`pantometry-scene` cannot carry it: a conversion is a choice about presentation, and a layer that
knows no domain cannot know that this one field wants 273.15 subtracted and the pressure beside
it does not.

**Not fixed, and the shape of the answer is not obvious.** `ScalarField::unit` at least makes it
*detectable*: a panel now says `K` and a reading says `C`, and a test compares them by converting
from the two declared units rather than assuming they match. That is the honest state and it is
better than the old one, which was a silent relabelling — the offset and the label were applied
in the same expression, so nothing could have disagreed with anything.

What it wants is a view-level unit conversion with the *dimension* known, so that asking for
celsius from a pressure field fails rather than subtracting 273.15 from a pascal. The library has
that machinery — `Temperature` is a dimensioned type — and what it does not have is a way for a
`&dyn ScalarField` to say which dimension it returns. Recorded rather than guessed at.

## 23. The capture layer described a plane and called it a field

`Extent::samples` was `(usize, usize)` and `pantometry-scene`'s sampler built its position as
`(u, v, 0)`. Every field the workspace had was a line or a plane, so for six domains that was
exactly right and nobody looked at it twice.

Then `Solid3D` arrived with a field that is genuinely a volume. Nothing failed. The block was
captured as its `z = 0` face, and a 9x9 plane of a 9x9x9 aluminium block **is a perfectly
plausible picture of a block** — a hot spot in the middle, cooling at the edges, no artefact
anywhere. The report would have drawn it, the filmstrip would have drawn it, the JSON would have
carried it, and the only thing wrong was that two thirds of the samples never existed.

**Fixed.** `samples` is a triple, `PanelData::Field` carries `nz`, `Extent::volume` constructs a
box, and `Panel::slice` hands out one plane at a time so a two-dimensional view has to *ask* for
slice zero rather than get it by taking the first `nx*ny` entries.

The type system did the rest. Adding a field to a struct variant broke all three view sites --
`filmstrip`, `report` and `to_json` -- and each had to decide rather than default. The filmstrip
draws the middle slice and labels it `z-slice 5/9` in the caption; the report draws every slice
as a montage; the JSON carries all of it. A slice presented as the whole is the failure worth
preventing, and the compiler is what made three separate places confront it.

What this says about the layer split is not that `Extent` was badly designed. It is that **a
layer's assumptions are only visible from below**. `pantometry-scene` was written to name no domain
and it succeeded at that; it could not have discovered that it assumed flatness, because
everything it had ever been handed was flat. The seventh domain found it in an afternoon.

## What this says about the exercise

Forty-five findings, and the source has shifted ten times — the table below has eleven rows and
that sentence said "seven" through four of them.

| how many | where they came from |
| --- | --- |
| 1–12 | writing the application against the API |
| 13–16 | **running the two subagents built out of what the first twelve taught** — one hunting outcomes that come out empty, one building against the *published* 0.1.0 rather than the working tree |
| 17, 23 | adding a domain the library did not have, and finding the *new* API had the old shape |
| 18–22 | **splitting the application into layers**, which turns assumptions into statements |
| 24 | being asked for something the format could not express — every material rather than nine |
| 25, 26, 29 | **writing a scene about a real object** — a power module, a part under a lid — and finding the format could not pose the question |
| 27, 28, 32 | **building the browser**, which is direct manipulation: the gaps arrive in the order a user meets them |
| 30, 31 | **making an unreachable domain reachable**, which is where the layers above it show what they assumed |
| 33 | **the audit refusing three correct runs in one sitting**, all with the same shape |
| 34 | **reading a scene's own output at a scale nobody had run before** — nanoseconds and picojoules |
| 35–45 | **auditing the shipped scenes against the physics they claim** — asking of each one whether it sets up a condition anybody would recognise, rather than whether it runs |

**Splitting into layers** and **making an unreachable domain reachable** are the two rows a reader
should take away, because neither is "use the API and see what hurts". Building the next domain and
pulling out a layer are both cheap, and each finds a
class of thing the other cannot: a domain finds what the layers above it assumed, and a layer
finds what one crate doing everything had hidden.

**The last row is a third, and it is the cheapest of all.** Eight findings came from reading the
shipped scenes and asking, one at a time, whether each poses a question somebody would recognise.
Every one of those scenes ran, conserved to 1e-9, and passed a closed-form check. A power module
whose junction temperature was 32.5% wrong, a housing that was a vacuum, a shot that was a quarter
of one, a notch nobody could measure the grid of: none of them was a bug, and all of them were
wrong. What a suite checks is that the arithmetic is consistent with the file; nothing in it asks
whether the file describes anything.

Thirty-eight are fixed. That line said "ten" until a test counted them, and "twenty-eight" for
seven findings after that — the test counts the *summary* at the top of the file, and this sentence
is below it, which is the failure
`prose-auditor` exists for and the second time this file has been the one carrying it — and the
count is now checked by `friction_counts.rs` against the headings, because the author evidently
cannot do it reliably and a reader cannot do it at a glance.

Findings 11 and 22 are the same event seen twice. Pulling the scene layer into its own crate
**paid** 11, which had sat unfixed for months, and **created** 22 in the same edit: a conversion
that was invisible while one crate did everything became a thing somebody has to declare. That is
what a layer boundary does — it turns assumptions into statements, and some of the statements
turn out to be missing.

Finding 23 is the third source again, and the cheapest of the four to run: **build the next
domain**. A layer that names no domain still carries assumptions about every domain it has met,
and it cannot audit those from where it stands. `pantometry-scene` assumed fields were flat, honestly
and invisibly, until a domain with a volume was written. Nothing found that by reading.

Finding 13 is the one that changes the ledger on this exercise. Every earlier finding was
ergonomic or a defect in a domain; that one is a first-order accuracy defect in the *kernel's*
scheduler, in the schedule whose whole purpose is accuracy, invisible to the audit, and it was
found by comparing a coupled run against the closed form of its own recursion. Not by reading
the code — the code is doing exactly what it documents.

Findings 1, 2, 3 and 7 were the same shape. **The API was comfortable when the set of domains
was known at compile time and awkward the moment it was not** — and that was never a decision
anybody made. It is the shape that falls out of writing a library with no consumer, where
`&'static str` costs nothing because every name is a literal in a test.

It has been made deliberately now, in the other direction, and the cost was a tenth of the
argument for keeping it: **no existing call site changed**, five test comparisons did, and the
application lost its leak, both of its downcast matches and about forty lines.

Finding 6 is the one that matters most, and the one nothing inside could have produced. A
first-order startup error survived a second-order interior, a second-order wall fix, exact
energy conservation and 345 passing tests — two of which had turned the bug into the
specification by asserting what the implementation did. It took an outside program comparing a
released mode against `|cos(2 pi f t)|` at four grid resolutions. Nothing in the library was
checking a *rate*.

That is the case for building a consumer early, and it is stronger than the ergonomic half.
None of this was visible from inside.

## 24. A world made of nine substances

The scene format could name nine materials and no others, and this is the fifth finding of the same
shape: **the library is comfortable when the set of parts is known at compile time and awkward the
moment it is not.** Findings 1, 2, 3 and 7 were domains; this one is matter.

It was not obvious from inside, because from inside the answer looked done. `Substance` has derived
`Serialize` and `Deserialize` since the crate split, `check` validates one, and `any_material.rs`
demonstrates a datasheet material working in three domains. Every piece was there. What was missing was
a *place in the file to put one*, and nobody writing the library needed one — a test says
`Substance::bulk(...).with_thermal(...)` and the question never comes up.

Two things went wrong quietly on the way to noticing, and the second is the reason this is a finding
rather than a feature request.

`MATERIALS` was a hand-written copy of the catalogue's spelling, eight names beside nine constructors,
and the missing one was `water`. It had been unnameable from a scene since the format learned to name a
material at all — v0.3.0 through v0.13.0, **eleven releases**. Nothing could have noticed: a name absent
from a lookup is not a wrong answer, it is a substance that never appears. The lookup lives beside the
catalogue now and `MATERIALS` is an alias, and a test checks both directions, because a constructor with
no name is a defect no error message can ever report.

The other is what happens when a declaration is a mistake, and it took writing the tests to see it.
`material` on a block is **optional and defaults to aluminium**, so a scene that declares a wax and then
does not name it — a key left off, or spelled right in the declaration and wrong at the use site — runs
as a block of metal. It runs, it audits, it renders, and it answers about the wrong substance with
nothing anywhere saying so. Two hundred times the conductivity, and a picture that looks like a working
simulation throughout.

**Fixed.** `Substance::CATALOGUE` and `Substance::from_name` put the spelling beside the material, where
a consumer writing their own data-driven front end can reach it instead of retyping a nine-arm match.
`Scene.materials` is a `BTreeMap<String, Substance>` — a map so a name resolves, ordered so an error
message is the same on every platform. `Palette` refuses three things: an impossible substance,
by `check`, before any domain is built; a name that shadows the catalogue, because two files saying
`"copper"` have to mean the same copper or no comparison between two runs means anything; and a
declaration nothing used, for the reason above. "Used" has exactly one definition — it went through the
resolver — so there is no second list of the places a material name can appear, which is the defect
`MATERIALS` was just cured of one level up.

The physics is checked where the physics is. `substances_from_a_file.rs` marches gallium and
n-octadecane, declared as JSON text, against Neumann's exact solution, **and marches ice beside them
through the identical harness**: 0.039% worst for the two declared, 0.035% for the catalogue's own, over
a 21× range of Stefan number. Ice sits inside the declared range at every undercooling, which is the
only form in which the claim means anything.

---

## 25. A region could say what a box was made of and not how hot it started

Found writing a scene for a hot part under a cooled lid, which is the commonest thermal question
there is. `regions` states a material per box and `initial_c` is a property of the whole block, and
the bus — deliberately — carries an amount and no location, so heat arriving there spreads to a
uniform rise. Between them the format could not say "this corner starts at 300 °C" **in either
direction**.

The scene that exposed it ran, conserved and answered a question about a block that was uniformly
warm. Nothing was wrong; there was simply no way to pose the problem.

**Fixed.** `Region::initial_c`, applied after the fill and refused on a void region — nothing has no
temperature to start at. `verify`'s refinement carries it unscaled, because a temperature is not a
length: the same box starts at the same degrees whatever the grid.

---

## 26. A grid had no word for nothing, so a clearance had to be a bad conductor

`ARCHITECTURE.md` had already named this and it was still true in the format: a part in a box was
surrounded by another material, and insulating it meant a substance with a low conductivity — which
still conducts, still stores heat, and still sets a stability limit.

Measured on three copper bars differing only in one cell: the far end warms 50 K through copper,
still warms through the catalogue's poorest insulator, and moves only by what radiation carries
across nothing.

**Fixed.** `"material": "void"` on a region, and `void` is **reserved rather than resolved** — a
scene that declared a material by that name would be solid in one file and empty in another.

---

## 27. A scene said *where* a part's bytes were, when it should have said *which* bytes

`parts` named a path and the builder called `std::fs::read`. That is a sentence with no meaning in a
browser: there is no filesystem in a tab, and the page already **has** the bytes because somebody
dropped a file on the window. So the web editor could open a scene, run it, verify it and draw it,
and could not do the first thing anybody tries.

**Fixed.** A `Parts` trait: `World::build` reads from a disk and `World::build_with` reads from
whatever it is given. The trait is not the interesting part — the assertion beside it is, that the
same STL voxelises to the same block **cell by cell** from either source. Without that the browser
is a demo, which is a thing that looks like the product and answers a slightly different question.

The general form is worth keeping: **an interface that names a location has assumed a machine.**

---

## 28. Nothing chose the cell size, and it turned out not to need a guess

The last of `ARCHITECTURE.md`'s three assembly gaps, left open on purpose because picking a cell
"would be the first place in the workspace that *guesses* — so it has to guess visibly".

It does not have to guess. `Voxels::loss` already measures what a cell size cost, so a proposal can
rasterise at every candidate and report what happened. What is left is the *rule* for which row to
recommend, and that is one sentence a reader can disagree with.

Predicting would have been wrong in a way `Loss` documents about itself: `volume_error` is a
lattice-point count after the bulges and the cuts cancel, and a sphere at 2.5, 2.0 and 1.5 mm gives
`+4.9%, +5.8%, −2.3%` — the first refinement makes it worse and the second changes its sign.

**Fixed.** `pantometry_world::fit`, and the ladder is a statement about the assembly rather than about
millimetres: the thinnest dimension of any part gets 1, 2, 4, 8 … cells.

---

## 29. A scene could say how much heat there was and not where it was made

Every source in this format hands watts to the bus, and the bus carries an amount and no location.
That is right for what the bus is and wrong for every real thing that dissipates — a die, a winding,
a brake disc, a laser absorber all do it *somewhere*, and the gradient between there and the
heatsink is the entire question a thermal model is asked. `Solid3D::deposit` could put a joule in a
cell and nothing in the format could reach it.

**Fixed.** `dissipation`, a list of boxes and their watts, symmetric with `cooling`: one takes energy
out at a face, the other puts it in at a region. The watts are the box's **total, not a figure per
cell**, so the answer does not move when the grid does — measured, the junction moves by 2e-5 °C
when every grid doubles.

---

## 30. A body could be pushed and pulled and could not want to be a different size

`pantometry-elastic` could be loaded, clamped, pressed and prescribed, and had no way to say that a
piece of it would be larger if nothing were holding it. That is thermal expansion, and it is the
thing standing between a temperature field and a stress: the platform could compute a power
module's temperature to four figures and do nothing with it.

**Fixed.** `Block::stress_free_strain`, taking an **eigenstrain** rather than a temperature — so
swelling, curing shrinkage and a phase change are the same statement, and the domain never has to
depend on whatever computed it.

---

## 31. Four of the eleven domains could not be reached from a scene at all

Measured rather than noticed: `elastic`, `em`, `fluid` and `quantum` were referenced **zero** times
in `pantometry-world/src`. They existed as libraries with their own tests and could not be touched from
a scene file, the CLI or the browser — which is to say the platform could not ask them anything.

`scenes/README.md` had been saying "seven of the library's eleven domains" the whole time, correctly,
and nobody had read it as a gap.

**Fixed.** Four new domain kinds — `structure`, `channel`, `cavity`, `well` — each with a scene
written around a closed form the domain's own docs name. All eleven reach a scene now.

---

## 32. A caller could read seven domains and write none of them — finding 18, again

`Domain::as_any_mut` is opt-in and returns `None` by default. `pantometry-elastic::Block` implemented
`as_any` and not its mutable twin, so a coupling that wrote a temperature into it did **nothing**
and the scene reported zero strain, zero stress and zero strain energy — which reads as *no thermal
stress* rather than as *not connected*, and is the more believable of the two. Seven domains were in
that state.

This is finding 18 recurring after it was fixed, in a different set of domains, and the recurrence
is the finding: **an opt-in method with a silent default is one every later domain will forget.**

**Fixed.** All seven implement it, and `World::build` now **probes the coupling and refuses** rather
than letting a future domain fail the same way in silence.

---

## 33. A domain given energy from outside has to count it, and three of them did not

The audit stopped three correct runs in one sitting, each with the same shape: a domain holding a
conserved quantity that something outside the simulation had added to.

```text
  a block with a heat source     0 became -3.7e-11        a relative change of 1.0
  a body with an eigenstrain     5.056732 became 0.231394 a relative change of 7.1e-3
  a channel driven by a pump     0 became 2.9e-12         a relative change of 1.0
```

Two of the three had a near-zero opening balance, and that is the second half: `Ledger::add` raises
an entry's *scale* to the largest thing added to it, and the audit judges a change against that.
Writing `stored + lost − supplied` as one contribution throws it away, so the first joule of
rounding is a hundred-percent error.

**Fixed.** Each domain counts what it was given — `Solid3D::supplied`, `Block::received`,
`Channel::driven` — and each ledger adds the parts separately rather than their sum. All three are
readings too: a source nobody can see in the report is a run where 45 W and 45 mW look the same.

---

## 34. Fixed-decimal output is only readable at the scale it was chosen for

The CSV wrote every value as `{:.9}`. A cavity holding 3.2e-10 J had its entire energy history
written as a column of `0.000000000`, beside a run that had just reported the field at 921 V/m.

Then the same mistake again, an hour apart and by the same hand: the values moved to scientific
notation and the *time* column was left fixed, on the reasoning that a reader scans it for a frame.
On a 4 ns run **every timestamp printed as one of two values**, and a frequency measured 15% wrong
for no other reason.

**Fixed.** Both are `{:.9e}`. A scene format spanning nanoseconds to hours and picojoules to
megajoules has no scale a fixed format could have been chosen for.

---

## 35. A joint thinner than a cell could not be stated, and the module scene was 32.5% low for it

Every real thermal design is a chain of **contact resistances**: a die soldered to a substrate, a
substrate bonded to a baseplate, a baseplate bolted to a heatsink through grease. All of them are
tens of microns thick in a part discretised at a millimetre, and this format could state none of
them.

The reason is not an oversight, it is a floor. A `regions` entry is a box of **cells**, so the
thinnest resistance it can express is `dx/k`. Refining the mesh only lowers the floor; it does not
approach the right answer, because there is no right answer being approached — the layer's stated
thickness is whatever the grid can say. And no grid a real part can afford reaches 100 µm: on a
12 mm module that is 120 cells a side, a million times the work to carry a layer with no
interesting field inside it.

So the scenes wrote what they could. `24-a-power-module-junction-to-ambient` had its 100 µm solder
as a 1.5 mm region — **fifteen times** its own resistance — and its 0.63 mm DBC ceramic as another.
Meanwhile the **largest** resistance in a real junction-to-ambient path was absent altogether,
because a mounting is on the one face that has no cell on the other side of it and nothing in the
format described that face except the film in front of it.

The two errors pushed opposite ways and left a plausible number:

| the scene said | with the joints stated |
| --- | --- |
| 3.1379 K/W, junction at 181.19 °C | 4.1594 K/W, junction at 227.15 °C |

**32.5% low on the number the scene is named after**, and every check it had agreed with it to
`1.1e-4` — because they were all checking the same stack the file described.

**Fixed**, in two halves that are the same physics on two kinds of face. `Solid3D::joined` puts a
conductance on an interior face, in series with the two half cells already there, so
`1/k_face = 1/k_series + 1/(dx·h)`; `Solid3D::mounted_on` puts one on an outer face, between the
half cell and the film. A scene spells them `contact` and `cooling`'s `contact_w_per_m2_k`. An
unstated joint returns the harmonic mean bit for bit, a joint of zero is a clearance, and a
negative one is refused by name rather than quietly modelled as an insulator.

The joint has a resistance and **no thickness**, so the cell size can go back to being chosen for
the part.

---

## 36. Six shipped scenes solve a field on a grid and answer it as a lump

Every finding the `verify` battery could raise was about **arithmetic**: determinism, a sweep, a
drift, a rasterisation loss. None of them asked whether the scene needed the arithmetic.

A block whose cells all hold the same number has been solved as a field and answered as a lump. It
is not wrong. It passes every check it has, converges perfectly, conserves to twelve digits, and
answers a question one ordinary differential equation answers — and the report reads exactly like
the report of a scene with a real gradient.

Measured across all thirty shipped scenes, of which eleven domains report both a peak and a
coldest. The number is the last frame's `peak − coldest` over the range that domain's readings
covered across the whole run:

```text
  0.00000  20-melting-a-block-of-ice        ice
  0.00000  21-a-wax-thermal-buffer          wax
  0.00000  22-wax-in-an-aluminium-matrix    buffer
  0.00000  30-two-phases-crossing           phase_a
  0.00000  30-two-phases-crossing           phase_b
  0.00134  19-a-coating-stops-the-heat      joint
  0.02108  29-a-designed-bracket            bracket
  ----------------------------------------------- 0.05
  0.10645  24-a-power-module                module
  0.13411  15-a-hot-spot-in-a-block         block
  0.25896  25-what-140-kelvin-does          module
  0.79150  23-a-part-radiating-to-its-lid   housing
```

`30-two-phases-crossing-at-a-clearance` states two busbars of thirty-two cells apiece whose spread
is **exactly zero** — uniform copper, uniform dissipation, uniform cooling, so the field is a
constant and always was. `29-a-designed-bracket-becomes-cells` rasterises **4 100 cells** from an
STL to hold a range of 0.067 K, a Biot number of 8.6e-4. And `19-a-coating-stops-the-heat`
asserted that the largest cell-to-cell step lands on the interface, which is true and was a claim
about **0.08 K** on a block whose excursion was 60 K.

**The battery raises it now**, with the corpus above written beside the threshold, and `scene.rs`
pins the set so an arrival or a departure is a failure rather than a quieter list. The table is
the corpus the threshold was chosen against, not the state of the tree.

**The coating scene is fixed and left the list.** Its pulse was one cell at +60 K — 0.145 J spread
over 1 458 cells — and it is the whole heated face now, eighty-one times the energy. The interface
step is **9.4656 K of a 19.1152 K rise**, half the profile on one face, with the second-steepest
step at 7.3315 K. It reads 0.1062.

That bought a sweep as well. `verify` refuses to refine a `hot_spot`, because one cell halves its
physical size when the grid doubles, and it says in as many words to state the initial condition
as a region instead; as a region the bounds double with the counts. The resolution sweep runs now
and measures **0.106%** on the peak. This scene had no resolution measurement at all.

**The bracket is fixed, and it needed a capability that was missing.** A bracket is bolted at
*pads*, and both `CoolingSpec` and `Solid3D::losing_from` took a face entire — so the only way to
state a mounting was to cool the whole footprint, and then every cell sheds where it stands and
there is no route along the shape at all. The route is the only reason the shape is in the file.

Stating a smaller `area_cm2` does not stand in for it. The area is divided among the cells on the
face, so a tenth of the area is a tenth of the conductance **spread over the whole face** — the
right total in the wrong place, which is exactly the difference a shape is for.
`Solid3D::losing_from_within` and `cooling`'s `from`/`to` name a box on a face instead. The
bracket carries a module's 20 W from the tip of one arm, round the corner, into a six-by-four bolt
pad at the tip of the other: **52.4420 °C at the module against 30.5449 at the bolts**, a 21.90 K
spread on a 32.44 K rise.

**The busbars could not see each other, and fixing that nearly silenced the measurement.**
`30-two-phases-crossing-at-a-clearance` states two bars four millimetres apart, blackened,
ε = 0.9. Two surfaces at 335.6 K and 319.0 K facing across that gap exchange **6.57 mW** over the
8 by 8 mm patch where they cross — 7.64% of what the cooler bar dissipates — and two `Solid3D`
domains exchange nothing but bus totals, which carry an amount and no location. The arrangement
the scene is named for was not modelled.

`Solid3D` already computes exactly that exchange, and has since `23-a-part-radiating-to-its-lid`:
`find_gaps` walks each grid line and pairs the solid cells at either end of a run of void. It does
it **within one block**, so the fix is to state the two bars as one assembly — which is also how
this format says *two parts that interact* everywhere else.

**And that would have made this finding stop firing.** A domain's `peak` and `coldest` are the
extremes of everything in it, so one block holding two objects at different temperatures reads as
a large spread — 15.77 K — with neither bar gaining a field. A measurement a refactor can silence
is not measuring what it says, so it is **per connected body** now, six-connected over the cells a
field panel reports as finite.

Measuring per body found a scene it had been blind to. `23-a-part-radiating-to-its-lid` read 0.79
as one domain and passed; its part and its lid are two isothermal objects trading radiation, and
apart they read 0.052% and 0.049%.

**Seven bodies across five scenes are left**, and `pantometry verify` exits 1 on each. Three are
melting, where the temperature is the melting point everywhere and the answer lives in a phase
fraction the readings report only as a total — so the measurement cannot see their structure, and
the finding says as much in its own words rather than being quietly suppressed for them. Two are
that part and that lid. Two are the busbars, and **that is the physics rather than a defect**:
with every watt leaving at the ends, a 32 mm copper bar 8 mm square generating 0.344 W varies by
`P·L/(8kA)` = **53.6 mK** along its length, 0.24% of its own rise. It measures 0.585 mK. No
cooling arrangement makes that object have a field, and a report that says so is doing its job.

---

## 37. A part inside a housing could not lose heat to the air in it

`losing_from` reaches a block's six **outer** faces. A part rasterised inside a block, or a bar
with a clearance beside it, has surfaces that are interior to the grid — and nothing could reach
them. Such a part shed heat by radiating to whatever faced it across the gap, and by nothing else.

So a scene describing a housing was describing an **evacuated** one, and it failed by producing
nothing: state a film on a face made of void and `cells_on` counts no solid cell there, so the film
is charged to nobody. Measured — a copper bar walled in by void, with a film stated on the face its
neighbours occupy, sat at its initial 200 °C for the whole run, and the run completed and reported
four figures.

`23-a-part-radiating-to-its-lid` is that scene. Its part settles at **536.22 K** as shipped and at
**471.22 K** with still air in the cavity, over the same 600 s: **65 K**, on a scene whose title
says housing.

**Fixed.** `Solid3D::air_in` fills a block's void with air at a stated temperature and film, and
every solid face touching a void cell sheds `h · dx² · (T − T∞)` to it. A scene spells it `air`.
The area is the **grid's** rather than the caller's: a `cooling` entry states what an outer face
exposes because a rasterised part covers less of it than the grid does, and the faces touching an
internal void are exactly the ones the grid has.

**Convection only.** A transparent gas does not radiate, and what a surface facing a clearance
exchanges with the surface across it is the pairing `find_gaps` already computes; the two run
beside each other because both paths are real and in parallel. Scene 23's closed form is that pair
of ODEs with one term added to each body, and it agrees to **0.12%** where it agreed to 0.20%
before.

**What it does not model** is air with a state of its own: it does not warm, it does not move, and
two cavities in one block share it. That is the same model `cooling` uses for the air outside, and
it is wrong for a sealed cavity small enough that the part heats its own air.

**And it is not a replacement for a stated area.** `30-two-phases-crossing-at-a-clearance` was
tried this way and is **worse** for it: its bars are in open room air, and attaching their side
area to the block faces they touch carries convection *and* radiation over that area, where air
carries only convection. The bars came out 33.0 K against 22.0. Air models a closed cavity; a
stated area models a surface open to the room, and the Biot number of 8e-5 across a copper bar is
what makes attaching it to any subset of that bar's cells exact.

---

## 38. A design answer is a steady-state answer, and the only way to get one was to march

Every reading a designer asks for is a number a part **settles at**: a junction temperature, a
margin, a rise above ambient. There was one way to get one — run for several time constants and
look at the last two frames to decide whether it had stopped moving.

That is expensive, and it is easy to get wrong in the direction that looks like success.
`29-a-designed-bracket-becomes-cells` needed **900 s** from a cold start; nine test binaries walk
every shipped scene, so the scene walk went from 65 s to 534 s and the app gate stopped fitting
inside ten minutes. It was shortened by starting the run near the answer — which avoids the
question rather than answering it.

**Fixed.** `Solid3D::steady_state` solves the balance the march converges to, and `settle` puts a
block at it. The residual is `flux_at`, the function the sweep itself marches with, plus the source
and the clearance pairs the sweep applies beside it — so it finds the march's own fixed point
rather than a second opinion about where it is.

On the bracket: the marched answer is 52.4420 °C and the solve moves it **0.0465 K** further, in
**1.15 s** against the 470 s that march cost each of nine test binaries.

**Successive over-relaxation with Newton on the diagonal**, matrix-free.
`ThermalNetwork::steady_state` builds a dense Jacobian and factors it, which is right for a handful
of nodes and impossible here — the bracket's would be 66 million entries.

Three things it needed that plain relaxation does not, each measured failing first:

**Over-relaxation.** Gauss-Seidel's spectral radius on an `n`-cell chain is `cos²(π/2n)`, so its
sweeps grow as `n²`: a sixteen-cell bar took over two thousand and hit the bound.

**A uniform correction per body.** Relaxation converges at the rate of the *smallest* eigenvalue,
and for a block whose conduction dwarfs what it loses that eigenvalue belongs to moving the whole
field together. A 4×4×4 aluminium block losing through one face has 0.0128 W/K against 0.668 W/K a
face inside it, and its step fell from 0.104 K to 0.031 K over eighteen hundred sweeps. Face terms
cancel under a uniform shift, so what resists one is exactly what leaves the block — one scalar
equation — and a block with a clearance in it has one such mode **per connected body**.

**Damping.** The slope of a `T⁴` term is a bad guide far from the answer: a bar whose only path out
is a clearance asked for a first step of 93 000 K and the solve reported diverging. Each step is
held to a quarter of its body's own level.

**And the arrival is a measurement rather than a finding**, which took measuring to establish. A
run that stopped on the way to its answer is quoting a number its own length chose, and that reads
like a finding — it was written as one, at a percent, and then measured across the shipped scenes:

```text
   0.000%  24-a-power-module              0.000041 K
   0.020%  30-two-phases-crossing         0.004503 K
   0.211%  29-a-designed-bracket          0.046502 K
   0.540%  25-what-140-kelvin-does        0.757912 K
  10.385%  19-a-coating-stops-the-heat   18.693843 K
  13.256%  15-a-hot-spot-in-a-block       7.953824 K
  63.596%  23-a-part-radiating-to-its-lid 178.069152 K
```

The last three are transient **on purpose** — `19`'s claim is the interface step before the front
reaches the glass, and running it to steady state moves the steepest step off the interface and
breaks it. No threshold separates those from a run cut short, because what separates them is the
question the scene is asking and this format does not carry one. So `verify` reports the number and
leaves the judgement to a person.

---

## 39. The elastic model knew where it stops applying, and nothing checked a scene against it

`pantometry-elastic` is linear. It has no plasticity, it says so, and it says what that costs, in
its own words:

> `Substance` says where a material **stops coming back**; this type has no yield and no
> plasticity, so it cannot represent that and does not pretend to. A solve past yield returns a
> displacement that is arithmetically correct and physically meaningless, and nothing in the
> answer says which.

`Elastic::from_substance` drops the yield strength on the way in, so by the time a body is being
solved the number that would say has already gone. Nothing above it looked.

`25-what-140-kelvin-does-to-the-solder` ships **5.201×** past it. SAC305 is assembled at its 217 °C
reflow and sits at 40, so its free strain is `2.15e-5 × (40 − 217)` = **0.3806%** against a yield
strain of `30 MPa / 41 GPa` = **0.0732%**. The scene reports 0.2314 J of strain energy relaxing to
0.0274, which is arithmetic on a constitutive law the solder left long before — a real joint would
have yielded and crept.

**Fixed, as a measurement and a finding.** `World` keeps each element's yield strain beside the
expansion coefficients it already kept for the same reason, tracks the worst
`|free strain| / yield strain` any element reaches over the run, and `verify` reports it and raises
it above one.

**One is a physical boundary and not a chosen threshold**, which is what separates this from the
arrival measurement two findings up. There is no corpus to calibrate against and none is needed:
the model says where it stops.

The free strain is an **upper bound** on what an element carries — fully constrained it takes all of
it, free it takes none — so under one is a body certainly inside the model and over it is a body
that may not be. That asymmetry is the useful direction for a warning, and the finding says "up to".

**And a second thing fell out of it.** `DomainSpec::refined` passed a structure through unchanged,
with a comment saying it was "reported as unswept". It was not. A structure that `follows` a block
moves with it, so the block doubled, the structure did not, and `World::build` refused the pair —
the sweep reported the whole refined scene as **refused**, and this scene had been failing its own
resolution sweep since it shipped. It skips with a reason now. Refining the two in step is a
statement about two domains at once, which `DomainSpec::refined` cannot make and `Scene::refined`
could; that is left undone rather than done wrongly.

---

## 40. A reading cannot say whether it is an answer or a statement about the solve

`Reading` carries a domain, a label, a value and a unit. Nothing in it separates *what the
simulation computed* from *how well it computed it*, and both come back in the same `Vec` from the
same method. A consumer that compares readings between two runs — which is exactly what a
convergence study is — has no way to tell them apart, and will compare them all.

Mine did. `pantometry verify --deep` reruns a scene at twice the grid under a heading that says
"what moved is discretisation", and on `17-a-busbar-with-a-notch` it printed:

```text
  busbar   residual   0.000000 -> 0.000000  (281492.026%)
```

Every part of that row is wrong. The two values are 7.793e-13 and 6.644e-13 — both converged,
printed as zero by a `{:.6}` that was chosen for temperatures. The denominator is the **4.08e-17**
the residual wobbled by between the first and second frame of the base run, which is conjugate
gradients stopping at a different iterate. And the header attributes the result to the grid. The
same sweep also measured the residual "converging at order **-0.77**", in a column beside four real
ones.

**Flooring the denominator does not fix it, and measuring says why.** The wobble is 5.2e-5 of the
residual's own magnitude — far above any rounding floor, because it is real variation in a real
quantity. The quantity is simply not one that converges to anything: refining the grid changes the
iteration count, and the residual is wherever the iteration crossed its tolerance. The
discriminator is not numerical, it is what the number *means*.

`Conductor`'s own documentation is clear about which it is, and is right to keep it:

> **The residual is a reading**, not merely an internal number, and that is the point of having it
> here: an iterative solve that quietly stopped early produces a field shaped like an answer, and
> the only thing that would ever say otherwise is a column somebody can look at.

**Worked around, in the consumer, with the cost written down.** `verify::DIAGNOSTICS` names the
five labels the shipped scenes emit that describe the solve — a residual, a `divergence` a
projection removes, a `div B` a Yee grid preserves, a wavefunction norm, and a cell Reynolds
number. They print their two values and no percentage, and they are kept out of `Sweep::worst`,
which the window sweep raises a finding on: a residual that moved would otherwise have fired it
with a message about the scene's answer depending on `frames`.

A list is a shape that goes stale in silence, and it took three pins to make one safe. The scene
walk collects every `(label, unit)` the thirty scenes emit and pins all **47** against
`is_diagnostic`; it pins the length of `DIAGNOSTICS` itself, because that first pin constrains only
the *intersection* — adding `"flux"` and `"peak "`, the second a trailing-space near-miss of a live
answer label, left the whole walk green; and it pins the **six** domain/label pairs that carry a
diagnostic, because `is_diagnostic` keys on the label alone and a new domain reporting `norm` as its
*answer* would otherwise be dropped from the sweep with nothing to say so.

**Not fixed in the library, and the reason is the field list.** `Reading`'s fields are public, so
adding one is a breaking change to a published crate for a flag that five readings in forty-seven
need. The right shape is probably a constructor — `Reading::diagnostic(...)` beside
`Reading::new(...)` — with the flag behind a method, which is additive; it is worth doing at the
next break rather than on its own.

---

## 41. A scene named for an espresso shot pulled a quarter of one

`18-an-espresso-shot` ran for **eight seconds** and delivered 3.28 g from a 4.29 g dose: **4.81%
extraction**, where a shot is pulled to 18–22%. Under that the cup is sour and thin; 4.81% is not
a weak shot, it is the first fifth of one.

**Nothing failed, and the reason is the shape of every check it had.** All of them are about the
*comparison* between two baskets — a flow ratio, a strength ordering, a ring against a core — and
Darcy at a fixed pressure through a fixed bed gives a constant flow, so every one of them holds at
any length. A scene can be internally consistent, conserve to 1e-9, agree with a closed form to six
digits, and describe something nobody would recognise.

**Fixed**: 25 s, which is a shot. 10.26 g from 4.29 g is **19.90%** at 8.31% TDS, a 2.4:1 ratio. The
comparison it was built for is unchanged — the flow ratio measures 1.786775 either way — and the
channelling reads harder for it: the gapped basket pours **18.34 g at 4.82% TDS against 10.26 g at
8.31%**, nearly twice the liquid at 58% of the strength.

The check now asserts the yield lands in 18–22%, as a band rather than a number, because what makes
this a shot is that it is in the range a barista pulls to and a drift out of it either way is worth
knowing about.

---

## 42. A notch's sweep refused to run, for a reason that was true of the wrong refinement

`DomainSpec::refined` doubles every grid in a scene so `verify` can say how much of an answer is
discretisation. On a conductor with blocked cells it refused, and the refusal was reasoned:

> a blocked cell is a one-cell notch, so refining shrinks it — a different geometry, not a finer
> one

That is true of a refinement that keeps the **indices** and of no other. A cell at `(6, 0, 0)` on a
1 mm grid occupies `x ∈ [6, 7] mm`; at 0.5 mm those same millimetres are indices 12 and 13. Each
blocked cell becomes its **eight children** and the notch is the same notch — which a test now
asserts in millimetres, not in indices, because a count of eight is not the claim.

**What the refusal cost was the measurement.** `17-a-busbar-with-a-notch` is titled "the resistance
the shape actually has", and every check it had was an inequality or an identity: more than
`ρL/A`, more than a series estimate, Tellegen, the books closing. All hold at any grid. Nobody
could ask how much of the resistance was grid, because the one thing that would have asked refused
to run.

It is **4.65%**, and it is not waiting for a finer grid:

```text
      h        in-plane     R at t = 5 mm     above the limit    three-point order
    1 mm        12 x 5      1.2391920e-05        +4.6525%
    0.5         24 x 10     1.2059472e-05        +1.8449%           1.33534
    0.25        48 x 20     1.1927724e-05        +0.7323%           1.33296
    0.125       96 x 40     1.1875426e-05        +0.2906%           1.33318
    0.0625     192 x 80     1.1854670e-05        +0.1153%           1.33367
    0.03125    384 x 160    1.1846434e-05        +0.0458%
```

Current crowds into a re-entrant corner the way stress does. The conducting wedge at the notch root
has an interior angle of `3π/2`, so the potential goes as `r^λ` with `λ = π/ω = 2/3` and the flux
as `r^(-1/3)`, unbounded. A resistance is a **quadratic** functional of that field, so it converges
as `h^(2λ) = h^(4/3)` — **1.33333**, against four measured orders that bracket it rather than
descend to it.

The slot is one cell wide at the shipped grid, so the two notch roots are `h` apart there and the
corner is not resolved at all. That is where the 4.65% is.

**Fixed as a measurement, not as a grid**, and the grid runs out first anyway. `R·t` is exactly
independent of the thickness — `6.195960000000e-08` at `nz` = 1, 2, 3, 5 and 8, to twelve digits —
so the plane can be refined cheaply: 0.0625 mm costs 1.1 s that way against the 201 s the same
plane took at the shipped `nz = 5`. But 0.015625 mm is **refused**, on the audit rather than the
clock: the scene's own 1e-9
drift budget reads 1.094e-9 after 95 s, floating-point noise over 246 k cells crossing a tolerance
written for a scene a thousandth the size.

So the scene keeps its 1 mm cells and gains a check against the Richardson limit at the measured
`p = 4/3` — `1.184101e-05 ohm`, seven digits being what the finest three pairs earn — and the
discretisation error it ships with is stated rather than discovered. A separate test holds the
`4/3` itself, which is a **rate** and so the one kind of check a second copy of the same arithmetic
cannot satisfy.

**And the reason first written for it was a theorem about a different method.** `physics-checker`
confirmed `λ = 2/3`, the rate, the limit and every percentage — and found that "a resistance is an
energy-norm quantity" is the conforming-Galerkin identity, which by Dirichlet's principle makes a
voltage-driven resistor read **low** while every grid above reads high. `Conductor` is a cell-centred
conductance network, not a Galerkin method; its face currents are conservative per cell, so the
argument is the dual one and Thomson's principle makes every grid an *upper* bound. The numbers were
right, the check passed, and the sentence explaining it predicted the opposite sign from the one
every measurement showed.

---

## 43. Every scene ran at a constant drive, and no design question is constant

Thirty scenes, and every one of them turns its source on at `t = 0` and leaves it there. The
questions a design actually asks are not that shape: the junction temperature of a module under a
duty cycle, a motor at start-up, an espresso pulled with a pre-infusion. Averaging the power answers
a different question, and the difference is not small — a part that survives 50 W forever can fail
at 100 W half the time, because the peak is set by the thermal mass near the source and not by the
average.

Measured, on a heater into a bar. The same thousand joules, at two rates:

```text
  100 W for 10 s   peak 677.6442 C   final mean 350.6878 C   spent 1000 J
   50 W for 20 s   peak 521.4072 C   final mean 350.6878 C   spent 1000 J
```

The means agree to nine digits, because energy is energy. The peak is **30.0% higher**.

**Fixed**: `Scene::stages`, a load profile keyed by the domain it drives, reaching the three kinds
that have a drive and can be reached today — a `heater`'s `watts`, a `conductor`'s `volts`, a
`puck`'s `bar`. `11-motor-thermal-network` uses it, and stopped being a scene that could only
climb: at 36 W for 300 s and 12 W after, its winding peaks at **72.82 °C** and settles to 59.64,
where the constant run reported 55.04 — **17.8 K low** for choosing an insulation class. A
first-order network under a constant source cannot overshoot at all, so the peak being above the
end is a claim the key is required for.

### Three things it cost, and the one that would have been silent

**A stage between two steps.** A run advances in whole steps, so a stage asked for at 10.5 s on a
1 s step lands at 11 — half a second of extra heating, an audit that closes because the joules that
were paid were taken, and every reading a correct run of a different experiment. Refused at build.

**The refusal named a window that did not work.** It derived one as
`duration / ceil(duration / at_s)` = 10 s, and `World::steps` raises the step count to `frames`, so
the step stayed 1 s and 10.5 still landed nowhere. Searched now, and a test feeds the suggestion
back in. A refusal naming a number that does not help is worse than one naming none.

**The profile hung off `run` and not off `advance`**, which would have given the batch path one
experiment and the editor's streaming path another — the divergence `a_streamed_run_reads_back`
already exists for, one level up. Its doc comment claimed a step-by-step caller "therefore drives
the profile too", which was false when written. `World` carries the clock now.

`World::steps` also kept its own copy of the step count, and the check needs that number before a
world exists. Two copies would have drifted into exactly the failure the check is for, so there is
one and `World::steps` calls it.

### What it does not reach

`beam` and `light` have a `watts` and no `as_any_mut`, so a profile cannot reach them without an
additive change to `pantometry-optics`; `winding` has an `amps`. Each is a line of library and a
line here, left until a scene wants one rather than added on speculation.

---

## 44. The editor ran a different experiment from the CLI, and its own doc said it did not

`editor_core::run_streaming` is what the editor drives so a run can be watched while it happens.
Its documentation says:

> the last payload is **byte-identical** to what `run` returns for the same text, which the tests
> pin

Nothing pinned it. The only test on that path checked that its JSON *reads back* — a different
claim, and the one `a_streamed_run_reads_back` exists for.

What the gap hid is one line: the streaming path took `duration_s / frames` as its step, which
ignores `window_s` — **the key that exists so the step can be shorter than a frame**. Measured on a
scene whose window is a quarter of its frame:

```text
  batch     15459 bytes
  streamed  15463 bytes
  first frame  106.555400 C either way, and 106.555400 C with no window at all
```

The third line is the proof: the streaming path's answer was identical to the *unwindowed* scene,
so the step had never been shortened at all.

**No shipped scene could show it.** All thirty have `steps == frames`, measured, which is why
thirty scenes through this path in CI on every commit said nothing. A scene using `window_s` for
what it is for was one file away.

**Fixed** by giving the streaming loop `World::run`'s schedule step for step — `steps` whole steps
of `duration_s / steps`, photographed at `ceil(i * steps / frames)` by the same integer arithmetic.
`the_two_run_paths_are_one_run` compares the bytes on three scenes.

**The third of those was added because a sabotage passed.** Swapping `div_ceil` for plain integer
division left the first two green: 80 steps over 20 frames is an exact four, so the two roundings
agree there. A scene with 20 steps over 7 frames tells them apart — 3, 6, 9, 12, 15, 18, 20 against
2, 5, 8, 11, 14, 17, 20, every picture one step early and the last one right. A test that pins the
*step* is not a test that pins the *schedule*.

This is the second time this session that a load-bearing claim about these two paths was written in
a doc comment and held by nothing; `stages` was the first, and it was caught while being written
rather than after. The pattern is not the paths, it is that a sentence in a doc comment costs
nothing to write and reads exactly like a guarantee.

---

## 45. The verification battery verified a run nobody wrote

`run_measured` is the instrumented loop `pantometry verify` marches a scene through, reading the
ledger and the stability limits between advances. Its documentation says:

> The loop is `World::run`'s — advance, close the feedback, capture

It advanced `world.sim`, the `Simulation` inside the world. That does **one** of the five things
`World::advance` is: it steps the domains. It does not apply a `stages` profile, hand each
structure the stress-free strain its block's temperature implies, or solve the structure
afterwards.

Measured on a motor with a start-up load:

```text
  the battery      winding 44.801909 C
  the run          winding 59.637440 C
```

**33% apart, and the battery's number is self-consistent.** With no profile the element ran at its
start-up 36 W for the whole run, emptied its 28800 J tank at 800 s and cooled from there — every
audit margin healthy, the determinism digest stable, no finding raised. A verification battery
reporting a clean bill on a scene nobody wrote is the worst shape this repository has, because it
is the one thing that is supposed to catch the others.

**Fixed** by calling `World::advance`. A loop that *is* `World::run`'s is the only version of that
sentence which cannot go stale, and `the_battery_measures_the_run_the_world_performs` now compares
the two reading by reading.

### The third time this loop drifted

It computed its own `duration / frames` until `window_s` existed — recorded in its own comment —
and now this. The pattern across findings 43, 44 and 45 is one thing said three ways in one
session: **a sentence in a doc comment costs nothing to write and reads exactly like a
guarantee.** `stages` claimed a step-by-step caller drove the profile and was caught while being
written; `run_streaming` claimed byte-identity "which the tests pin" and nothing pinned it; this
claimed to be `World::run`'s loop while being a third of it.

What separates the three is only when they were caught, and that was luck rather than method.

---

## What this report does not cover

**All eleven domains have scenes** — findings 31 and after closed the last four. What is left is
smaller and more specific.

**`TreeNBody`, `RigidBody` and the rest of mechanics.** Four types took `as_any` in this pass
but only `NBody` and `ContactSystem` have scene variants, so Barnes-Hut and rigid rotation are
still driven only from inside.

**A room in three dimensions.** `Room` is two-dimensional by construction and the crate says
why, so this is a limit of the physics rather than of the scene format — but it does mean the
only fields anyone can draw are flat, while every *body* domain is now drawn in space.
