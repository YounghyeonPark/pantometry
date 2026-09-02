# What is checked, and what checking it found

Every claim in this workspace is checked against a **closed form** — an exact limit, a conservation
law, a convergence rate — or against an independent computation. Never against a second copy of the
same idea, which verifies nothing. This file is the long form of that: the invariants, the places
two domains meet, and the defects that came out of taking the rule seriously.

It was the body of `README.md` until that file was cut back to what it is and how to run it.

## Two invariants

Both are now enforced rather than promised, and both survived being generalised past
optics.

**Nothing is created or destroyed without being noticed.** A `Ledger` is what a
process claims to hold, `audit` is the check, and a `Violation` names what went
missing and where. Energy crossing between domains goes over an `Exchange`, which
compares what was published against what was consumed — because each domain
conserving energy internally says nothing about the interface between two
discretisations of the same surface, and that interface is where it actually leaks.
`SurfaceOptics` still cannot be written down in a form that returns more light than
reached it.

**Nothing is random.** `Rng::for_index(seed, index)` hashes a work item into its own
stateless stream, so ray 10 000 can be drawn before ray 3 and neither result
changes. That is what lets a simulation be parallel *and* bit-reproducible — a
single shared generator loses reproducibility precisely when the run gets big enough
to need it. `rng::tests::the_stream_is_pinned` fixes the generator's output as a
constant, and changing that constant is never the fix.

The parallel half of that is executed rather than argued. `TreeNBody::with_threads`
changes how many threads evaluate forces, and
`tree::tests::parallel_and_sequential_agree_bit_for_bit` asserts the answers are
identical across one, two, four, eight and sixteen of them — through a whole
integration, not just one evaluation. It holds because each thread owns a disjoint
range of the output: there is no reduction, so there is no summation order to vary.

## Running several domains at once

Domains do not agree on how big a step is: an explicit FDTD solver on a nanometre
grid is stable to about 10⁻¹⁷ s, heat conduction to 10⁻⁹ s, rigid contact to 10⁻⁴ s,
and a thermal drift that defocuses an instrument plays out over seconds. Stepping
everything at the smallest limit integrates the slow domains ten billion times for
nothing. Three mechanisms deal with that:

- `Kind::QuasiStatic` — a domain with no state to roll forward, re-solved on demand
  rather than stepped. Light crosses an instrument in nanoseconds; against a thermal
  timescale that is zero, so optics is never integrated. The largest single saving
  available.
- `Schedule::Multirate` — each evolving domain takes as many equal substeps as its
  own stability limit needs. Integer counts from a fixed limit, never adaptive, so
  two runs follow the same arithmetic path.
- `Schedule::Iterative` — repeats the sweep until the residuals settle. Necessary
  because a staggered coupling can be unstable no matter how small the step is; the
  standard example is fluid-structure interaction at comparable densities. Failing
  to converge is a `Violation`, not a result — an unconverged coupling produces
  numbers that look like physics.

## Where two domains meet

Sharing a clock is not enough; they also have to share a *place*. `Exchange::publish`
carries an amount and nothing else, so "a 1 mm beam on a 20 mm mirror" and "one watt
spread over the whole mirror" are the same message — and a coating fails at its hot
spot, not at its average.

An `Interface` is a boundary cut into faces that both sides address, and a `Flux` is a
quantity spread over them. Both sides share one discretisation on purpose: interpolating
between two meshes is where energy quietly goes missing, so a face-count disagreement is
refused rather than papered over, and a caller who genuinely needs to cross grids says so
with `Flux::resample`, which conserves the total by construction. Spatial channels are
audited **face by face** — a redistribution that keeps the total but moves it to the wrong
end of a mirror is exactly the bug a total-only check cannot see.

Three things fell out of building it that are worth knowing.

An enthalpy's reference point is arbitrary, so it should be chosen for precision, and the
obvious choice is the bad one. Measured from absolute zero the aluminium bar in these
tests holds 1.42 kJ, so a millijoule arriving is a change in the seventh significant figure,
and differencing two such numbers leaves a rounding floor of a few times 10⁻¹² J whatever
the transfer was — worse as the grid is refined, because there are more absolute
temperatures to add up: 1.6×10⁻¹² J at 41 cells against 7.3×10⁻¹² J at 161. Measured from the initial temperature, the number being summed *is* the change.

An insulated bar under a continuous beam does not flatten out. It settles into a fixed
shape and rides upward on a mean that climbs forever, so the hot spot is permanent and the
lumped model's error in kelvin never shrinks; what shrinks is the *fraction*. The
intuition that says it evens out eventually is about a bar heated once.

And how often the domains talk is itself an error source — here the largest one, and the
only one no conservation check can catch. A 10 ms coupling window on a bar whose own step
is 0.44 ms delivers exactly the right joules and reads the peak temperature 12% low,
because the bar spends most of each window relaxing with nothing arriving. Nothing is
lost, nothing is created, both domains stay inside their stability limits, and the answer
is still wrong. The only thing that finds it is a solution that never went through the
coupling: `tests/beam_heats_where_it_lands.rs` computes one by quadrature from the
steady-state energy balance, checks *that* against the exact `ṫL²/12α` for a point
source, and then watches the coupled answer converge on it as the window shrinks.

## Closed form where there is one

`Motion` is a function of `t`: ask for the world at 0.7 s and you get it, without
having computed 0.6 s first. An exposure can therefore be sampled at seven instants
for motion blur, and frame 7 of a recording does not depend on frame 6.

That is not available in general — three bodies under gravity have no closed form,
and neither do contact or heat — so `Integrator` rolls those forward instead, under
three rules that keep the reproducibility: fixed steps, no wall clock, ordered
reduction. And for anything conservative it is `velocity_verlet` rather than `Rk4`:
RK4 is more accurate per step and steadily *dissipates*, while the symplectic method
is second order and holds its energy within a bound. The test module proves that on a
harmonic oscillator against the closed-form energy.

## What implementing the field interface found

The kind of thing that only surfaces when something uses a trait, and more of it once there
were two implementors that disagreed.

**`ScalarField::at` takes a time, and a marched domain has only "now".** A closed-form field
like `Motion` answers for any instant; `Bar1D` holds one state and no history, so the
argument is ignored. That would have made the default `rate` — a difference across two
times — read zero, which is worse than unavailable, since the bar is visibly heating. The
fix is not a workaround: a diffusive field's time derivative *is* `α∇²T`, so the governing
equation supplies from the present state exactly what the finite difference wanted history
for. And it is not approximately right. The explicit update is `T += α·dt·∇²T` on the same
stencil, so the rate the field reports is bit-for-bit the step the domain is about to take.

**A field's gradient cannot always be the derivative of its own values.** `at` interpolates
linearly between cell centres, so its exact second derivative is zero between nodes and
infinite at them — a Laplacian read off the interpolant is useless, and must come from the
cell stencil instead. The cost is that `gradient` is the derivative the *scheme* uses rather
than of `at` exactly. That is a property of sampling a discrete field, not a rough edge that
could be polished out, and it is documented where a caller meets it.

**`rate` is not always about the next step.** `Bar1D` marches forward, so the rate it reports
is exactly the change the next step will make. `Room` is a leapfrog and stores its velocities
half a step behind its pressures, so the velocity for the next update does not exist yet —
what it can report exactly is `(pⁿ − pⁿ⁻¹)/h`, the step just taken. Being half a step behind
is the better trade rather than a compromise: a centred difference at the midpoint is second
order where a forward one at the present instant is first.

**And a claim of mine that was simply wrong.** I wrote that the mirrored stencil gives a zero
gradient at an insulated wall. For `Bar1D` it does not: its grid is cell-centred, so the first
sample sits half a cell inside the wall where the temperature really is still changing.
`Room`'s grid is node-centred, a node sits *on* the wall, and there the gradient is zero to
the last bit. Same physics, different sampling — and only the second implementor could have
shown that the first one's behaviour was about the grid rather than about the boundary.

Two smaller things a caller will meet. A grid-backed field interpolates in `at` but snaps to
the nearest node for its derivatives, because the derivatives are the scheme's own stencils
and the scheme has values only at nodes — so combining them at an arbitrary point divides a
number computed here by one computed slightly over there. And a signed field needs a
different colour map from a positive one: a temperature climbs from a floor, a pressure swings
about zero, and a ramp built for the first draws a room at rest as though half of it were
cold. The report chooses between them from the run's own range — `pantometry-view::ramp` — and
puts the neutral colour at the value zero rather than at the middle of the range, which for
−100 to +300 are a quarter of the scale apart.

## A boundary defect, and the thing that found it

Both `Tube` and `Room` put their pressure samples *on* the walls rather than half a cell
inside, so a wall sample owns half a cell. Both divided its divergence by the whole `dx`. The
walls came out twice as heavy as they are, and every mode read low: 1.4% on an 89-cell room,
5.4% on a coarse one.

A percent and a half looks like discretisation error and invites a loose tolerance. What
showed it was not: refining the grid **halved** the error instead of quartering it. A scheme
second order in the interior converging at first order overall means the boundary is first
order, and that is a wrong condition rather than a coarse one. The order of convergence found
it; no single run could have.

Two independent confirmations that the factor of two was the right one. The mirrored
five-point stencil — the standard second-order treatment of a zero-gradient boundary — gives
`2(p₁ − p₀)/dx²` at a wall, and the scheme was producing exactly half of that; with the fix
they agree, which is now asserted directly by stepping once from rest and comparing against
`ScalarField::laplacian` at every node. And the conservation audit still holds at 1e-15,
because the energy carries the same half weights the update does. Those two had to move
together: changing the update alone breaks the audit, and changing the audit alone hides the
bug.

The fix cost something real, which is the interesting part. `Tube`'s absorbing ends had been
built against the unweighted update, and at the full CFL limit the corrected boundary drains
a half-width cell with a factor of exactly **−1** — it inverts the wave instead of swallowing
it. Not a divergence: a perfectly stable run that quietly reflects. So `max_stable_dt` now
reports the impedance's own limit as well as the wave's, `Z·dx/2ρc²`, which halves the step
for a matched end and leaves a closed one alone.

### The same crate, the same shape, a second time

There was another one, in the boundary in *time*, and it was still there after all of the
above. `released_from` set the velocity to zero at `t = 0`, but a staggered leapfrog carries
velocity at `t = −h/2`, so the first velocity update travelled a whole step where the initial
condition entitled it to half. `O(h)`, permanent, and `h` follows `dx` — first order again,
from a scheme whose interior and whose walls were now both second order.

Found the same way, by the rate. Not by anything in the library: it took `pantometry-world`, the
first consumer, checking a released mode against `|cos(2πft)|`. The worst departure over 20 ms
went from 0.0528 to 0.00238 at 31 cells once the first update took half a step.

Two things fell out of it that are worth more than the fix.

**A test had turned the bug into the specification.** One step from rest was asserted to move
the pressure by `h²c²∇²p` — but from rest `ṗ(0) = 0`, so Taylor gives `½h²c²∇²p`. The test was
missing the half because the scheme was, and `Tube` carried the identical pair. Both were
written by reading the implementation. A test written from the closed form cannot do that,
which is why the conventions say to write it that way.

**And the old startup conserved energy exactly, while the correct one does not.** With `v = 0`
read as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first step's energy
change cancel to the last bit. Starting correctly breaks that cancellation by `O(h²)` — 0.42%
at 31 cells, quartering on refinement, and only at the first step; from there the invariant
holds to 1e-15. So the old code had bought exact bookkeeping by making the scheme first order.
That is the trap this repository already had written down — *the functional and the update
were consistent with each other and both wrong* — turning up a second time in the same crate,
and the audit that was supposed to catch it was the thing keeping it in place. The energy is
now measured against the released state as its datum, with the difference reported by
`Room::startup_adjustment` and bounded so a real first-step bug cannot hide in it.

## The domain whose answers are distributions

Every domain but one is checked against a value: a Fresnel coefficient, a mode frequency, a
temperature rise. `pantometry-molecular` mostly cannot be. A hundred atoms in a box are chaotic by
construction — perturb one coordinate in its last bit and the trajectories separate completely
inside a few hundred steps — so "close enough" is not available as a notion, and the physics
lives in what the trajectory averages to.

There turns out to be plenty of that, and all of it exact:

| Claim | Exact form |
| --- | --- |
| Equipartition | `⟨KE⟩ = (3N − 3) k_BT / 2` |
| Ideal gas, dilute | `PV = N k_BT`, departing as `1 + B₂(T)ρ` |
| Virial pressure, any density | `PV = N k_BT + ⟨Σ f·r⟩/3` |
| Lennard-Jones well | `−ε` exactly, at `2^(1/6)σ` |
| Momentum | To the last bit, by Newton's third law |
| Energy, unthermostatted | Bounded rather than drifting, because Verlet is symplectic |

The `− 3` in the first row is not a rounding. Momentum is conserved and the drift was removed,
so the centre of mass is frozen and three degrees of freedom are gone; counting all `3N` makes
every reported temperature low by `1/N`, which is a percent at a hundred atoms and looks
exactly like statistics.

Two things this domain leans on harder than any other. `Rng::for_index` keys the Langevin noise
on `(seed, step, particle)`, which is the *only* route back to a reproducible run in a chaotic
system — there is no "close enough" to fall back on. And the symplectic argument behind
`velocity_verlet`, which the kernel proves on a harmonic oscillator and which is relied on here
for a many-body potential with a truncated force.

The `melting` example is where this pays off. Three state points are equilibrated and the
radial distribution function decides which phase each one is, rather than anyone declaring it:
`g(r)` peaks at 5.1 for the crystal, 2.9 for the liquid and 1.4 for the gas, while the *number*
of first neighbours barely moves between solid and liquid — 12.4 against 12.8. Melting costs the
order and keeps the packing, which is why a liquid is nearly as dense as its solid and a gas is
a thousand times thinner.

The crystal panel is checked exactly before anything harder is attempted: an fcc lattice's
neighbour shells are at `1 : √2 : √3 : 2` holding 12, 6, 24 and 12 atoms, and the histogram
reproduces all four to 1e-12. Combinatorics rather than measurement.

And one thing the data corrected. Long-range order should be checked at the *third* shell, not
the fourth. The fourth was the obvious choice and reads 1.045 for the crystal against 0.999 for
the liquid — no discrimination at all — because thermal broadening washes out a thin shell at a
large radius long before a fat one nearer in. The third holds 24 atoms, twice as many as the
fourth, and separates them cleanly at 2.81 against 1.27.

One test in the crate is worth describing because of how it was wrong first. The ideal-gas check
asserted that doubling the density doubles the departure, ran on one seed, and passed. Across
four seeds the ratio came out 1.35, 1.61, 2.34 and 2.92 — averaging to 2.06, but landing
anywhere. A hundred atoms is a small sample and two thousand correlated snapshots are far fewer
independent ones than they look. The test now averages over seeds, and the fix for a noisy
statistical test is more samples rather than a wider tolerance.

## Watts, photons, and the chain that closes

A spectrum is a shape until it is integrated against something. `SpectralPower`
carries a shape and a total wattage, and answers two different questions from the
same distribution: `through()` gives watts and `photon_rate()` gives photons per
second. They are not proportional — a photon at 450 nm carries 1.44 times the energy
of one at 650 nm, so a milliwatt of blue is *fewer* photons than a milliwatt of red,
and every silicon detector responds to the count. Filtering moves the mean wavelength
as well as the total, which is why scaling an unfiltered photon rate by a power
fraction is wrong by nearly three times for a tungsten lamp.

That number is also the seam between the two domains, and the integration test
follows it end to end: a dichroic under a 5 W lamp absorbs 96.2 mW, which warms a
25 mm lens 5.69 K above ambient, which grows its 100 mm mount by 4.04 µm — 46% of the
depth of focus at NA 0.25 and two thirds of it at NA 0.30. Radiometry comes from the
optics crate, heat capacity and expansion from the kernel's `Substance`, depth of
focus from diffraction; the dimensions are what let them compose, and the `Exchange`
audit is what proves nothing leaked on the way.

Those figures used to be 10.0 K, 7.10 µm and 81%, and the difference is **radiation**.
`LumpedMass::equilibrium_rise` was `P/(hA)` — convection only — while `step` always
included the radiative term, so two public functions disagreed with the crate's own
physics for as long as they existed. `Environment::loss_from`'s documentation had
been explaining why that would be wrong the whole time: a black surface at room
temperature radiates about 6 W·m⁻²·K⁻¹, the same order as still-air convection. It was
reported from outside, by somebody who built a conclusion on the old number and had to
retract the mechanism, and the error was always in the reassuring direction.

## Numerical against analytic, on purpose

Two ways of computing the same optics, kept side by side so each checks the other.
`diffraction` answers what a perfect system does, in closed form from Bessel
functions. `wavefront` answers what an imperfect one does, by transforming a pupil.
Set the aberrations to zero and the second must reproduce the first — and it does, in
three independent places that share no code:

| | |
| --- | --- |
| Airy profile | The transformed pupil agrees with `[2J₁(v)/v]²` to under 0.5% of the peak, out to 3 λ/D |
| Ideal MTF | The transform of the PSF agrees with `(2/π)(arccos s − s√(1−s²))`, the closed-form autocorrelation of a disc |
| Strehl ratio | Small aberrations follow `exp(−(2πσ)²)`, the Maréchal approximation, whichever mode produced the error |

The comparison is what makes either side trustworthy. A pupil transform has a great
deal of room to be subtly wrong — a sampling factor, a sign, a normalisation — and
none of it shows up as anything but a plausible picture.

It also catches mistakes in the *questions*. The first Airy zero is at 1.22 λ/D and
at 0.61 λ/NA; those are the same zero and the factor of two is the numerical aperture.
Asking for encircled energy at the wrong one of them gives 59% instead of 84%, which
is exactly what happened the first time.

## What the audit tolerance is really measuring

A conservation audit needs a tolerance, and the right one is a property of the
*integrator*, not a concession to sloppiness. Three cases in this workspace, all
different:

- **Momentum under N-body gravity**: 1e-11 and could be tighter. Newton's third law
  is structural here, so the only loss is turning each force into an acceleration and
  back.
- **Energy across an optics-to-thermal coupling**: 1e-9. Both sides are closed-form
  evaluations, and nothing is being integrated.
- **Energy through a penalty contact**: 2e-2. Semi-implicit Euler is symplectic, so
  its energy error is meant to stay bounded — but that guarantee needs a smooth
  potential, and a contact that switches on the moment two things touch is not one.
  Each transition shifts the shadow Hamiltonian the method actually conserves, so the
  true energy takes an `O(dt)` step at every bounce and those accumulate. Resolving
  the contact more finely shrinks each step, which is why `ContactSystem` reports a
  stability limit a hundred times smaller than bare stability needs, but no step size
  recovers the bound.

The third is the interesting one, because a tolerance chosen without knowing that
would either fail on correct code or hide a real leak.
