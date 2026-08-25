# A GPU accelerator, with the CPU domain as the reference

```sh
cd app
cargo test --release -- --nocapture
```

`--release`. This README used to say `cargo test` without it, which turns out **not** to be why its
figures were wrong: measured at 64³, debug costs the CPU column 6.7× and the device column 5.9×, and
the ratio survives (37× against 33×). The per-step host work — encoder, submit, uniform write — is
Rust too.

`GpuSolid` runs `Solid3D`'s seven-point stencil as a WGSL compute shader and implements `Domain`,
so it drops into a `Simulation` like anything else.

## The one rule

**`Solid3D` is the answer and this is a cache of it.** Where the two disagree the CPU is right.
That is what makes `tests/against_the_cpu.rs` the point of the crate rather than a nicety, and it
is the only arrangement under which the library's promises survive a GPU at all.

## What it costs and what it buys

`how_much_faster` measures **one** grid, named by `PANTOMETRY_SWEEP`, and prints the adapter it ran
on. The table is a shell loop over it:

```sh
for n in 16 24 32 48 64 96 128; do
  PANTOMETRY_SWEEP=$n cargo test --release --test against_the_cpu how_much_faster -- --nocapture
done
```

One process per grid, and that is the finding rather than a detail — see below. On an **NVIDIA
GeForce RTX 4090 Laptop GPU (Vulkan)**, release, same deposit in the same cell, each row the best of
three:

| grid | steps | cpu threads | gpu s/step | cpu s/step | speedup |
| --- | --- | --- | --- | --- | --- |
| 16³ | 400 | 1 | 2.2–2.6e-5 | 2.2e-5 | **0.87–0.98×** |
| 24³ | 400 | 1 | 2.8e-5 | 7.6–8.0e-5 | 2.7–2.9× |
| 32³ | 400 | 1 | 2.4–2.7e-5 | 1.8e-4 | 6.9–7.9× |
| 48³ | 400 | 1 | 2.7–3.4e-5 | 6.1–6.2e-4 | 18–23× |
| 64³ | 400 | 2 | 2.9–4.0e-5 | 1.1–2.3e-3 | **33–67×** |
| 96³ | 100 | 4 | 7.1–7.2e-5 | 3.4–3.6e-3 | 49–50× |
| 128³ | 100 | 7 | 1.6e-4 | 6.7–6.9e-3 | 41–43× |

**The bands are the CPU column's, not the device's.** 64³ was run eight times: the device held
`2.9e-5` to `4.0e-5` s a step — ±20% — while this laptop's CPU stencil ran anywhere from `1.1e-3` to
`2.3e-3`, a factor of two, and the ratio followed it from 33× to 67×. The other rows are two runs
each, so read them as the shape rather than as six significant figures. The device column is the one
this crate controls and it is the steady one.

**It peaks around 64³ and comes back down.** The device column is flat from 16³ to 64³ — bound by
dispatch, not by the stencil — so up to there the ratio is the CPU's `n³` growing away from a
constant. Past 64³ the device starts paying for cells too, and the CPU's sweep picks up threads
(2 → 4 → 7), so both columns grow and the ratio settles in the forties.

**At 16³ it is a wash**, and the first row says so rather than rounding `0.87` to `1×`. By 24³ it is
2.8×, so the crossover is between those two.

**The cpu threads column is what the CPU actually used**, not what the machine has: `threads_for` is
`√(cells/39000)` capped at the cores, so 48³ and under run on one thread on a 32-core machine. A
"32 threads" claim was in the first draft of this table and was wrong.

### One process per grid, because seven in one measured the order instead

The first replacement for the 191× was a sweep over all seven sizes in one test. It did not work, and
how it failed is the useful part.

**This machine slows under sustained load, and both columns slow together** — so it is neither the
device nor the allocator. The drift is seconds fast. Whatever size went last was penalised by a
factor of two to three: run backwards, 128³ read **42×** where forwards it read **21×**, and 16³ went
from `2.2e-5` to `4.7e-5` s a step for identical work. Best of three did not help, because all three
reps of a row sit at the same place in the sweep. Round-robining the reps did not help either — it
only made the spread column honest about the size of the drift, which is 100% to 400%.

What is reproducible is **the first thing measured in a fresh process**, to ±5% across runs. Hence
the loop. Ascending and descending now agree row for row, which is the check that the instrument is
fixed and not merely quieter.

Two things were chased before the machine was suspected, and neither was the cause — but both were
worth keeping and are in `src/lib.rs`:

- Every `GpuSolid::new` created its own `wgpu::Instance`, requested its own adapter and compiled the
  shader again. That is now one `Shared` per process behind a `OnceLock`, and it fixed a **real**
  bug: `cargo test --release` — the command this README gives — never finished, because creating
  those concurrently blocks. It ran in 32 s with `--test-threads=1` and not at all in ten minutes
  without. It is 12 s now.
- The bind group was rebuilt every step. Nothing in it changes except which of the two cell buffers
  is read, so there are two and they are built once.

Two attempts to make the device reclaim eagerly both **wedged the suite past ten minutes** and are
recorded in the source as such: `Maintain::Poll` after every submit, and `Maintain::Wait` in a
block's drop — the second because `Wait` waits for the whole device, and other tests were still
submitting to it.

### The 191× that stood here

It said **191× at 64³**, and `ARCHITECTURE.md` said it twice and the root README once. **Nothing
measured it.** This test ran at one hard-coded grid size and printed one row; the four-row table was
a hand measurement with no test, no adapter named and no build profile.

Both columns have since got faster, the CPU far more so, which is the whole of why the ratio fell:
the old table's CPU figure is 2.85e-2 s/step at 64³ against 1.1–2.3e-3 now, and its device figure is
1.5e-4 against 2.9–4.0e-5 — so the accelerator is **4–5× faster per step** than the number that was
used to advertise it. `Solid3D` stopped cloning the whole grid every step and gained a threaded sweep; the
kernel became the conductance form, which reads five coefficient arrays where the old one read a
scalar. Beyond that the difference is not attributable, because there is no measurement of the old
code to attribute it to. That is the actual lesson.

Accuracy, after 60 steps with the field still structured:

| | |
| --- | --- |
| worst relative difference from the reference | `1.07e-7` |
| what single precision predicts over 60 steps | `7.7e-7` |
| conservation drift, CPU (`f64`) | `9.1e-15` |
| conservation drift, GPU (`f32`) | `1.45e-11` |

## How a run asks for this

A scene says so, in the domain that wants it:

```json
{
  "title": "a block that says where it runs",
  "duration_s": 0.02,
  "frames": 3,
  "conservation_tolerance": 1e-4,
  "domains": [
    { "kind": "block", "name": "part", "cells": [64, 64, 64], "cell_mm": 1.0,
      "material": "aluminium", "initial_c": 20.0, "device": "gpu" }
  ]
}
```

`conservation_tolerance` is the scene's, not the block's, and `1e-4` rather than the default `1e-9`
because single precision cannot hold that on a long run. `the_readme_scene_parses` parses this exact
block out of this file — the first draft of it put the tolerance inside the domain, where
`deny_unknown_fields` would have refused it and the README would have gone on saying so.

`pantometry-world` **refuses** that on its own, by name, and says what would honour it. It has no
device and cannot acquire one: its workspace is thirteen licence-gated crates that compile to
`wasm32` and to Rust 1.78, and a GPU stack is none of those things. So the scene carries the request
and an application honours it:

```rust
let mut world = World::build_with_accelerator(scene, &OnDisk, &pantometry_gpu::OnTheGpu)?;
```

**Never inferred, and never fallen back from.** A device is a lower-precision computation, not a
faster one, so choosing it is choosing what the run may lose — and a heuristic on grid size would
have moved half the shipped scenes onto a different physics, or silently back off the device again.
A block the device has no pass for is an error naming the reason, the same way the run file's reader
refuses a panel kind it does not know rather than skipping it.

`tests/a_scene_that_says_where_to_run.rs` runs one such scene **both ways and compares the
readings** — which is how `coldest` was found missing from the device's four scalars, so a run on
the device had been writing a CSV with a column absent and nothing saying it was absent.

## Why this is a different computation, not a faster one

WGSL has no `f64`. This is single precision against the domain's double, so the two cannot agree —
the useful question is by how much, and the tests measure it rather than asserting it away.

The consequence is not cosmetic: `Simulation`'s conservation audit defaults to a relative `1e-9`,
and `f32` cannot hold that on a long run. A scene using `GpuSolid` needs
`conservation_tolerance_for(quantity::ENERGY, ..)` loosened, and choosing that number is choosing
what the run is allowed to lose. `GpuSolid` also declines `books_balance` for the same reason.

## The finding: single precision was not the problem

The first version stored absolute kelvin and diverged from the reference by `1.4e-3` after two
hundred steps — a thousand times what accumulation predicts. The cause was not accumulation.

The update **then** was `centre + F·(sum − 6·centre)`; it is the conductance form now, for a
different reason — see below. On absolute temperatures near 293 K that `sum` is
about 1759, where `f32`'s resolution is `1.2e-4`, and the difference being extracted from it is of
order `1e-3` K. Subtracting two numbers that agree to five digits **keeps less than one digit of
the answer**, every step, forever.

The buffer holds `T − T₀` now. The same numbers are near 1 K, resolution `1.2e-7`, and the
subtraction keeps about four digits:

| | absolute | deviation | |
| --- | --- | --- | --- |
| divergence, 200 steps | `1.449e-3` | `8.7e-7` | 1660× better |
| conservation drift | `7.4e-7` | `1.2e-10` | 6300× better |

The stencil is linear, so subtracting a constant commutes with it exactly and the fix cost
nothing. Single precision was adequate; spending it on an offset nobody needed was not.

## Determinism, and what is deliberately still on the CPU

Rule 5 of `ARCHITECTURE.md` is that results are bit-for-bit across platforms and thread counts.
A GPU cannot promise that in general, and the parts of this that could break it are handled
separately:

- **The stencil is safe.** Each cell reads six neighbours and writes itself, in no particular
  order. There is no reduction and therefore no ordering to depend on.
- **Reductions are not on the GPU.** A mean summed with atomics depends on the order workgroups
  finish, and floating-point addition is not associative, so the answer would change between runs
  on one machine. `ledger` reads the grid back and sums it on the CPU in index order.

`Ensemble` solved the same problem on threads with fixed-size blocks and the same discipline would
work here. A readback is simpler and correct; a faster deterministic reduction is worth writing
when a grid is large enough to need it.

## Two buffers, not one

The stencil ping-pongs between two buffers. Reading and writing one array means some neighbours
are already updated and some are not — a Gauss-Seidel sweep pretending to be Jacobi, which is a
different scheme with a different stability limit and an update order nobody chose.

## No GPU, no test — and it says so

Every test skips with a printed reason when there is no adapter, which is the usual case on a CI
runner. A software rasteriser would be checking a different implementation than anyone runs, so
the CI job lints and builds and lets the tests skip loudly rather than passing for the wrong
reason.
