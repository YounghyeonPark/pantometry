# report-check

Runs the JavaScript that `pantometry-view::report` inlines into every HTML report, against a stub
canvas, and asserts on **what it drew**.

```sh
cargo run --release -p pantometry-world -- crates/pantometry-world/scenes/14-a-world.json /tmp/w.html
node tools/report-check/check.js /tmp/w.html
```

It needs node and nothing else — no npm, no `package.json`, no install step. CI runs it as its own
job over eight reports covering all six view kinds; see `.github/workflows/ci.yml`. It is **not**
part of the twenty-step gate in `CONTRIBUTING.md`, for the same reason `bindings/python` is not:
the gate is the toolchain a Rust change needs, and this needs another one.

## Why this exists

The report is about four hundred lines of JavaScript in a Rust string literal, and until this
directory existed **nothing had ever executed it**. Every test asserted on the HTML as a string:
the page contains a canvas with this `data-kind`, the page mentions that unit, the page has seven
cards. All of those pass just as well when the viewer throws on its first line and the reader gets
a page of empty boxes.

A string check also cannot tell *drew the field* from *drew nothing*. So the assertions here are:
every view got a context and made drawing calls; every view wrote its caption; every spatial axis
carries a unit; every colour bar carries numbers; hovering returns a value; the PNG button makes a
`data:` anchor; dragging a 3D view repaints it; a 2x display gets a 2x backing store at the same
page size; the arrow keys move the frame; and the animation loop keeps scheduling.

Two of those assertions were themselves wrong first, in the way an assertion usually is — by
passing. The frame-time measurement ran *after* the key presses, and a key press pauses playback,
so every report reported a median of 0.0 ms for an early return. The fix moved it earlier and
added the assertion that the frame counter actually moved. Then `step(1)` reset its clock on every
call, so eleven of twelve measurements timed nothing; the clock is monotonic now. A check that
cannot fail is worse than no check, because it is also a claim.

## The measurement that was not a measurement

The first version of this harness ran the viewer with `vm.runInContext`, which is the obvious way
to give a script a fake `document`. It reported the volume renderer at **118 ms a frame** and led
to an afternoon of optimising a renderer that was not slow.

`vm` contexts are the reason. One identical hot loop, same node, same machine:

| where it ran | time |
| --- | --- |
| main realm | 133 ms |
| `vm.createContext({})`, its own intrinsics | 3964 ms |
| `vm.createContext({ Math })`, the host's intrinsics | 3487 ms |
| defined in the context, called from outside | 3989 ms |

A 30× penalty that node's optimiser never pays back. Every timing taken through a `vm` is a timing
about node, and the two intermediate hypotheses that afternoon — that the cost was building a CSS
colour string per ray sample, then that it was writing into a host-realm typed array — were both
tested, both wrong, and both wrong *because the instrument was*.

So this harness uses `new Function` with the globals as named parameters. It stubs exactly as
completely and runs at the speed a browser would: the same volume view measures 4 ms.

The general form of that mistake is the one `CLAUDE.md` already names — read a result from the
thing that produced it. A harness is a thing that produces results, and it needs the same
suspicion as a test.

## What the stub is not

It counts calls; it does not rasterise. `fillRect` increments a counter and accumulates area, and
nothing checks that the pixels are the right colour, because there are no pixels. So this catches
a view that draws nothing, a view that throws, a missing label, an axis in the wrong units and a
frame that takes too long — and it does not catch a heatmap drawn upside down.

For that, open the file.
