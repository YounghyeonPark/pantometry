# Figures

Three pictures the library draws of itself, for the front page. All are the output of an example
that CI runs on every commit, and each example is a closed-form check rather than a demonstration
that something did not crash — so the claim in a caption is guarded where it is made, by the
example, and not here.

**CI runs the example and not the picture.** It passes no output path, so the committed SVG is
never regenerated there and nothing compares the two: a figure that stopped matching what its
example draws would age as quietly as the screenshot below, and only a caption that stopped being
*true* takes a failing example with it. That is the half worth having and it is not the whole.
Measured when this paragraph was written: regenerating all three gives files byte-identical to
the committed ones.

```sh
cargo run --release --example lens_spots      -- docs/lens-achromat.svg
cargo run --release --example beam_hot_spot   -- docs/beam-hot-spot.svg
cargo run --release --example optical_bench      docs/bench-3d.svg
```

`optical_bench` writes an SVG because a page takes static images and nothing else. The same run
gives `bench.html`, which rotates, and `bench.gltf`, which opens in somebody's renderer; a README
can show neither, and a still projected the way the HTML viewer projects is the nearest thing that
travels.

`editor.png` and `bench-app.png` are the two that are not an example's output. The first is the
editor's own window, taken by [`tools/screenshot/take.ps1`](../tools/screenshot/README.md); the
second is the viewer shading a run, which needs a GPU and no display:

```sh
cargo run --release --example optical_bench bench.json
cd app && cargo run --release -- view bench.json --snapshot ../docs/bench-app.png
```

`--snapshot` writes a PNG when the path says `.png` and a PPM otherwise, so nothing converts
anything.

**It has a caption a machine can read, for the reason the editor's does.** `bench-app.txt` is the
geometry the render is *of*: the two radii, the mesh counts, and the size of each flat. It is
written and compared by `optical_bench` itself on the argument-less run, which is the run CI makes
of every example on every commit -- so a bench that changed under the picture fails the example
rather than ageing quietly. It sees no pixels, which is the same trade the editor's dump makes and
the same one worth making.

```sh
cd app
cargo build --release --bin pantometry
powershell -File ../tools/screenshot/take.ps1
```

**It needs a display, which is why CI cannot refresh it — and why it has a caption a machine can
read.** The three figures above come from examples CI runs on every commit, so one that stopped being
true would take a failing example with it; the editor's picture would just quietly age, and every
change to its interface since it was taken would have left it wrong. So the script writes `editor.txt`
beside it: the same frame through `--ui-dump`, which is the same egui layout with no window and no
GPU. `the_screenshot_shows_the_editor_as_it_is` regenerates that text and compares it.

That guard sees the *frame*, not the pixels — a change to a colour, a font or the shaded pass moves
nothing in the dump. What it holds is the failure that was actually coming.

That is also how to refresh them. They are the one place in this repository where generated output
is tracked on purpose: `.gitignore` still refuses assets at the root, which is where a run leaves
them, and the rule those two lines exist for — a `git add -A` that put a 302 KB filmstrip and a
927 KB frame dump into history — is untouched. The three SVGs are 17, 52 and 40 KB, the two
renders 43 and 176, and all five are documentation.

They render on GitHub. They do not render on crates.io, which only shows images at absolute
`https` URLs, and `raw.githubusercontent.com` serves SVG as text rather than as an image — so
there is no URL that would work in both places, and the front page reads without them.
