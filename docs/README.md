# Figures

Seven. Three are pictures the library draws of itself: the output of an example that CI runs on
every commit, and each example is a closed-form check rather than a demonstration that something
did not crash — so the claim in a caption is guarded where it is made, by the example, and not
here. The other four need a GPU or a display, which CI has neither of, and each carries a
machine-readable caption instead for the reason set out below.

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

`editor.png`, `editor-protein.png`, `bench-app.png` and `protein-app.gif` are the four that are
not an example's output. The first two are the editor's own window, taken by
[`tools/screenshot/take.ps1`](../tools/screenshot/README.md); the other two are the viewer shading
a run, which needs a GPU and no display:

```sh
cargo run --release --example optical_bench bench.json
cd app && cargo run --release -- view bench.json --snapshot ../docs/bench-app.png
```

`protein-app.gif` is the third, and the same two commands, because the run it shows is an
animation and so is it:

```sh
cargo run --release --example ligand_binding closing.json
cd app && cargo run --release -- view ../closing.json --all-frames --snapshot ../docs/protein-app.gif
```

**`--all-frames` exists because of this figure**, and not only to save time. One invocation per
frame reparses the run and brings up a GPU each time — measured at 20 s a frame against an 18.2 MB
run, so forty-eight of them is sixteen minutes and a shell loop. Every other figure on this page
is refreshed by a command; a figure refreshed by a recipe is the one that goes stale. A path
ending `.gif` collects the frames and writes the animation instead of numbering PNGs.

`--snapshot` writes a PNG when the path says `.png`, a GIF when it says `.gif` and `--all-frames`
is given, and a PPM otherwise, so nothing converts anything.

**Each has a caption a machine can read, for the reason the editor's does.** `bench-app.txt` is
the geometry the render is *of*: the two radii, the mesh counts, and the size of each flat.
`protein-app.txt` is the same for the protein: the trace and tube counts, the tightest curve the
chain takes, the closest the chain comes to itself, and every constant that decides what the frame
looks like — the radius, the sides, the smoothing, the atom radii and their subdivision, and how
many frames there are and how far out they go. Each is written and compared
by its own example on the argument-less run, which is the run CI makes of every example on every
commit -- so a scene that changed under its picture fails the example rather than ageing quietly.
Neither sees any pixels, which is the same trade the editor's dump makes and the same one worth
making.

```sh
cd app
cargo build --release --bin pantometry
powershell -File ../tools/screenshot/take.ps1
powershell -File ../tools/screenshot/take.ps1 `
  -Scene pantometry-world/scenes/31-a-protein-shaking-at-body-temperature.json `
  -Out ..\docs\editor-protein.png
```

**The caption follows `-Out`.** It wrote `docs/editor.txt` whatever `-Out` said, so taking a
picture of a second scene silently replaced the first one's caption — the picture and the text
then described different scenes, and only the guard said so. It did, twice, while this second
figure was being taken.

**It needs a display, which is why CI cannot refresh it — and why it has a caption a machine can
read.** The three figures above come from examples CI runs on every commit, so one that stopped being
true would take a failing example with it; the editor's picture would just quietly age, and every
change to its interface since it was taken would have left it wrong. So the script writes `editor.txt`
beside it: the same frame through `--ui-dump`, which is the same egui layout with no window and no
GPU. `the_screenshots_show_the_editor_as_it_is` regenerates that text and compares it.

That guard sees the *frame*, not the pixels — a change to a colour, a font or the shaded pass moves
nothing in the dump. What it holds is the failure that was actually coming.

That is also how to refresh them. They are the one place in this repository where generated output
is tracked on purpose: `.gitignore` still refuses assets at the root, which is where a run leaves
them, and the rule those two lines exist for — a `git add -A` that put a 302 KB filmstrip and a
927 KB frame dump into history — is untouched. The three SVGs are 17, 53 and 41 KB, the bench still is
44 KB, and the editor's pair is 195 and 160.

**The animation is 2.2 MB and that is the largest thing in this repository by a wide margin**, so
it is worth saying what it buys and what it would take to shrink it. It is forty-eight frames of
1100x720, and `image`'s GIF encoder writes each one whole: it does not difference against the
frame before, which on a run that is mostly unchanging background is where the bytes are. Assembled
outside with a global palette and frame differencing the same forty-eight frames come to 1.21 MB,
measured — so the 0.87 MB is the price of not owning a differencing GIF encoder, and of the figure
being refreshable by the same one command as its siblings rather than by a tool that is not here.

They render on GitHub. They do not render on crates.io, which only shows images at absolute
`https` URLs, and `raw.githubusercontent.com` serves SVG as text rather than as an image — so
there is no URL that would work in both places, and the front page reads without them.
