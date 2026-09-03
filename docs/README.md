# Figures

Two pictures the library draws of itself, for the front page. Both are the output of an example
that CI runs on every commit, and each example is a closed-form check rather than a demonstration
that something did not crash — so the claim in a caption is guarded where it is made, by the
example, and not here.

```sh
cargo run --release --example lens_spots     -- docs/lens-achromat.svg
cargo run --release --example beam_hot_spot  -- docs/beam-hot-spot.svg
```

`editor.png` is the third and is not an example's output: it is the editor's own window, taken by
[`tools/screenshot/take.ps1`](../tools/screenshot/README.md).

```sh
cd app
cargo build --release --bin pantometry
powershell -File ../tools/screenshot/take.ps1
```

**It needs a display, which is why CI cannot refresh it — and why it has a caption a machine can
read.** The two figures above come from examples CI runs on every commit, so one that stopped being
true would take a failing example with it; this one would just quietly age, and every change to the
editor's interface since it was taken would have left it wrong. So the script writes `editor.txt`
beside it: the same frame through `--ui-dump`, which is the same egui layout with no window and no
GPU. `the_screenshot_shows_the_editor_as_it_is` regenerates that text and compares it.

That guard sees the *frame*, not the pixels — a change to a colour, a font or the shaded pass moves
nothing in the dump. What it holds is the failure that was actually coming.

That is also how to refresh them. They are the one place in this repository where generated output
is tracked on purpose: `.gitignore` still refuses assets at the root, which is where a run leaves
them, and the rule those two lines exist for — a `git add -A` that put a 302 KB filmstrip and a
927 KB frame dump into history — is untouched. These are 12 KB and 41 KB and are documentation.

They render on GitHub. They do not render on crates.io, which only shows images at absolute
`https` URLs, and `raw.githubusercontent.com` serves SVG as text rather than as an image — so
there is no URL that would work in both places, and the front page reads without them.
