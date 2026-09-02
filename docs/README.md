# Figures

Two pictures the library draws of itself, for the front page. Both are the output of an example
that CI runs on every commit, and each example is a closed-form check rather than a demonstration
that something did not crash — so the claim in a caption is guarded where it is made, by the
example, and not here.

```sh
cargo run --release --example lens_spots     -- docs/lens-achromat.svg
cargo run --release --example beam_hot_spot  -- docs/beam-hot-spot.svg
```

`editor.png` is the third and is not an example's output: it is the editor's own window, captured
with `PrintWindow` and `PW_RENDERFULLCONTENT`, which asks a window to draw itself into a bitmap
rather than reading the desktop in front of it. The window is found by walking the process's
top-level windows for the largest visible one — `MainWindowHandle` returned a 6x6 helper window on
two runs out of three, and a bitmap of that is not a screenshot. Taking it needs a display, so it
is refreshed by hand rather than by a command in this file.

That is also how to refresh them. They are the one place in this repository where generated output
is tracked on purpose: `.gitignore` still refuses assets at the root, which is where a run leaves
them, and the rule those two lines exist for — a `git add -A` that put a 302 KB filmstrip and a
927 KB frame dump into history — is untouched. These are 12 KB and 41 KB and are documentation.

They render on GitHub. They do not render on crates.io, which only shows images at absolute
`https` URLs, and `raw.githubusercontent.com` serves SVG as text rather than as an image — so
there is no URL that would work in both places, and the front page reads without them.
