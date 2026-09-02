# Figures

Two pictures the library draws of itself, for the front page. Both are the output of an example
that CI runs on every commit, and each example is a closed-form check rather than a demonstration
that something did not crash — so the claim in a caption is guarded where it is made, by the
example, and not here.

```sh
cargo run --release --example lens_spots     -- docs/lens-achromat.svg
cargo run --release --example beam_hot_spot  -- docs/beam-hot-spot.svg
```

That is also how to refresh them. They are the one place in this repository where generated output
is tracked on purpose: `.gitignore` still refuses assets at the root, which is where a run leaves
them, and the rule those two lines exist for — a `git add -A` that put a 302 KB filmstrip and a
927 KB frame dump into history — is untouched. These are 12 KB and 41 KB and are documentation.

They render on GitHub. They do not render on crates.io, which only shows images at absolute
`https` URLs, and `raw.githubusercontent.com` serves SVG as text rather than as an image — so
there is no URL that would work in both places, and the front page reads without them.
