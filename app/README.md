# The app: run, check, verify, view and edit

```sh
cd app
cargo run --release -- --help
cargo run --release -- run    pantometry-world/scenes/15-a-hot-spot-in-a-block.json out.html
cargo run --release -- view   out.json
cargo run --release -- edit   scene.json --run
```

One binary, named `pantometry`. It was four, in four workspaces: the scene format's CLI, a wgpu
viewer, an egui editor and a compute accelerator.

## Why one boundary and not four

The split existed for a reason that is still true, and the numbers are measured rather than argued:

| | external crates |
| --- | --- |
| the library, all seventeen published crates | **12** |
| `bindings/python`, which is split out for this | 15 |
| a wgpu window | 86 |
| a GUI shell | 371 |

`deny.toml` gates every one of the library's twelve by licence, CI builds it with `--locked`, and
the same crates compile to `wasm32` and to Rust 1.78. None of that can carry a GPU stack, a window
toolkit or a libpython link. **Nothing in `crates/` depends on anything here.** The arrow points one
way and always did.

That argument is about the boundary between the library and everything above it. It was never an
argument for three boundaries *above* it. The viewer, the editor, the CLI and the accelerator are
all things you do to a run; they share the camera, the colour scale and the scene format; and they
were held apart by nothing but the order they were written in.

### What it bought

**`"device": "gpu"` works from the command line.** `pantometry-world` refuses that key on its own —
it has no device and cannot acquire one without a dependency tree the library's promises forbid — so
honouring it needed an application that linked both the format and an accelerator. There was no such
application, so the key had never been used end to end by anything.

The first run of `one_scene_two_devices`, which is the test that could not exist before, found three
defects at once:

| | |
| --- | --- |
| the mirror was never synced | `ledger`, `readings` and `as_field` are all `&self` and `Simulation::advance` never hands back a `&mut`, so the device's cells were read back at construction and never again. `peak` read **120.0000011 °C at every one of five frames** while the CPU's spot diffused to 20.34 |
| `as_field` was missing | so `capture` made no panel: no report view, no CSV field, no glTF, no USD, no viewer. The CLI said "no field and no bodies — not drawn" about a block |
| the ledger's scale was its own net | which is `Σ cᵢ(Tᵢ − T̄)`, identically zero for equal capacities. A relative tolerance against zero is not a tolerance, and once the mirror moved it fired on a correct run at `1.487e0` |

**And three `--exclude` flags went away.** The `wasm`, `wasm-determinism` and `msrv` jobs each
carried `--exclude pantometry-world`, because WebAssembly, determinism and the 1.78 floor are
promises the *library* makes and an application should not be held to them. Everything left in that
workspace is published and every one of it is held. An exclusion is a promise with a hole in it; the
hole is a workspace boundary cargo enforces now, rather than a flag somebody has to remember on
three lines.

### What it cost, stated rather than hidden

**`pantometry-world` left the library's `deny.toml`, so this workspace grew one.** It was a member
there, so its dependencies were licence-gated with the library's twelve; it is here now beside
stacks that had never been gated at all. `app/deny.toml` gates all 374, and its allow-list was built
by censusing what is actually here rather than by copying the library's.

Its first run found **four advisories**, all reached through `eframe` and none fixable from here:
two denial-of-service vulnerabilities in `quick-xml 0.30` (via `accesskit_unix` → `atspi` →
`zbus_xml`, which parses D-Bus introspection XML from the *local* accessibility bus, and whose
requirement `^0.30` has no patched release in range), and `paste` and `ttf-parser` unmaintained —
a proc macro that never reaches a binary, and the parser for the fonts egui embeds. Each is in the
ignore list with that argument beside it rather than as a bare id. The library's own twelve remain
clean, and nothing in `crates/` depends on anything here.

**The twenty-nine scenes left the library's gate.** Their closed-form checks — scene 21's
`8.3e-15`, scene 24's `1.1e-4`, scene 25's 8.4× relaxation — are among the strongest physics tests
in the repository, and `cargo test --locked --workspace` in the root no longer runs them. They run
here. `CLAUDE.md` says both gates, not one, and CI has a job for each.

## The crates

| | |
| --- | --- |
| `pantometry-app` | the binary. `main.rs` dispatches; `cli.rs`, `view.rs` and `edit.rs` are what used to be three `main`s |
| `pantometry-world` | the scene format, the `World` builder, the verify battery, and the twenty-nine scenes |
| `pantometry-gpu` | `Solid3D`'s stencil as a compute shader, with the CPU domain as the reference. Its own README |
| `viewer-core` | reads a run **file** and turns it into something a renderer can draw. No GPU, no window |
| `editor-core` | the GUI-free half of the editor: check, place, run, verify, and the geometry the viewport draws |
| `editor-wasm` | the same editor machinery compiled for a browser |

### `viewer-core` does not link the library, and now a test says so

It is the load-bearing decision in that crate: if a viewer can be written against the file and
nothing else, the wire format is complete. It has been worth it once already — a field carried a
grid and not the extent it was sampled over, so a 9×9×9 block framed as a nine-metre cube at the
origin, and the fix went into the *format*.

That used to be a workspace boundary. In one workspace `pantometry = { workspace = true }` is one
line away in any manifest, and `editor-core` — which links the library deliberately — sits beside it
sharing the same dependency table. `the_wire_format_is_enough` checks the manifest: nothing from
`crates/`, and still two serde crates, because a viewer that had grown a dozen dependencies could
satisfy the first half and have stopped being an independent reading of the format.

## The gate

Two workspaces, two gates. This one:

```sh
cd app
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check
cargo build --locked -p editor-wasm --target wasm32-unknown-unknown
```

`cargo doc` is in it because this tree went unchecked for as long as it existed and had a broken
intra-doc link in `editor-core` from the day that file was written.

The device tests skip loudly on a machine with no adapter, and a runner has none. A software
rasteriser is deliberately not used: it would be checking a different implementation than anyone
runs. A skip that says why is a result; a skip that says nothing is the shape of a suite that has
stopped testing anything.
