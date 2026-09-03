# Working on pantometry

This file is loaded every session, so it holds only what is needed *every* session. Everything read
once — a release, a change to the bindings, a look at the viewer — lives beside its subject, and the
table at the bottom says when to go there.

For *using* the library, read [AGENTS.md](AGENTS.md). For where it is all going, read
[ARCHITECTURE.md](ARCHITECTURE.md): three layers, the state of each, and the rules that make "add a
physics" cost one crate.

## The gate, before any commit

**Save it and run the file, or chain every step with `&&`.** `set -euo pipefail` in a pasted block
does **not** protect here, measured: `( set -e; false; echo reached )` prints `reached` and exits 0. It
works in a fresh `bash -c`, so the option is right and the paste is what defeats it.

The surest thing is one check per command with its exit code read. That is what caught the sixth time
this gate reported a pass it had not earned — three shell guards had not.
[CONTRIBUTING.md](CONTRIBUTING.md#run-what-ci-runs) is the authority and has all eight with what each
cost. The seventh is not a shell problem: `gh run list --limit 1` after a push returns the *previous*
commit's run, and with `cancel-in-progress` on this workflow that run has probably just been
**cancelled** — neither a pass nor a failure. Select a run by `headSha`. The eighth is not either, and
points the other way: after this checkout moved from `C:\dev\pantometry-core`, the gate ran test binaries
with the old path still baked in through `env!("CARGO_MANIFEST_DIR")` — a loud failure only because
the old directory was gone; had it still existed, six tests would have **passed** against the wrong
tree. `cargo clean` after moving a checkout.

The ninth reached `main` and sat there: commit `4a3654f` shipped `heat_crosses_a_join_and_not_a_gap`
failing, and it was found a commit later by a clean worktree at that sha rather than by the gate.
The mechanism is the tool boundary, not the shell: `cargo test --locked --workspace` takes long
enough here to be **moved to the background**, and its output file is a transcript that is still
growing — a `grep` over it for `FAILED` finds nothing, because the failing binary has not run yet.
Reading a partial file is the same mistake as reading a roll-up. Wait for the exit code, and read
*that*. The claim that a test suite passed is the one claim in this repository that has never
survived being inferred.

**The tenth was not the shell, the tool boundary or the tree — it was the compiler.** Nothing here
pins a toolchain and CI takes `dtolnay/rust-toolchain@stable`, so **stable is the contract** and this
gate was two releases behind it: clippy **1.96** locally against CI's **1.98**. Two lints that do not
exist in 1.96 were failing every job that runs clippy, and both workspaces' gates were green.
`rustup update` before trusting a green. A gate on an older toolchain than CI's is a gate reporting
about a different repository, and it fails in the direction that looks like success.

```sh
# Correct in a file, inert when pasted -- see above.
set -euo pipefail

cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
# CI runs this as `test (release, LTO)` and CONTRIBUTING.md has always listed it. This block did
# not, for long enough that a gate built from it ran debug only and said "the gate passed".
cargo test --locked --workspace --release
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check
for e in beam_hot_spot airy_pattern detector_snr room_modes melting lens_spots \
         heat_in_three_dimensions room_in_three_dimensions busbar_rating optical_bench \
         espresso_shot portafilter_flow agents_quickstart readme_check; do
  cargo run --locked --release --example "$e"
done
# MSRV. No `--exclude` any more: everything in this workspace is published and every one of
# it is held to 1.78. The application lives in `app/`, which cargo keeps out.
cargo +1.78 build --locked --workspace

echo "the gate passed"     # and if this line does not appear, it did not
```

**The two wasm steps belong in this block, and did not have them.** A test that reads the
repository off a disk needs `#![cfg(not(target_family = "wasm"))]`, because a `wasm32` target has
none — `counts_in_prose` and `citation_is_valid` both carry that line. A new one shipped without
it, passed twenty-seven checks here and turned CI's `test (wasm32-wasip1, wasmtime)` job red. Add:

```sh
cargo test --locked --workspace --target wasm32-wasip1 --no-run
cargo build --locked --workspace --target wasm32-unknown-unknown
```

`--no-run` because *running* needs wasmtime, which the CI job installs and this machine does not.
That is the half CI keeps; compiling is the half that catches the missing `cfg`, and it is the
half that was missing here.

CI does **not** cover `bindings/python` or `app/` from this gate — each has its own job and its own
procedure.

**`app/` has its own gate and it is not optional.** Everything a person *runs* lives there — the
CLI, the viewer, the editor and the GPU accelerator — and so do the thirty scenes and their
closed-form checks, which used to run in the line above. Two workspaces, two gates:

```sh
cd app
cargo fmt --all --check
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo test --locked --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --locked --workspace --no-deps
cargo deny check
cargo build --locked -p editor-wasm --target wasm32-unknown-unknown
```

The device tests there skip loudly on a machine with no adapter. A skip that says why is a
result; a skip that says nothing is the shape of a suite that has stopped testing anything.

**Read a result from the thing that produced it, and read whether the check *ran* rather than what it
printed.** Six times this gate has said `ok` while failing, and once that reached `main`. A CI run's
roll-up has said `success` with a job still `queued`; ask each job for its own `conclusion`. A script
that edits several files can write the first and raise on the next, so check every anchor before
writing any of them.

## The five conventions worth knowing before you start

Stated in full in CONTRIBUTING.md. The compressed version:

1. **Check against a closed form, never against another implementation.** An exact limit, a
   conservation law, or a convergence rate. Comparing to a second copy of the same idea
   verifies nothing.
2. **A tolerance has to be earned.** It should trace to an effect — an integrator's order,
   `1/√N`, a discretisation. A tolerance that was loosened in the same change that made a test
   pass is the thing to look for.
3. **Nothing is random and nothing consults a clock.** `Rng::for_index(seed, index)`. Results
   are bit-for-bit identical across platforms, optimisation levels, WebAssembly and thread
   counts, and there is a pinned digest that says so.
4. **The kernel must never depend on a domain**, and no domain may depend on another. That
   claim is the reason for the crate split and has now been held through eleven domains.
5. **Every public item is documented.** `#![deny(missing_docs)]` in all **twenty-two** crates: the
   seventeen in `crates/` and the five libraries in `app/`. It said eighteen through the
   consolidation and one crate had never had the rule at all.

## Commit messages

Say what was **wrong**, not only what changed, and give the measurement rather than the adjective.
Several of the more useful commits here are corrections to a mistaken assumption, and the reasoning
is the part worth keeping. Backticks in a message break `git commit -m` under some shells — write the
message to a file and use `git commit -F`.

**Commit as `ypark.dev@gmail.com`.** This is a personal project and every commit says so. Set it per
repository with `git config user.email ypark.dev@gmail.com`, which does not touch any global identity —
and a fresh clone starts from the global one, so set it again before the first commit there.

The whole history was rewritten to that address on 2026-08-12: all 150 commits carry it as author and
committer, in the message bodies as well as the address fields. The first pass changed only the fields
and left two messages quoting the old address in prose, which is the kind of thing to check for
directly rather than to assume a rewrite covered.

**No `Co-Authored-By` trailer.** An agent's default is to add one and this repository does not want it:
the commits are the maintainer's, and a tool that wrote some of the text is not a contributor. The
seventy that carried one no longer do — they came off in the same rewrite that changed the author
address, which is the only reason the price was worth paying: one force-push and one set of twelve
recreated tags bought both. This line is the standing instruction so the next session does not put it
back.

## The subagent team

`.claude/agents/` holds eight reviewers, each built from a defect this repository actually shipped
rather than from a generic role. They are for developing pantometry and are useless to a consumer.

| agent | asks |
| --- | --- |
| `physics-checker` | Is this claim true? Finds an independent check for it |
| `numerics-reviewer` | Would this test *notice* if the code were wrong? |
| `silent-failure-hunter` | What here can come out *empty* and look fine? |
| `unearned-pass-hunter` | Did the check *run*, and against the thing it claims to test? |
| `consumer-advocate` | What is this API like from outside? |
| `invariant-guard` | Kernel purity, cross-domain deps, determinism, docs, licence, MSRV |
| `domain-builder` | Scaffolds a new physics crate on the kernel |
| `prose-auditor` | Do the numbers in the README and doc comments still match the code? |

Run `numerics-reviewer` and `physics-checker` on anything touching physics or tolerances, and
`invariant-guard` before a commit that adds a crate, a dependency, or anything with randomness in it.
`prose-auditor` before a release. **`unearned-pass-hunter` on anything that adds a check** — a test,
a harness, a gate step, a CI job, a script that verifies something. The table of disguises above is
its subject, and four of the nine instances in the session that produced it were caught by the gate
rather than by anybody reading the output.

**Verify what they report by reproducing it.** Their findings have been wrong in both directions: a
seed-spread reported as 0.66% measured 0.96%, which would have produced a tolerance that looked
earned and was not.

## Where the rest of it is, and when to go there

| read | before |
| --- | --- |
| [RELEASING.md](RELEASING.md) | any release. Cadence, the **nine** places a version lives, the crate order, the wheel, what the pipeline has actually been run through — and what it costs to change the project's name, which is not a rename because a published name is permanent |
| [EVIDENCE.md](EVIDENCE.md) | claiming that something is checked. Ten sections of what the closed forms are and what checking against them found, including two defects the conservation audit called clean |
| [EXAMPLES.md](EXAMPLES.md) | adding or changing an example. What each one demonstrates and what it checks itself against |
| [CONTRIBUTING.md](CONTRIBUTING.md) | changing a test or a tolerance. The authority on the gate and on the five conventions in full |
| [app/pantometry-world/FRICTION.md](app/pantometry-world/FRICTION.md) | changing the public API. Thirty-four findings from using the SDK as a stranger, five of them the same underlying decision |
| [bindings/python/README.md](bindings/python/README.md) | touching the bindings. Its own cargo workspace, its own gate, and the two boundary decisions not to relitigate |
| [app/README.md](app/README.md) | touching anything a person *runs*. One workspace, one binary: the CLI, the viewer, the editor and the GPU accelerator. What the merge bought, what it cost, and the gate it has of its own |
| [app/viewer-core/README.md](app/viewer-core/README.md) | touching the viewer. Why it does not link `pantometry`, and the test that holds that now the workspace boundary does not |
| [app/editor-core/README.md](app/editor-core/README.md) | touching the editor. Why it *does* link `pantometry`, the shaded viewport, and the two halves the platform rules keep apart |
| [tools/screenshot/README.md](tools/screenshot/README.md) | changing the editor's interface. `docs/editor.png` is the one figure no command in CI can refresh, so it carries `docs/editor.txt` — the same frame through `--ui-dump` — and a test that fails when the picture is stale |
| [tools/presets/README.md](tools/presets/README.md) | adding or removing a scene, or touching the New-project screen. It writes `presets.rs` and the 27 tiles the chooser draws, both of which are **committed** — the exception to "nothing generated is committed", and what holds them against the scenes |
| [tools/report-check/README.md](tools/report-check/README.md) | touching the HTML report's viewer. It is four hundred lines of JavaScript in a Rust string and this is the only thing that executes it — plus why a `vm.runInContext` harness measured a renderer 30x slower than it is |
| [.claude/agents/README.md](.claude/agents/README.md) | adding a reviewer |

Two of those are separate workspaces because a dependency tree does not belong in the library's
lockfile. Measured: the library resolves **12** external crates, `bindings/python` **15**, a wgpu
window **86**, and a GUI shell **371**. `deny.toml` gates every one of the library's twelve, CI
builds with `--locked`, and the same crates go to `wasm32` and Rust 1.78 — none of which can carry
a GPU stack, a libpython link, or a window toolkit.

It was **four** workspaces above the library and is now one. That argument is about the boundary
between the library and everything above it, and it was never an argument for three boundaries
*above* it — the viewer, the editor, the CLI and the accelerator are all things you do to a run and
share the camera, the colour scale and the scene format. `app/README.md` has what the merge bought
and what it cost.

`app/pantometry-world` is the scene format and the first consumer: an application, `publish = false`,
whose purpose is to use the SDK the way a stranger would and report back. Anything you find yourself
adding to it that a *consumer* would want is in the wrong crate.

## What is deliberately not here

No GPU in the library, no implicit solvers, no mesh generation, no unstructured grids, no FEM beyond
the trilinear element `pantometry-elastic` uses. Adding physics means a new crate on the kernel, not a new
branch inside an existing one. If a change would make the kernel know about a domain, stop and
reconsider — see `domain-builder`, which is instructed to stop and report in exactly that case.
