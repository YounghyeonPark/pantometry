# `tools/prose-gate` — the gate for a change that touches no code

```sh
bash tools/prose-gate/run.sh
```

Both full gates took about 40 minutes by the time CI took 9, and on 2026-09-29 they ran twice over
commits that changed nothing but a workflow file and some prose. They
compiled the same code the previous commit had and ran the same tests against it. This gate runs
**what can have changed**, and decides for itself whether that is all that changed.

## When it applies, which it decides

Every changed path — staged, unstaged or untracked, and both names of a rename — has to be one of:

- a `.md` file outside `bindings/`
- anything under `.claude/` or `.github/`
- `CITATION.cff`

Three more refusals, each with the reason it gives:

- **A file staged and then changed again.** The tests read the working tree and `git commit` takes
  the index, so the two would be different texts: a wrong count staged and then put right only on
  disk passes `counts_in_prose` and is committed wrong. Stage all of it or none.
- **A document compiled into a build.** `bindings/python/README.md` is the Python module's
  `__doc__`. A plain `include_str!("path")` is resolved against the file it is in, so editing one
  README is not refused because another crate includes a different one. Any other shape, such as
  `concat!(env!(...))` or a call rustfmt split over lines, is refused if a changed file's name is on
  its line or the three after, because this cannot resolve it and must not guess.
- **A path git has to quote** (a space, non-ASCII). It keeps its quotes, matches nothing, and is
  refused, which is the safe way to be wrong.

Refusal is exit 2 with the list. The full gate in [CLAUDE.md](../../CLAUDE.md) is then the one to
run, and it is always correct to run it instead.

## What it runs

**Found from the tests, not from a list.** Every file under `crates/*/tests` and `app/*/tests` that
names a changed file's basename, or the directory it sits in — so a change to
`.claude/agents/prose-auditor.md` reaches the tests that walk `.claude/agents`, which never name
the agent. And always, whatever changed, the five that hold the documents against each other:
`counts_in_prose`, `documented_version`, `the_documents_link_to_things_that_exist`,
`citation_is_valid` and `friction_counts`, because an edit to one document can break a claim
written in another.

It over-includes. A change to `ci.yml` runs `every_scene_can_be_drawn`, whose comments mention the
file. That costs minutes and is the right direction for the error to be in: the failure this gate
must not have is passing over the one test that reads the sentence.

Each package's tests run in one `cargo test`, and the gate reads the exit code from that run. It
also requires every named target to report tests that ran and passed, because a target whose
tests were all ignored prints `running 0 tests` and exits 0.

## What "passed" does not cover, said before it

A `.md` file is read by the always-run walks, and `CITATION.cff` by `citation_is_valid`. **Nothing
else in a prose-only change is read by any test.** A name in a test is not a read: `ci.yml` is named
in two, both times in a comment. So every other changed path is listed under `not checked by any
test:` before the last line. A workflow is parsed as YAML first and fails the gate if it does not
parse. A workflow that parses and deletes a job still gets through, and the listing is how you know
it has not been checked.

## What it was checked against

Run against real changes, each restored by bytes afterwards:

| change | result |
| --- | --- |
| nothing changed | refused, exit 2 |
| a trailing space in CLAUDE.md | the five document tests and `threads_do_not_change_the_answer`, which names CLAUDE.md; passed |
| `twenty-four` crates written `twenty-five` in CLAUDE.md | `counts_in_prose` failed; the gate FAILED |
| that wrong count staged, then put right on disk only | refused: staged, then changed again |
| a newline added to an agent file | the tests that walk `.claude/agents`; passed |
| `52 answers` written `53` in EVIDENCE.md | found `scene.rs`, where that count is measured; failed after 264 s |
| a comment appended to `ci.yml` | passed, with `ci.yml: parses as YAML; no test reads what it says` above it |
| an unclosed list appended to `ci.yml` | failed: does not parse |
| a newline added to `app/viewer-core/README.md` | not refused; its tests ran and passed |
| a newline added to `bindings/python/README.md` | refused: bindings have their own gate |
| a newline added to `pantometry-units/src/lib.rs` | refused as not prose |

`unearned-pass-hunter` read the first version. It found the staged-then-changed hole and the
"passed" over an unread workflow. It also found three latent ones: a `grep -r` read error under
`pipefail` that switched the include check off, an include regex that missed `concat!` and split
calls, and zero-test targets counting as passes. All of them are fixed above.

What it does not do: build anything, run clippy, or check formatting. There is no Rust in the
change, which is the condition it checks first.
