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

**A directory is looked for as a path, not as a word**: quoted or slashed, as `"app"`, `"app/`,
`/app"` or `/app/`. `app/README.md` sits in `app`, and as a word those letters are in 143 of the 196
test targets (`apply`, `approximately`), which would be the whole suite run over a README. As a
path they are in four targets and one module. Over every directory a tracked document sat in at
`f796d37`, what that left out was comments and messages that mention one, never a read. The
patterns reach `git grep` on its standard input. Git Bash rewrites an argument that starts with
`/` into a Windows path, so `-e /app/` matches nothing, and a file of patterns named by its path is
not found at all once `MSYS_NO_PATHCONV=1` turns that rewrite off.

**A file under `tests/` is not always a test.** Cargo's targets are `tests/NAME.rs` and
`tests/NAME/main.rs`. Anything else there is a module, compiled into each target that declares it
with `mod NAME;`, and a hit in one runs those targets. A module no target declares is refused.
Resolving by directory holds only while modules are declared the plain way, so a module's hit is
also refused while any test uses `#[path]`, `cfg_attr(..., path = ...)` or `include!`, while the
crate's tests declare that module inline as `mod NAME { ... }`, or while its manifest has `[[test]]`
or `autotests`. Each of those can compile a file into a target that names neither the file nor its
directory. The tree had none of them on 2026-11-14, so these refuse nothing today.

It over-includes. A change to `ci.yml` runs `every_scene_can_be_drawn` and
`every_scene_is_offered_as_a_start`, whose comments mention the file. This named the first alone,
and the second's comment was there when it was written. That costs minutes and is the right
direction for the error to be in: the failure this gate must not have is passing over the one test
that reads the sentence.

Each package's tests run in one `cargo test`, and the gate reads the exit code from that run. It
also requires every named target to report tests that ran and passed, because a target whose
tests were all ignored prints `running 0 tests` and exits 0.

`PROSE_GATE_LIST=1` prints what it would run and runs nothing, and exits 4, because a list is not a
pass. An error the script did not expect prints the line it stopped on and exits 3. The rest are
0 for a pass, 1 for a failed test and 2 for a refusal.

## What "passed" does not cover, said before it

A `.md` file is read by the always-run walks, and `CITATION.cff` by `citation_is_valid`. Every
other changed path is listed before the last line, under `what reads these, as far as a search of
the tests can say:`. Each comes with every line of a test that has its whole path outside a comment,
or with `no test has its path outside a comment`, and the listing says where rather than deciding
what those lines check. For `ci.yml` it is the line of `every_example_is_run_and_only_examples_are`
that reads its example loop, and no test reads any other line of it. A workflow is parsed as YAML
first and fails the gate if it does not parse. A workflow that parses and deletes a job still gets
through, and the listing is how you know how little of it was checked.

**The listing errs toward saying less was checked**, the other way from the search that chooses
the tests. It looks for the path, not the name, so another file's `config.yml` is not taken for
this one's, and a test that joins the name to a directory or walks the directory is run and not
listed. A line that starts with `//`, `/*` or `*` is a comment. Three things still count as code: a
comment after code on the same line, a line inside a block comment that does not start with `*`,
and the path in a string such as an assertion's message. It asks git for no column, because
`grep.column` in a git config puts one before the text.

**This section said no test read anything else.** It said `ci.yml` was "named in two, both times in
a comment", and the script printed "no test reads what it says". All of it was true on 2026-09-29,
when it was written. On 2026-10-02 `7f1c3c5` made `every_example_is_run_and_only_examples_are`
read `ci.yml`'s example loop, and all of it stood until 2026-11-17. So the listing is found by a
search, not written down.

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

### And the version that resolves modules, 2026-11-14

The first version took every file under `tests/` for a target, because in a git pathspec `*`
crosses a `/`. Given a newline at the end of `app/README.md`, it printed one line from `sed`,
`can't read crates/pantometry-forcefield/tests/Cargo.toml`, and exited 2, the refusal's code. The
word `app` had matched eight of the twelve modules under `tests/*/`, `tests/pinned/mod.rs` among
them, and the gate went looking for a manifest beside one. Fixing only that would have run 143
test targets over a README.

Run against the version above on 2026-11-14, each change restored by bytes. The ones that need a
committed `.rs` change ran in a throwaway worktree, listed and not run, and each guard's case was
also run through the version without that guard:

| change | result |
| --- | --- |
| a newline added to `app/README.md` | eleven targets in five packages, one of them `benzene_bound_in_tip3p`, because the `amber` module it declares quotes an `openmm/app/data` path; passed in 164 s, listed alone in 39 s |
| a module, `tests/pinned/mod.rs`, naming a changed file | the one target that declares `mod pinned;` |
| a module no target declares, naming it | refused: no target declares it |
| a directory target, `tests/NAME/main.rs`, naming it | target `NAME` |
| `tests/helpers/reads.rs` naming it, declared by `mod helpers;` in one target and compiled into another by `include!`, by `cfg_attr(all(), path = ...)`, by an inline `mod helpers { pub mod reads; }` or by a `[[test]]` | refused, each; without the guard, only the first target was listed |
| a module naming it while a `#[path]` sits in another crate's tests, or in a target that pulls a file out of a directory target | refused, each; without the guard, one target was listed |
| `git` not on `PATH` | stopped, with the line, exit 3 |
| a comment appended to `ci.yml` | listed, with the `parses as YAML` note; not run |
| nothing changed; a bindings README; a `.rs` file; an unclosed list in `ci.yml`; a count staged and then changed; a wrong count | the result in the table above, each |

The first guard on `cfg_attr` matched `cfg_attr\([^)]*path`, which stops at the `)` of `all()`: the
`cfg_attr` case came back listed rather than refused, and that is how it was found.

`unearned-pass-hunter` read this version. It found the six ways a module could be compiled into a
target that names neither it nor its directory, and the file of patterns that `MSYS_NO_PATHCONV=1`
made unreadable. It also found exit codes outside the scheme above, and counts in this section
that the measurements did not support. All of them are fixed above.

### And the version that says where a file is named, 2026-11-17

Run on 2026-11-17 by a script that held the new `run.sh` outside the tree, so the tree held only the
change each case made. Each change was restored by bytes and checked by SHA-256. The six that need
a committed `.rs` change ran in a throwaway worktree, listed and not run. Each committed one line
or a few to the end of `counts_in_prose.rs`, or the one change named, and then appended a comment
to the workflow:

| change | result |
| --- | --- |
| a comment appended to `ci.yml`, through the version before | listed, with `no test reads what it says` |
| a comment appended to `ci.yml` | passed in seven targets, with `its path is in a test outside a comment at crates/pantometry/tests/counts_in_prose.rs:1132`; the two comments that name the file are not listed |
| an example dropped from `ci.yml`'s loop | `every_example_is_run_and_only_examples_are` failed; the gate FAILED |
| a comment appended to `release-python.yml` | listed, with `no test has its path outside a comment` |
| an unclosed list appended to `ci.yml` | failed: does not parse |
| `ci.yml`, with line 1132 of `counts_in_prose.rs` made a comment | `no test has its path outside a comment` |
| `release-python.yml`, with `"app/tool/release-python.yml"` in code | `no test has its path outside a comment` |
| `release-python.yml`, with its path in a `/* */` comment and on a ` * ` line | `no test has its path outside a comment` |
| `release-python.yml`, with its path in a `//` comment and `grep.column=true` | `no test has its path outside a comment` |
| the same, through a copy without `--no-column` | that comment listed as where its path is: the case reaches the guard |
| `release-python.yml`, with `// .github/workflows/release-python.yml` after a line of code | that line listed: the stated limit, a comment after code taken for code |

`unearned-pass-hunter` read the first draft, which searched for the name. It found that
`grep.column` put a column before the text and turned both of `ci.yml`'s comments into reads; that
a search built to over-include could name a read that was not there, through another file's
`config.yml` or a `/* */` comment; and that this README wrote `1132` down in the change whose
point is that it is found. The path, the wider comment rule and `--no-column` are the fixes. A path
in a string is still taken for a read, which the section above says, and so is a `//!` behind a
byte-order mark, which no test file has.

What it does not do: build anything, run clippy, or check formatting. There is no Rust in the
change, which is the condition it checks first.
