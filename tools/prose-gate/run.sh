#!/bin/bash
# The gate for a change that touches no code: prose, the agent team, CI's own files, the citation.
#
# It decides for itself whether it applies, and runs every test that reads what changed. Run it
# from the repository root, as a file:
#
#   bash tools/prose-gate/run.sh
#
# The last line says `the prose gate passed`, or it did not pass. Anything it refuses goes through
# the full gate in CLAUDE.md, which is the default and always correct.
#
# Exit codes: 0 passed, 1 a test failed, 2 refused, 3 this script itself went wrong, and 4 for
# `PROSE_GATE_LIST=1`, which prints the tests it would run and runs none -- which is not a pass.
set -Eeuo pipefail

# **A failure in this script is not a refusal.** A `sed` over a path that was not there once ended
# the run here with exit 2, the refusal's code, and printed nothing at all. Anything this script
# did not expect now says where it stopped, under a code of its own.
trap 'echo "the prose gate stopped on an error it did not expect, at line $LINENO: $BASH_COMMAND"; exit 3' ERR

top=$(git rev-parse --show-toplevel)
cd "$top"

refuse() {
  echo "this change is not for the prose gate; run the full gate in CLAUDE.md:"
  printf '  %s\n' "$@"
  exit 2
}

# --- Does it apply? --------------------------------------------------------------------------
# Every changed path, staged or not, tracked or not; a rename counts as both of its names. A path
# git has to quote (a space, non-ASCII) keeps its quotes here and so matches nothing below: it is
# refused, which is the safe way to be wrong.
status=$(git status --porcelain=v1 -uall)
# A command substitution and not `< <(...)`, whose failure nothing would see.
paths=$(printf '%s\n' "$status" | sed -E 's/^.. //' | sed -E 's/ -> /\n/' | sed '/^$/d' | sort -u)
changed=()
[ -z "$paths" ] || mapfile -t changed <<< "$paths"
if [ ${#changed[@]} -eq 0 ]; then
  echo "nothing has changed, so there is nothing for this gate to say"
  exit 2
fi

# **What is tested has to be what is committed.** The tests read the working tree and `git commit`
# takes the index, so a file staged with one text and then edited again on disk (`MM`) would be
# checked as one thing and committed as the other -- a wrong count staged, "fixed" only on disk,
# passes `counts_in_prose`. Measured by `unearned-pass-hunter` in a scratch repository.
mapfile -t twice < <(printf '%s\n' "$status" | grep -E '^[MARC][MD] ' | sed -E 's/^.. //' || true)
[ ${#twice[@]} -eq 0 ] || refuse "${twice[@]/%/ (staged, then changed again: stage all of it or none)}"

refused=()
for p in "${changed[@]}"; do
  # `bindings/` has its own workspace and gate.
  if [[ "$p" == bindings/* ]]; then refused+=("$p (bindings have their own gate)"); continue; fi
  if [[ "$p" == *.md || "$p" == .claude/* || "$p" == .github/* || "$p" == CITATION.cff ]]; then
    continue
  fi
  refused+=("$p")
done
[ ${#refused[@]} -eq 0 ] || refuse "${refused[@]}"

# **A document pulled into a build is code**, whatever its extension -- `bindings/python/README.md`
# is the Python module's `__doc__`. Asked of tracked Rust with `git grep`, whose status is read
# rather than lost in a pipe: a `grep -r` that hit an unreadable file under `pipefail` switched this
# check off without a word.
#
# A plain `include_str!("path")` is resolved against the file it is in and compared with the
# changed paths, so editing one README is not refused because another crate includes a different
# README. Any other shape -- `concat!(env!(...), ...)`, or rustfmt's call split across lines -- is
# refused if a changed file's name appears on its line or the three after, since this cannot
# resolve it and must not guess.
#
# Each status is taken with `|| code=$?`, not inside `set +e`: the trap above fires on a failure
# that `set +e` lets through, and `git grep` finding nothing is exit 1.
plain_code=0
plain=$(git grep -nE 'include_(str|bytes)!\("[^"]*"\)' -- '*.rs') || plain_code=$?
other_code=0
around=$(git grep -nE -A3 'include_(str|bytes)!' -- '*.rs') || other_code=$?
# The plain includes are resolved above; what is left is every other shape, with its next lines.
other=$(printf '%s\n' "$around" | grep -vE 'include_(str|bytes)!\("[^"]*"\)' || true)
[ "$plain_code" -le 1 ] && [ "$other_code" -le 1 ] || refuse "git grep failed ($plain_code, $other_code); cannot say what is compiled in"
while IFS= read -r hit; do
  [ -n "$hit" ] || continue
  file=${hit%%:*}
  lit=$(printf '%s' "$hit" | sed -nE 's/.*include_(str|bytes)!\("([^"]*)"\).*/\2/p')
  target=$(realpath -m --relative-to=. "$(dirname "$file")/$lit")
  for p in "${changed[@]}"; do
    [ "$target" = "$p" ] && refused+=("$p (compiled in by $file)")
  done
done <<< "$plain"
for p in "${changed[@]}"; do
  base=$(basename "$p")
  if printf '%s\n' "$other" | grep -qF "$base"; then
    refused+=("$p (named near an include this cannot resolve)")
  fi
done
[ ${#refused[@]} -eq 0 ] || refuse "${refused[@]}"

# --- Which tests read what changed? ------------------------------------------------------------
# Found from the tests themselves, not from a list somebody keeps: every file under a crate's
# `tests/` that names a changed file, or the directory it sits in. A walk over `.claude/agents`
# never names the agent that changed, so the directory counts too. Over-inclusion costs minutes;
# under-inclusion is a gate that passes over the one test that reads the sentence.
#
# **A file under `tests/` is not always a test.** Cargo's targets are `tests/NAME.rs` and
# `tests/NAME/main.rs`; anything else there is a module, compiled into whichever targets declare
# it. This used to take every hit for a target -- `crates/*/tests/*.rs` lets `*` cross a `/` -- and
# went looking for `tests/Cargo.toml` beside `tests/pinned/mod.rs`. So a hit is resolved to the
# targets that compile it, and one that resolves to none is refused rather than guessed at.
#
# Resolving by directory holds only while modules are declared the plain way. A `#[path]`, a
# `cfg_attr(..., path = ...)` or an `include!` anywhere in the tests can compile a file into a
# target that names neither it nor its directory; `mod NAME { ... }` inline can pull in
# `tests/NAME/` with no `mod NAME;` to find; and `[[test]]` or `autotests` in a manifest makes
# targets of its own. Any of those and a module's hit is refused. The tree has none of them at
# 2026-11-14, so this refuses nothing today and is what keeps "resolved" from becoming a guess.
tests=(':(glob)crates/*/tests/**/*.rs' ':(glob)app/*/tests/**/*.rs')
opaque_code=0
# `cfg_attr` is matched to the end of its line: `cfg_attr\([^)]*path` stopped at the `)` of
# `cfg_attr(all(), path = ...)` and let it through, which is how the case that checks this failed.
opaque=$(git grep -lE '#[[:space:]]*\[[[:space:]]*path|cfg_attr.*path[[:space:]]*=|(^|[^_a-z])include![[:space:]]*\(' -- "${tests[@]}") || opaque_code=$?
[ "$opaque_code" -le 1 ] || refuse "git grep failed ($opaque_code) looking for #[path] and include!"
resolved=()
why=
resolve() { # a file under some crate's tests/ -> `resolved` holds "crate_dir|target", or `why` says
  local f=$1 crate_dir rest sub decl inline code t
  resolved=(); why=
  crate_dir=${f%%/tests/*}; rest=${f#*/tests/}
  if [[ "$rest" != */* ]]; then resolved=("$crate_dir|${rest%.rs}"); return 0; fi
  sub=${rest%%/*}
  if [ -n "$opaque" ]; then
    why="${opaque%%$'\n'*} compiles a file by a path of its own (#[path], a cfg_attr path or include!), so no module's directory says which targets compile it"
    return 0
  fi
  if grep -qE '^[[:space:]]*(\[\[test\]\]|autotests)' "$crate_dir/Cargo.toml"; then
    why="$crate_dir/Cargo.toml declares test targets of its own"
    return 0
  fi
  if [ -f "$crate_dir/tests/$sub/main.rs" ]; then resolved=("$crate_dir|$sub"); return 0; fi
  code=0
  inline=$(git grep -lE "(^|[^A-Za-z0-9_])mod[[:space:]]+$sub[[:space:]]*\{" -- ":(glob)$crate_dir/tests/*.rs") || code=$?
  [ "$code" -le 1 ] || { why="git grep failed ($code) looking for mod $sub {"; return 0; }
  if [ -n "$inline" ]; then
    why="${inline%%$'\n'*} declares mod $sub inline, and what that pulls in from tests/$sub/ cannot be read off the directory"
    return 0
  fi
  code=0
  decl=$(git grep -lE "(^|[^A-Za-z0-9_])mod[[:space:]]+$sub[[:space:]]*;" -- ":(glob)$crate_dir/tests/*.rs") || code=$?
  [ "$code" -le 1 ] || { why="git grep failed ($code) looking for mod $sub"; return 0; }
  while IFS= read -r t; do
    [ -n "$t" ] && resolved+=("$crate_dir|$(basename "$t" .rs)")
  done <<< "$decl"
  [ ${#resolved[@]} -gt 0 ] || why="no target in $crate_dir/tests declares mod $sub"
  return 0
}

declare -A wanted
# The documentation checks run whatever changed: they hold counts, versions, links and the
# citation across every document at once, so an edit to one can break a claim in another.
always=(crates/pantometry/tests/counts_in_prose.rs crates/pantometry/tests/documented_version.rs
        crates/pantometry/tests/the_documents_link_to_things_that_exist.rs
        crates/pantometry/tests/citation_is_valid.rs app/pantometry-world/tests/friction_counts.rs)
for t in "${always[@]}"; do
  [ -f "$t" ] || { echo "a test file this gate names is gone: $t"; exit 3; }
  resolve "$t"
  [ ${#resolved[@]} -eq 1 ] || { echo "$t is not a test target of its own: $why"; exit 3; }
  wanted[${resolved[0]}]=1
done
# **A directory is looked for as a path, not as a word.** `app/README.md` sits in `app`, and those
# three letters are in `apply` and `approximately`: 143 of the 196 test targets, which is the whole
# suite run over a README. Quoted or slashed -- `"app"`, `"app/`, `/app"`, `/app/` -- they are in 4
# targets and 1 module, and over every directory a tracked document sat in at f796d37, what that
# left out was comments and messages that mention one, never a read. The patterns reach `git grep`
# on its standard input: Git Bash rewrites an argument that starts with `/` into a Windows path, so
# `-e /app/` matches nothing, and a file of patterns named by its path is not found at all once
# `MSYS_NO_PATHCONV=1` turns that rewrite off.
unread=()
for p in "${changed[@]}"; do
  base=$(basename "$p")
  dir=$(dirname "$p")
  code=0
  hits=$(git grep -lF -e "$base" -- "${tests[@]}") || code=$?
  [ "$code" -le 1 ] || refuse "git grep failed looking for $base"
  if [ "$dir" != "." ]; then
    code=0
    more=$(printf '"%s"\n"%s/\n/%s"\n/%s/\n' "$dir" "$dir" "$dir" "$dir" | git grep -lF -f - -- "${tests[@]}") || code=$?
    [ "$code" -le 1 ] || refuse "git grep failed looking for $dir as a path"
    hits+=$'\n'"$more"
  fi
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    resolve "$f"
    [ ${#resolved[@]} -gt 0 ] || refuse "$f names $p, and $why"
    for k in "${resolved[@]}"; do wanted[$k]=1; done
  done <<< "$hits"
  # A `.md` is read by the always-run walks and `CITATION.cff` by `citation_is_valid`. Anything
  # else is said to be unread, rather than covered by the word "passed" -- *whatever* the grep
  # found, because a name in a test is not a read: `ci.yml` is named in two tests, both in
  # comments, and neither opens it.
  if [[ "$p" != *.md && "$p" != CITATION.cff ]]; then unread+=("$p"); fi
done

# A workflow no test reads can at least be read as YAML. It is still not *checked*: a valid file
# that deletes a job parses fine, which the line before `passed` says.
notes=()
for p in "${unread[@]}"; do
  if [[ "$p" == .github/*.yml || "$p" == .github/*.yaml ]] && [ -f "$p" ]; then
    # A missing parser is this script's problem, not the file's: it would otherwise be reported as
    # a workflow that "does not parse".
    python -c "import yaml" 2>/dev/null || { echo "cannot read $p as YAML: no python with PyYAML here"; exit 3; }
    python -c "import sys, yaml; yaml.safe_load(open(sys.argv[1], encoding='utf-8'))" "$p" \
      || { echo "$p does not parse as YAML"; exit 1; }
    notes+=("$p: parses as YAML; no test reads what it says")
  else
    notes+=("$p: no test reads it")
  fi
done

# --- Run them, one cargo invocation per package, every exit code read. --------------------------
declare -A groups count
for k in $(printf '%s\n' "${!wanted[@]}" | sort); do
  crate_dir=${k%%|*}; name=${k#*|}
  [ -f "$crate_dir/Cargo.toml" ] || { echo "the prose gate found the test $name in $crate_dir, which has no Cargo.toml"; exit 3; }
  ws=.
  if [[ "$crate_dir" == app/* ]]; then ws=app; fi
  pkg=$(sed -nE 's/^name = "([^"]+)"/\1/p' "$crate_dir/Cargo.toml" | head -1)
  [ -n "$pkg" ] || { echo "no package name in $crate_dir/Cargo.toml"; exit 3; }
  groups["$ws|$pkg"]+=" --test $name"
  count["$ws|$pkg"]=$(( ${count["$ws|$pkg"]:-0} + 1 ))
done

if [ "${PROSE_GATE_LIST:-}" = 1 ]; then
  for key in $(printf '%s\n' "${!groups[@]}" | sort); do
    echo "== ${key%%|*}: -p ${key#*|}${groups[$key]}"
  done
  [ ${#notes[@]} -eq 0 ] || { echo "not checked by any test:"; printf '  %s\n' "${notes[@]}"; }
  echo "listed, and nothing run: this is not a pass"
  exit 4
fi

log=$(mktemp)
trap 'rm -f "$log"' EXIT
failed=0
for key in $(printf '%s\n' "${!groups[@]}" | sort); do
  ws=${key%%|*}; pkg=${key#*|}
  echo "== $ws: -p $pkg${groups[$key]}"
  # The exit code is cargo's own, taken from the run itself and not from a pipe's last stage.
  # shellcheck disable=SC2086
  code=0; (cd "$ws" && cargo test --locked -p "$pkg" ${groups[$key]}) > "$log" 2>&1 || code=$?
  grep -E "^test result|FAILED|panicked" "$log" || true
  # **And every target ran something.** One whose tests were all ignored or compiled out prints
  # `running 0 tests` and exits 0, which is a pass over nothing.
  ran=$(grep -cE '^test result: ok\. [1-9][0-9]* passed' "$log" || true)
  echo "   exit $code, $ran of ${count[$key]} targets ran tests and passed"
  { [ "$code" -eq 0 ] && [ "$ran" -eq "${count[$key]}" ]; } || failed=1
done

[ "$failed" -eq 0 ] || { echo "the prose gate FAILED"; exit 1; }
[ ${#notes[@]} -eq 0 ] || { echo "not checked by any test:"; printf '  %s\n' "${notes[@]}"; }
echo "the prose gate passed"
