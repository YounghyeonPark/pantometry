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
set -euo pipefail

cd "$(git rev-parse --show-toplevel)"

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
mapfile -t changed < <(printf '%s\n' "$status" | sed -E 's/^.. //' | sed -E 's/ -> /\n/' | sed '/^$/d' | sort -u)
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
set +e
plain=$(git grep -nE 'include_(str|bytes)!\("[^"]*"\)' -- '*.rs'); plain_code=$?
around=$(git grep -nE -A3 'include_(str|bytes)!' -- '*.rs'); other_code=$?
set -e
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
# Found from the tests themselves, not from a list somebody keeps: every test file that names a
# changed file, or the directory it sits in. A walk over `.claude/agents` never names the agent
# that changed, so the directory counts too. Over-inclusion costs minutes; under-inclusion is a
# gate that passes over the one test that reads the sentence.
declare -A wanted
# The documentation checks run whatever changed: they hold counts, versions, links and the
# citation across every document at once, so an edit to one can break a claim in another.
always=(crates/pantometry/tests/counts_in_prose.rs crates/pantometry/tests/documented_version.rs
        crates/pantometry/tests/the_documents_link_to_things_that_exist.rs
        crates/pantometry/tests/citation_is_valid.rs app/pantometry-world/tests/friction_counts.rs)
for t in "${always[@]}"; do wanted[$t]=1; done
unread=()
for p in "${changed[@]}"; do
  terms=("$(basename "$p")")
  dir=$(dirname "$p")
  [ "$dir" != "." ] && terms+=("$dir")
  found=0
  for term in "${terms[@]}"; do
    set +e
    hits=$(git grep -lF "$term" -- 'crates/*/tests/*.rs' 'app/*/tests/*.rs'); code=$?
    set -e
    [ "$code" -le 1 ] || refuse "git grep failed looking for $term"
    while IFS= read -r t; do
      [ -n "$t" ] && { wanted[$t]=1; found=1; }
    done <<< "$hits"
  done
  # A `.md` is read by the always-run walks and `CITATION.cff` by `citation_is_valid`. Anything
  # else is said to be unread, rather than covered by the word "passed" -- *whatever* the grep
  # found, because a name in a test is not a read: `ci.yml` is named in two tests, both in
  # comments, and neither opens it.
  : "$found"
  if [[ "$p" != *.md && "$p" != CITATION.cff ]]; then unread+=("$p"); fi
done

# A workflow no test reads can at least be read as YAML. It is still not *checked*: a valid file
# that deletes a job parses fine, which the line before `passed` says.
notes=()
for p in "${unread[@]}"; do
  if [[ "$p" == .github/*.yml || "$p" == .github/*.yaml ]] && [ -f "$p" ]; then
    python -c "import sys, yaml; yaml.safe_load(open(sys.argv[1], encoding='utf-8'))" "$p" \
      || { echo "$p does not parse as YAML"; exit 1; }
    notes+=("$p: parses as YAML; no test reads what it says")
  else
    notes+=("$p: no test reads it")
  fi
done

# --- Run them, one cargo invocation per package, every exit code read. --------------------------
declare -A groups count
for t in "${!wanted[@]}"; do
  [ -f "$t" ] || { echo "a test file this gate names is gone: $t"; exit 1; }
  crate_dir=$(dirname "$(dirname "$t")")
  ws=.; [[ "$t" == app/* ]] && ws=app
  pkg=$(sed -nE 's/^name = "([^"]+)"/\1/p' "$crate_dir/Cargo.toml" | head -1)
  [ -n "$pkg" ] || { echo "no package name in $crate_dir/Cargo.toml"; exit 1; }
  groups["$ws|$pkg"]+=" --test $(basename "$t" .rs)"
  count["$ws|$pkg"]=$(( ${count["$ws|$pkg"]:-0} + 1 ))
done

log=$(mktemp)
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
rm -f "$log"

[ "$failed" -eq 0 ] || { echo "the prose gate FAILED"; exit 1; }
[ ${#notes[@]} -eq 0 ] || { echo "not checked by any test:"; printf '  %s\n' "${notes[@]}"; }
echo "the prose gate passed"
