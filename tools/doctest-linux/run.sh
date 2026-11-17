#!/bin/bash
# app/'s doctests on Linux, from a Windows machine, in WSL. See README.md beside this file.
#
#   bash ../tools/doctest-linux/run.sh      # from app/, as the gate runs it; any directory works
#
# Exit 0 only when cargo exited 0, no doctest failed and at least one passed. 1 for a failed or an
# empty run, 3 when it could not run at all, and cargo's or wsl.exe's own code otherwise.
main() {
  if [ "${1-}" = "--inside" ]; then
    shift
    inside "$@"
    exit
  fi
  cd "$(dirname "$0")/../../app" || { echo "REFUSED: no app/ beside tools/"; exit 3; }
  case "$(uname -s)" in
    MINGW* | MSYS* | CYGWIN*) windows ;;
    *)
      # WSL is Linux, and a shell in it that finds this line has skipped the host's toolchain.
      if [ -n "${WSL_DISTRO_NAME-}" ]; then
        echo "REFUSED: this is WSL called as a Linux machine. Run this from Git Bash on Windows,"
        echo "which runs it in here with the host's toolchain."
        exit 3
      fi
      echo "not Windows: running app/'s doctests here, on $(uname -s), where cfg(windows) is false"
      judged cargo test --locked --workspace --doc
      ;;
  esac
}

# Runs "$@" with its output streamed and kept, and holds it to the result. `--doc` runs nothing but
# doctests, so every result line is one of theirs, and over no doctests it prints `0 passed` for
# every crate and exits 0.
judged() {
  local out code passed failed
  out=$(mktemp) || exit 3
  "$@" 2>&1 | tee "$out"
  code=${PIPESTATUS[0]}
  passed=$(grep -cE '^test result: ok\. [1-9][0-9]* passed' "$out")
  failed=$(grep -c '^test result: FAILED' "$out")
  rm -f "$out"
  echo "doctest-linux: exit $code; $passed result lines with a doctest passed, $failed failed"
  [ "$code" -eq 0 ] || exit "$code"
  [ "$failed" -eq 0 ] || exit 1
  [ "$passed" -ge 1 ] || { echo "no doctest passed, so nothing was checked"; exit 1; }
}

windows() {
  command -v wsl.exe > /dev/null || { echo "REFUSED: no wsl.exe, so nothing here can run Linux"; exit 3; }
  local rustc cargo app self
  rustc=$(rustc -V) || { echo "REFUSED: no rustc on this machine"; exit 3; }
  cargo=$(cargo -V) || { echo "REFUSED: no cargo on this machine"; exit 3; }
  app=$(pwd -W) || exit 3
  self=$(cd ../tools/doctest-linux && pwd -W)/run.sh
  # There is no time limit here, because the gate's timeout is enough: it kills wsl.exe, and WSL
  # ends what that wsl.exe started. Measured on 2026-11-17: a 41 s sleep behind this pipe and a
  # 43 s one alone were gone at the first look, and a cargo cut off mid-build was gone within
  # 0.3 s, with the next run compiling what it had not finished.
  export MSYS_NO_PATHCONV=1
  judged wsl.exe ${PANTOMETRY_WSL:+-d "$PANTOMETRY_WSL"} --exec /bin/bash -c \
    'exec /bin/bash "$(wslpath -u "$1")" --inside "${@:2}"' _ "$self" "$rustc; $cargo" "$app"
}

# The Linux half. It enters app/ before it asks for the toolchain, so that an override set for
# that directory is the one compared.
inside() {
  local host=$1 app=$2 dir rustc cargo
  dir=$(wslpath -u "$app") && cd "$dir" || { echo "REFUSED: cannot enter $app from WSL"; exit 3; }
  export PATH=$HOME/.cargo/bin:$PATH
  rustc=$(rustc -V) || { echo "REFUSED: no rustc in WSL"; exit 3; }
  cargo=$(cargo -V) || { echo "REFUSED: no cargo in WSL"; exit 3; }
  if [ "$rustc; $cargo" != "$host" ]; then
    echo "REFUSED: the toolchains differ, Windows '$host' and WSL '$rustc; $cargo'."
    echo "CI takes stable, so run 'rustup update stable' on whichever is behind."
    exit 3
  fi
  # Not app/target, which Windows builds in, because two hosts in one directory would each take the
  # other's build scripts and fingerprints for their own. And one directory per tree, because
  # cargo names a unit by its path relative to the workspace, so a second checkout could reuse
  # the first one's build of a crate whose files are older than that build.
  export CARGO_TARGET_DIR
  CARGO_TARGET_DIR=${XDG_CACHE_HOME:-$HOME/.cache}/pantometry/doctest-linux/$(printf %s "$dir" | sha256sum | cut -c1-12)
  echo "WSL: $(. /etc/os-release && echo "$PRETTY_NAME"), $(uname -r); $rustc; $dir; $CARGO_TARGET_DIR"
  # CI's `the app` job sets this for the step whose `cargo test` runs these doctests.
  export CARGO_PROFILE_DEV_OPT_LEVEL=2
  exec cargo test --locked --workspace --doc
}

main "$@"
exit
