# `tools/doctest-linux` — `app/`'s doctests, on Linux, from Windows

```sh
cd app
bash ../tools/doctest-linux/run.sh
```

`cargo test` in `app/` compiles the workspace's doctests for the machine it runs on, and CI's `the
app` job runs `cargo test` on `ubuntu-latest`. On Windows that leaves one mistake for CI to find
first: a doctest that uses an item a `cfg` holds to Windows, which compiles here and not there.
Nothing else in `app/`'s gate reads a doctest for Linux:

- clippy does not lint doctests;
- rustdoc compiles them only under `cargo test`, so the two `--target x86_64-unknown-linux-gnu`
  lines pass over them.

The library's doctests had the same gap and it was closed already. They also run under
`wasm32-wasip1`, where `cfg(windows)` is false too, and `app/` has no wasip1 run.

This runs `cargo test --locked --workspace --doc` in the machine's WSL distribution, on the same
working tree and at CI's opt-level. `CONTRIBUTING.md` said until 2026-11-17 that it would need a
Linux linker. This machine had one the whole time, in WSL's Ubuntu.

## What it checks, before and after

- **The toolchain is the host's.** It passes the host's `rustc -V` and `cargo -V` into WSL and
  refuses with exit 3 when WSL's are different. It asks for them after entering `app/`, so an
  override set for that directory is the one compared. CI takes `stable`, and `CLAUDE.md` records
  a gate that ran clippy 1.96 against CI's 1.98. When this was written, WSL had 1.98.1 and the host
  1.100.0. The check shows that the two sides agree. It does not show that either is stable.
- **A run in which no doctest passed fails.** Over a workspace with no doctests, `--doc` prints
  `test result: ok. 0 passed` for each crate and exits 0. Today there is one doctest, in
  `pantometry-gpu`.
- **It uses its own target directory**, `~/.cache/pantometry/doctest-linux/<tree>` in WSL. `<tree>`
  is the first twelve hex digits of the SHA-256 of the tree's path. It does not use `app/target`,
  where Windows builds, because two hosts in one directory would each take the other's build
  scripts and fingerprints for their own. It keeps one directory per tree because cargo names a
  unit by its path relative to the workspace. With one directory shared, a second checkout could
  reuse the first one's build of a crate whose files are older than that build.
- **It has no time limit of its own.** The gate's `timeout` kills `wsl.exe`, and WSL ends what that
  `wsl.exe` started; the cases below measure it. The plan was a limit inside WSL as well, on the
  belief that a cargo there would outlive its `wsl.exe`, and measured it does not.

`PANTOMETRY_WSL` names a distribution; without it, the default is used. On Linux or macOS it runs
`cargo test --locked --workspace --doc` itself and holds it to the same result. Run from a shell
inside WSL it refuses: that shell is Linux, but it skipped the comparison with the host's toolchain.

## What it was checked against

Run on 2026-11-17 in WSL's Ubuntu 26.04.1 LTS, with `rustc 1.100.0` on both sides. Each sabotage
edited `app/pantometry-gpu/src/lib.rs` by exact anchors. Each was restored by bytes and checked by
SHA-256, with a fresh mtime so that cargo would build again from the restored file.

| change | result |
| --- | --- |
| none, cold | passed in 56 s, 50 s of it cargo's; 104 crates, 1.5 GB in WSL; `OnTheGpu (line 1038) - compile ... ok` |
| none, warm | passed in 4 to 8 s |
| a `#[cfg(windows)]` function appended to the file, and a line calling it in the existing doctest | `cargo test --doc` passed here, as did clippy and rustdoc for Linux; this failed with `E0425: cannot find function held_to_windows` |
| a second doctest, `assert!(cfg!(windows))` | `2 passed` here; this failed at run time, `test result: FAILED. 1 passed; 1 failed` |
| the one doctest made `text` | cargo exited 0, and this exited 1: `no doctest passed, so nothing was checked` |
| a `rustc` on `PATH` reporting 1.99.0; then a `cargo` doing the same | each refused before cargo started, exit 3 |
| run from a shell inside WSL | refused before cargo started, exit 3 |
| run in WSL with `WSL_DISTRO_NAME` unset, WSL's toolchain on `PATH`: a Linux machine, as near as this one has | ran the doctests itself, `1 result lines with a doctest passed`, exit 0 |
| `timeout 5` around it, after editing `pantometry-world` | exit 124. At 0.1 s cargo was still in WSL; at 0.2 s nothing was. The next run compiled the three crates the killed one had not finished |

The last row was run eleven times with this version:

- five at `timeout 3`;
- four at `timeout 5`;
- two at `timeout 7`.

The two at 7 s landed after the build, so they show nothing. In the other nine, nothing of the run
was left in WSL 0.2 s after the kill, and the next run compiled all three crates. Three of the nine
had started compiling. The other six were stopped before their first `Compiling`. Two of those
had cargo at process id 353 and 362, so WSL had only just started.

Measured the same day:

- **Sleeps:** a 41 s `sleep` behind the same pipe, and a 43 s one alone, were both gone at the first
  look.
- **Signals:** a process killed by a signal inside WSL comes back through `wsl.exe` as the signal's
  number. KILL, TERM and SEGV gave 9, 15 and 11.

No case reaches the check that fails a run with a `FAILED` line in it after cargo exited 0. Every
run that printed one had exited non-zero first.

## What it does not do

- **Nothing runs it but the gate.** CI does not need it: its `the app` job runs these doctests on
  Linux already, and that is what this reproduces.
- **WSL's Ubuntu is not `ubuntu-latest`.** The glibc, the kernel and the installed packages differ.
  That does not matter for the error this exists for, which comes from the compiler and not from
  anything the doctest runs. A doctest that passed here and failed on a runner for another reason
  would still be CI's to find.
- **It runs the doctests and nothing else.** `app/`'s tests on Linux are CI's alone.
- **It does not check WSL's clock.** Cargo compares modification times. If WSL's clock ran ahead of
  Windows', an edit made just after a build would look older than the build, and the old build
  would be used. On 2026-11-17 WSL read 9 ms behind Windows, inside a 138 ms round trip.

`unearned-pass-hunter` read the first version and found four faults:

- off Windows it exited 0 on a sentence it never checked;
- one target directory served every checkout;
- it asked for the toolchain before entering `app/`, and asked for `rustc` alone;
- the numbers in this file came from earlier versions of the script.

All four are fixed above.
