# The one figure CI cannot take

`docs/editor.png` is a photograph of a window. The other two figures under `docs/` are an example's
output, and CI runs those examples on every commit — so a figure that stopped being true would take
a failing example with it. This one needs a display and a GPU, which CI has neither of.

```sh
cd app
cargo build --release --bin pantometry
powershell -File ../tools/screenshot/take.ps1
```

It dumps the frame first and reads the size it was laid out at. Then it opens the editor on
`29-a-designed-bracket-becomes-cells.json` with `--run` and waits until the window's title says the
run finished. It gives the shaded pass three seconds to draw and checks the window is the dump's
size. Only then does it capture the window, close it, and write two files: `docs/editor.png` and
`docs/editor.txt`.

## Why it captures the way it does

Each of these is a thing that went wrong once.

**`PrintWindow` with `PW_RENDERFULLCONTENT`**, not a screen grab: it asks the window to draw itself
into a bitmap, so what is captured is the application rather than whatever happens to be in front of
it on the desktop.

**The window is found by walking the process's top-level windows for the largest visible one.**
`Process.MainWindowHandle` returned a 6x6 window on two runs out of three — winit owns a small helper
window beside the real one and either can win the race, and a bitmap of that is not a screenshot.

**`[Console]::OutputEncoding` is UTF-8 around the dump.** PowerShell decodes a child process's stdout
with the console's code page, which turned the dump's em dashes into three bytes of mojibake and made
the stored file differ from a fresh one on every run. A guard that can never be satisfied is worse
than no guard.

**The text is written through `System.IO.File`**, because `Set-Content -Encoding utf8` in PowerShell
5.1 writes a byte-order mark and the file then differs from a fresh dump in its first three bytes.

**The script declares itself DPI-aware before it measures a window.** PowerShell is not, so on a
display scaled to 150% `GetWindowRect` returned the window's size divided by 1.5 while the editor,
which is DPI-aware, drew at full size: the bitmap was the window's top-left two thirds, enlarged,
and the script printed `wrote` as for any other picture. Aware, a scaled display gives a larger
picture of the same frame — both figures are 2272x1481 for that reason, where `docs/editor.png`
was 1516x989 when it was taken at 100%. It is aware *per monitor*, as the editor is. System-aware,
a window on a monitor scaled differently from the primary one would have its rectangle reported at
the primary's scale, and the size check divides that by the window's own.

**A failure stops it, and `wrote` is measured.** PowerShell 5.1 ends only the *statement* when a
.NET method throws, inside `try` as well. So a picture that failed to save went on to write a fresh
caption, and the script printed `wrote` and exited 0. An old picture under a new caption is the one
state the guard test cannot see. The script now runs under `$ErrorActionPreference = "Stop"`. It
reads the dump's exit code, which `Stop` does not reach, and refuses an empty dump rather than
writing it over the caption. It prints `wrote` only after the file is there and is not empty. Each
was sabotaged: a picture into a missing directory, a caption path that is a directory, a dump that
exits 3, and one that prints nothing. Each now exits 1 with the caption left as it was. The old
script exited 0 on the two it could reach.

**A scene that no longer runs is an exit code too.** `--ui-dump --ran` used to print
`the run refused: …` as though it were the frame, and exit 0. That walked past both checks above,
and the guard then agreed with the refusal, because a fresh dump of the same broken scene prints
the same line. It exits 1 now, with the refusal on stderr. Retaking a figure of a scene a format
change has broken stops at `--ui-dump exited 1` and leaves the caption alone. `-Out` must end in
`.png`: one ending `.txt` is its own caption.

**`-Out` may be absolute.** Both paths are resolved against the shell's location once, at the top.
The caption was `Join-Path (Get-Location) $caption`, and with an absolute `-Out` that is not a
path, so taking a picture anywhere outside `docs/` failed at the caption.

**The picture is taken of the state and the size the caption is.** Nothing tied the two. The
script slept fifteen seconds and photographed whatever was drawn, under a caption of a finished
run. Run with no wait, the old script saved a run in progress, at frame 6 of 13 with "streaming
— a prefix of the run" on screen, beside that caption, and exited 0. The size agreed only
because `1500x950` was written twice, once as the window's size and once as the dump's default.
Now:
- It is written once, as `edit::WINDOW`.
- The editor's title gains `, ran` when a run *finishes*, and not when one is stopped or fails.
  The script waits for that, for up to `-RunSeconds`.
- Before writing anything, the script measures the window's client area in points and requires it
  to equal the dump's `size=`.
- It reads the title again just before and just after the capture. The editor watches its file,
  and a write to it during the wait drops the run.
- It requires the editor to be the foreground window. A covered window is not painted, so its
  title can say `ran` over a frame from mid-run, and the answer to `SetForegroundWindow` used to
  be thrown away.

Each piece was sabotaged, and each stopped the script with nothing written:
- an expected size off by one point;
- a title that never says `ran`;
- a write to the scene during the wait.

**What the size check cannot see is egui's zoom.** The client area in points is winit's size, and
egui lays out at that over its zoom factor. A fresh window starts at a zoom of 1, because nothing
persists it, so the script, which always opens a fresh one, is safe. A window someone had pressed
Ctrl+= in would pass at the wrong layout.

**The script is ASCII.** PowerShell 5.1 reads a `.ps1` without a byte-order mark as the system code
page, so an em dash in a comment is a parse error at the line that contains it.

## The second file is the point

`docs/editor.txt` is the same frame through `--ui-dump` — the same egui layout, built with no window
and no GPU. `the_screenshots_show_the_editor_as_it_is` regenerates it and compares, so a change to
the interface fails a test instead of ageing a picture in silence. It sees the frame and not the
pixels: the same strings in the same places, the same viewport rect, the same counts.

**Agreeing is not enough, so the frame must also read `overlap=0` and `cut=0`.** A caption taken
of a broken frame agrees with a fresh dump of it. `docs/editor-protein.png` showed its colour
bar's two ends run together as `3.9674e-119.5054e-11` for exactly that reason. Its caption had
recorded the collision faithfully, as `overlap=0` at the time, because the overlap check shrank
rects sideways too.

When it fails, retake the picture with the command above. Editing `editor.txt` alone would restore
the green and leave the PNG exactly as stale as it was, which is why one script writes both.
