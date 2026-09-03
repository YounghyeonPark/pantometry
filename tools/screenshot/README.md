# The one figure CI cannot take

`docs/editor.png` is a photograph of a window. The other two figures under `docs/` are an example's
output, and CI runs those examples on every commit — so a figure that stopped being true would take
a failing example with it. This one needs a display and a GPU, which CI has neither of.

```sh
cd app
cargo build --release --bin pantometry
powershell -File ../tools/screenshot/take.ps1
```

It opens the editor on `29-a-designed-bracket-becomes-cells.json` with `--run`, waits for the run to
finish and the shaded pass to draw, captures the window, closes it, and writes two files:
`docs/editor.png` and `docs/editor.txt`.

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

**The script is ASCII.** PowerShell 5.1 reads a `.ps1` without a byte-order mark as the system code
page, so an em dash in a comment is a parse error at the line that contains it.

## The second file is the point

`docs/editor.txt` is the same frame through `--ui-dump` — the same egui layout, built with no window
and no GPU. `the_screenshot_shows_the_editor_as_it_is` regenerates it and compares, so a change to
the interface fails a test instead of ageing a picture in silence. It sees the frame and not the
pixels: the same strings in the same places, the same viewport rect, the same counts.

When it fails, retake the picture with the command above. Editing `editor.txt` alone would restore
the green and leave the PNG exactly as stale as it was, which is why one script writes both.
