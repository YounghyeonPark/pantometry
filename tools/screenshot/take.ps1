# Take `docs/editor.png`: open the editor on one scene, capture that window, close it.
#
# Run from `app/`, after `cargo build --release --bin pantometry`:
#
#     powershell -File ../tools/screenshot/take.ps1
#
# **This needs a display, which is why CI does not run it.** Everything else in `docs/` is one
# command from an example that CI runs on every commit; this one is a photograph of a window. What
# stands in for it is `docs/editor.txt` -- the same frame as text, from `--ui-dump`, regenerated and
# compared by `the_screenshots_show_the_editor_as_it_is`. When that test fails the picture is
# stale and this script is how it is retaken.
param(
    [string]$Scene = "pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json",
    [string]$Out = "..\docs\editor.png",
    # How long the run may take. The script waits for the window to say the run finished rather
    # than for a fixed time, and gives up -- writing nothing -- after this.
    [int]$RunSeconds = 120,
    # After the run has finished, for the shaded pass to draw what it produced.
    [int]$SettleSeconds = 3
)

# **A failure stops the script.** PowerShell 5.1 treats an exception from a .NET method as ending
# that statement and no more, inside `try` as well, so a picture that failed to save went on to a
# fresh caption, printed `wrote` and exited 0 -- and an old picture under a new caption is the one
# state `the_screenshots_show_the_editor_as_it_is` exists to catch and cannot.
$ErrorActionPreference = "Stop"

# Both paths resolved against this shell's location once, here. `Bitmap.Save` resolves a relative
# path against the process's directory rather than PowerShell's, and `Join-Path` of the location
# and an absolute `-Out` is a path that does not exist.
$Out = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($Out)
$caption = [System.IO.Path]::ChangeExtension($Out, ".txt")
# A `-Out` ending `.txt` is its own caption: the picture was written, reported, then overwritten.
if ([System.IO.Path]::GetExtension($Out) -ne ".png") { throw "-Out must end in .png; got $Out" }

# `wrote` only for a file that is there and has something in it, measured after the write.
function Wrote([string]$path) {
    $f = Get-Item -LiteralPath $path
    if ($f.Length -eq 0) { throw "$path is empty" }
    "wrote $path ($($f.Length) bytes)"
}

Add-Type -AssemblyName System.Drawing

Add-Type @"
using System;
using System.Runtime.InteropServices;
using System.Text;
public class Win {
    public delegate bool EnumProc(IntPtr h, IntPtr p);
    [DllImport("user32.dll")] public static extern bool EnumWindows(EnumProc cb, IntPtr p);
    [DllImport("user32.dll")] public static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] public static extern bool PrintWindow(IntPtr h, IntPtr dc, uint f);
    [DllImport("user32.dll")] public static extern bool GetWindowRect(IntPtr h, out RECT r);
    [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
    [DllImport("user32.dll")] public static extern bool SetProcessDpiAwarenessContext(IntPtr c);
    [DllImport("user32.dll")] public static extern IntPtr GetForegroundWindow();
    [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr h, out RECT r);
    [DllImport("user32.dll")] public static extern uint GetDpiForWindow(IntPtr h);
    [DllImport("user32.dll", CharSet = CharSet.Unicode)] public static extern int GetWindowText(IntPtr h, StringBuilder s, int n);

    public static string Title(IntPtr h) {
        StringBuilder s = new StringBuilder(512);
        GetWindowText(h, s, s.Capacity);
        return s.ToString();
    }
    [StructLayout(LayoutKind.Sequential)]
    public struct RECT { public int Left, Top, Right, Bottom; }

    // The largest visible top-level window belonging to `pid`, or zero.
    //
    // Not `Process.MainWindowHandle`: that returned a 6x6 window on two runs out of three, because
    // winit owns a small helper window beside the real one and either can win the race. A bitmap
    // of that is not a screenshot.
    public static IntPtr Biggest(uint pid) {
        IntPtr best = IntPtr.Zero;
        int area = 0;
        EnumWindows(delegate(IntPtr h, IntPtr _) {
            uint owner;
            GetWindowThreadProcessId(h, out owner);
            if (owner != pid || !IsWindowVisible(h)) return true;
            RECT r;
            if (!GetWindowRect(h, out r)) return true;
            int a = (r.Right - r.Left) * (r.Bottom - r.Top);
            if (a > area) { area = a; best = h; }
            return true;
        }, IntPtr.Zero);
        return best;
    }
}
"@

# **DPI-aware before any window is measured.** PowerShell is not, so on a display scaled to 150%
# Windows hands it the window's rectangle divided by 1.5 while the editor, which is DPI-aware,
# draws at full size: the bitmap was the window's top-left two thirds, enlarged, and it was saved
# as a screenshot without complaint. Aware, the rectangle is the window's real one in pixels, and
# a scaled display gives a larger picture of the same frame rather than a crop of it.
# **Per monitor, as the editor is** (-4 = DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2). System-aware,
# as this was, a window on a monitor scaled differently from the primary one has its rectangle
# scaled to the primary's DPI while `GetDpiForWindow` reports its own, and the size check below
# divides one by the other.
if (-not [Win]::SetProcessDpiAwarenessContext([IntPtr](-4))) { throw "could not make this process DPI-aware per monitor -- the picture would be a crop" }

$exe = Join-Path (Get-Location) "target\release\pantometry.exe"
if (-not (Test-Path $exe)) { throw "no binary at $exe -- build it first" }
if (-not (Test-Path $Scene)) { throw "no scene at $Scene -- run this from app/" }

# **The caption first, and nothing written until the picture is known to be of it.** The caption
# was dumped after the picture and the two were tied by nothing: the picture was whatever the
# window drew fifteen seconds in, at whatever size the window had, and the caption a finished run
# at `--ui-dump`'s default. They agreed while two copies of 1500x950 agreed and the run beat the
# clock. Now the dump says what it is -- its `size=` line -- and the picture is taken only of a
# window that is that size and says its run finished.
#
# Written through .NET rather than `Set-Content -Encoding utf8`, which in PowerShell 5.1 writes a
# byte-order mark: the file then differs from a fresh dump in its first three bytes and the test
# fails on a screenshot that is perfectly current. Newlines are the dump's own, so a checkout with
# `core.autocrlf` does not change what is compared either -- the test strips carriage returns for
# the same reason.
# PowerShell decodes a child's stdout with the console's code page, which turns the dump's em
# dashes into three bytes of mojibake and makes the stored file differ from a fresh one forever.
# Told that the child speaks UTF-8, it reads them back as they were written.
$was = [Console]::OutputEncoding
[Console]::OutputEncoding = New-Object System.Text.UTF8Encoding $false
# **The caption follows -Out.** This wrote `..\docs\editor.txt` whatever -Out said, so taking a
# picture of any other scene silently replaced the committed caption for this one -- the picture
# and its caption then described different scenes, and only the guard test would have said so.
# `$caption` is set at the top, from the resolved `-Out`.
$dump = & $exe --ui-dump $Scene --ran | Out-String
$code = $LASTEXITCODE
[Console]::OutputEncoding = $was
# **The dump's exit code, read.** `Stop` does not reach a native command's exit code in 5.1, and a
# dump that failed is an empty string: written, it replaced the caption with nothing and said so.
if ($code -ne 0) { throw "--ui-dump exited $code -- $caption left as it was" }
if ($dump.Trim().Length -eq 0) { throw "--ui-dump printed nothing -- $caption left as it was" }
if ($dump -notmatch '(?m)^size=(\d+)x(\d+)\r?$') { throw "the dump has no size= line -- $caption left as it was" }
$wantW = [int]$Matches[1]
$wantH = [int]$Matches[2]

# What the window calls itself once the run it was opened with has finished: the editor adds
# `, ran` to its title then, and only for a run that finished -- not one stopped or refused.
$ranTitle = "pantometry editor " + [char]0x2014 + " " + [System.IO.Path]::GetFileName($Scene) + ", ran"

$p = Start-Process -FilePath $exe -ArgumentList @("edit", $Scene, "--run") -PassThru
try {
    $hwnd = [IntPtr]::Zero
    $r = New-Object Win+RECT
    $deadline = (Get-Date).AddSeconds(40)
    do {
        Start-Sleep -Milliseconds 300
        $hwnd = [Win]::Biggest([uint32]$p.Id)
        if ($hwnd -ne [IntPtr]::Zero) { [void][Win]::GetWindowRect($hwnd, [ref]$r) }
    } while ((($r.Right - $r.Left) -lt 400) -and (Get-Date) -lt $deadline)
    if (($r.Right - $r.Left) -lt 400) { throw "no window over 400 points wide in 40 s" }
    "window $hwnd is $($r.Right - $r.Left)x$($r.Bottom - $r.Top)"

    # **Until the run has finished, not for a fixed time.** Fifteen seconds was a guess that held
    # for scene 29 on this machine; past it, the picture was a run in progress under the caption of
    # a finished one.
    $deadline = (Get-Date).AddSeconds($RunSeconds)
    $t0 = Get-Date
    while ([Win]::Title($hwnd) -ne $ranTitle) {
        if ($p.HasExited) { throw "the editor exited before its run finished" }
        if ((Get-Date) -gt $deadline) { throw "the run had not finished after $RunSeconds s; the title is '$([Win]::Title($hwnd))' -- nothing written" }
        Start-Sleep -Milliseconds 200
    }
    "the run finished after $([int]((Get-Date) - $t0).TotalSeconds) s"

    Start-Sleep -Seconds $SettleSeconds
    # **In front, measured.** A covered window is not painted -- eframe calls `logic` instead of
    # `ui` for it -- so its title can say `ran` over a last painted frame that is mid-run. Brought
    # forward it paints, but Windows may refuse to bring it forward, and the answer was discarded.
    if (-not [Win]::SetForegroundWindow($hwnd)) { throw "Windows would not bring the editor forward -- nothing written" }
    Start-Sleep -Seconds 2
    if ([Win]::GetForegroundWindow() -ne $hwnd) { throw "the editor is not the foreground window -- nothing written" }

    # **The window is the size the caption is.** The client area in pixels, over the window's own
    # scale, is the size egui laid the frame out at. Anything else -- a display too small to hold
    # the window, or the window size and the dump's default parted -- is a picture of another
    # layout under this caption.
    $c = New-Object Win+RECT
    if (-not [Win]::GetClientRect($hwnd, [ref]$c)) { throw "GetClientRect refused" }
    $dpi = [Win]::GetDpiForWindow($hwnd)
    if ($dpi -eq 0) { throw "GetDpiForWindow refused" }
    $gotW = ($c.Right - $c.Left) * 96.0 / $dpi
    $gotH = ($c.Bottom - $c.Top) * 96.0 / $dpi
    if ([math]::Abs($gotW - $wantW) -ge 1 -or [math]::Abs($gotH - $wantH) -ge 1) {
        throw "the window is ${gotW}x${gotH} points and the caption ${wantW}x${wantH} -- nothing written"
    }
    "the window is ${gotW}x${gotH} points at $dpi dpi, as the caption is"

    [void][Win]::GetWindowRect($hwnd, [ref]$r)
    $w = $r.Right - $r.Left
    $h = $r.Bottom - $r.Top

    # **Still `ran`, on both sides of the capture.** The title was read once and the picture taken
    # five seconds later: a write to the scene in between -- the editor watches it -- rechecks and
    # drops the run, and the picture would have been of no run under a finished one's caption.
    if ([Win]::Title($hwnd) -ne $ranTitle) { throw "the title changed to '$([Win]::Title($hwnd))' before the capture -- nothing written" }
    $bmp = New-Object System.Drawing.Bitmap($w, $h)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $dc = $g.GetHdc()
    # 2 = PW_RENDERFULLCONTENT: ask the window to draw itself, so what is captured is the
    # application and not whatever happens to be in front of it on the desktop.
    $ok = [Win]::PrintWindow($hwnd, $dc, 2)
    $g.ReleaseHdc($dc)
    $g.Dispose()
    if (-not $ok) { throw "PrintWindow refused" }
    if ([Win]::Title($hwnd) -ne $ranTitle) { throw "the title changed to '$([Win]::Title($hwnd))' during the capture -- nothing written" }

    $bmp.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png)
    $bmp.Dispose()
    Wrote $Out
}
finally {
    if (-not $p.HasExited) { $p.Kill(); "closed the editor" }
}

[System.IO.File]::WriteAllText($caption, $dump, (New-Object System.Text.UTF8Encoding $false))
Wrote $caption
