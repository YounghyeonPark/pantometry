# Take `docs/editor.png`: open the editor on one scene, capture that window, close it.
#
# Run from `app/`, after `cargo build --release --bin pantometry`:
#
#     powershell -File ../tools/screenshot/take.ps1
#
# **This needs a display, which is why CI does not run it.** Everything else in `docs/` is one
# command from an example that CI runs on every commit; this one is a photograph of a window. What
# stands in for it is `docs/editor.txt` -- the same frame as text, from `--ui-dump`, regenerated and
# compared by `the_screenshot_shows_the_editor_as_it_is`. When that test fails the picture is
# stale and this script is how it is retaken.
param(
    [string]$Scene = "pantometry-world/scenes/29-a-designed-bracket-becomes-cells.json",
    [string]$Out = "..\docs\editor.png",
    # Long enough for the run to finish and the shaded pass to draw. Scene 29 is 4 100 cells over
    # 13 frames and settles well inside this.
    [int]$SettleSeconds = 15
)

Add-Type -AssemblyName System.Drawing

Add-Type @"
using System;
using System.Runtime.InteropServices;
public class Win {
    public delegate bool EnumProc(IntPtr h, IntPtr p);
    [DllImport("user32.dll")] public static extern bool EnumWindows(EnumProc cb, IntPtr p);
    [DllImport("user32.dll")] public static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
    [DllImport("user32.dll")] public static extern bool IsWindowVisible(IntPtr h);
    [DllImport("user32.dll")] public static extern bool PrintWindow(IntPtr h, IntPtr dc, uint f);
    [DllImport("user32.dll")] public static extern bool GetWindowRect(IntPtr h, out RECT r);
    [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
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

$exe = Join-Path (Get-Location) "target\release\pantometry.exe"
if (-not (Test-Path $exe)) { throw "no binary at $exe -- build it first" }
if (-not (Test-Path $Scene)) { throw "no scene at $Scene -- run this from app/" }

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

    Start-Sleep -Seconds $SettleSeconds
    [void][Win]::SetForegroundWindow($hwnd)
    Start-Sleep -Seconds 2
    [void][Win]::GetWindowRect($hwnd, [ref]$r)
    $w = $r.Right - $r.Left
    $h = $r.Bottom - $r.Top

    $bmp = New-Object System.Drawing.Bitmap($w, $h)
    $g = [System.Drawing.Graphics]::FromImage($bmp)
    $dc = $g.GetHdc()
    # 2 = PW_RENDERFULLCONTENT: ask the window to draw itself, so what is captured is the
    # application and not whatever happens to be in front of it on the desktop.
    $ok = [Win]::PrintWindow($hwnd, $dc, 2)
    $g.ReleaseHdc($dc)
    $g.Dispose()
    if (-not $ok) { throw "PrintWindow refused" }

    $bmp.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png)
    $bmp.Dispose()
    "wrote $Out"
}
finally {
    if (-not $p.HasExited) { $p.Kill(); "closed the editor" }
}

# And the frame as text, which is what notices when the picture above stops being true.
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
$dump = & $exe --ui-dump $Scene --ran | Out-String
[Console]::OutputEncoding = $was
[System.IO.File]::WriteAllText((Join-Path (Get-Location) "..\docs\editor.txt"), $dump, (New-Object System.Text.UTF8Encoding $false))
"wrote ..\docs\editor.txt"
