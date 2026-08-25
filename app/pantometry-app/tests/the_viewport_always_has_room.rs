//! **A 3D editor must always have somewhere to draw.**
//!
//! Every `SidePanel` is laid out before the `CentralPanel`, so the panels take what they ask for and
//! the viewport gets the remainder — including when the remainder is nothing. On the default window
//! it was exactly that: three panels with minimum widths of 170, 240 and 200 points on a 533-point
//! window, and the paint callback was handed a rect of **zero width**. It drew nothing, correctly,
//! and looked precisely like a renderer that did not work — an afternoon went into the renderer
//! before the rect was printed.
//!
//! Opening the window at 1500 × 950 hid that. This is the fix, and it is checked here rather than
//! by resizing a window, because the decision is a pure function of one number.

/// The binary, run with an argument that makes it print the decision and exit.
///
/// A subcommand rather than a unit test on a private function: `panels_that_fit` lives in a `bin`
/// crate, and the alternative — a `lib.rs` for one function — would be a crate boundary invented to
/// let a test in.
fn layout_at(width: f32) -> String {
    let mut p = std::env::current_exe().expect("the test binary knows where it is");
    p.pop();
    if p.ends_with("deps") {
        p.pop();
    }
    let out =
        std::process::Command::new(p.join(format!("pantometry{}", std::env::consts::EXE_SUFFIX)))
            .args(["--layout-at", &width.to_string()])
            .output()
            .expect("the binary runs");
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// The viewport keeps its floor at every width, including absurd ones.
#[test]
fn the_view_never_loses_all_of_its_width() {
    // Every width from a phone to a wall, in steps small enough to land on each threshold.
    let mut narrowest_view = f32::MAX;
    for w in (200..3200).step_by(17) {
        let width = w as f32;
        let text = layout_at(width);
        let view: f32 = text
            .split_once("view=")
            .and_then(|(_, rest)| rest.split_whitespace().next())
            .and_then(|n| n.parse().ok())
            .unwrap_or_else(|| panic!("the binary did not report a view width: {text:?}"));
        assert!(
            view >= 280.0 || width < 280.0,
            "at {width} points the viewport got {view}: {text}"
        );
        if width >= 280.0 {
            narrowest_view = narrowest_view.min(view);
        }
    }
    println!("  the narrowest the view ever got: {narrowest_view} points");
    assert!(narrowest_view >= 280.0);
}

/// The panels come back on their own when the window widens.
///
/// The `show_*` flags are what the reader asked for and are never written by the fitting; a window
/// that grew back would otherwise leave somebody hunting for a menu item they never touched.
#[test]
fn a_wider_window_gets_its_panels_back() {
    let narrow = layout_at(600.0);
    let wide = layout_at(1600.0);
    println!("  600: {}", narrow.trim());
    println!("  1600: {}", wide.trim());

    assert!(
        narrow.contains("inspector=off"),
        "600 points cannot hold all three and the view: {narrow}"
    );
    for panel in ["outliner=on", "text=on", "inspector=on"] {
        assert!(
            wide.contains(panel),
            "1600 points holds everything, but {panel} is not set: {wide}"
        );
    }
}

/// They go in the order they are worth least to somebody looking at a picture.
#[test]
fn the_outliner_is_the_last_to_go() {
    // Narrow enough that only one panel can stay beside the view.
    let text = layout_at(470.0);
    println!("  470: {}", text.trim());
    assert!(
        text.contains("outliner=on") && text.contains("text=off") && text.contains("inspector=off"),
        "the outliner is how you choose what to look at, so it survives longest: {text}"
    );
}
