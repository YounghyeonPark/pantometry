//! **A sweep from the command line and from the editor's window, and what each says.**
//!
//! `editor_core::sweep` is held to scene 28's eigenvalues in its own tests. This holds the two doors
//! a person uses — `pantometry sweep` and the Sweep window — to reaching it: the right number swept,
//! every value run, the readings the table and the CSV print, and refusals said rather than run.
//!
//! Scene 28 is an electron in the third state of a hard-walled well, and two of its readings have a
//! closed form at every width: the energy, the discrete Hamiltonian's eigenvalue, which goes as one
//! over the width squared; and the position expectation, the middle of the well, by the state's
//! symmetry.

#![cfg(not(target_family = "wasm"))]

const SCENE: &str = "28-an-eigenstate-that-does-not-move.json";

fn scene() -> String {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes")
        .join(SCENE)
        .to_string_lossy()
        .into_owned()
}

fn pantometry(args: &[&str]) -> std::process::Output {
    let none = std::env::temp_dir().join("pantometry-no-such-recent-list.json");
    std::process::Command::new(env!("CARGO_BIN_EXE_pantometry"))
        .args(args)
        .env("PANTOMETRY_RECENT", none.to_string_lossy().into_owned())
        .output()
        .expect("the binary runs")
}

/// The energy of the third state of a well `width` metres wide on two hundred points.
fn level(width: f64) -> f64 {
    const HBAR: f64 = 1.054_571_817e-34;
    const ELECTRON_KG: f64 = 9.109_383_701_5e-31;
    let dx = width / 201.0;
    let s = (3.0 * std::f64::consts::PI / 402.0).sin();
    2.0 * HBAR * HBAR / (ELECTRON_KG * dx * dx) * s * s
}

/// **`pantometry sweep` prints the table on stdout and the progress on stderr.** Three widths, three
/// rows, each energy its own width's eigenvalue to `1e-9` — the CSV carries every figure, which is
/// the point of reading a sweep's numbers out of its runs rather than out of the run files' seven.
#[test]
fn the_command_line_tabulates_each_widths_level() {
    let out = pantometry(&["sweep", &scene(), "particle.width_nm=5,10,20"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let said = String::from_utf8_lossy(&out.stderr);
    for (i, w) in [5, 10, 20].iter().enumerate() {
        assert!(
            said.contains(&format!("{} of 3: particle.width_nm = {w} ran", i + 1)),
            "{said}"
        );
    }
    let csv = String::from_utf8_lossy(&out.stdout);
    let mut lines = csv.lines();
    let header: Vec<&str> = lines.next().expect("a header").split(',').collect();
    let e = header
        .iter()
        .position(|h| *h == "particle.<E> [J]")
        .unwrap_or_else(|| panic!("no energy column: {header:?}"));
    assert_eq!(header[0], "particle.width_nm");
    assert_eq!(header.last(), Some(&"error"));
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split(',').collect()).collect();
    assert_eq!(rows.len(), 3, "{csv}");
    for row in &rows {
        let w: f64 = row[0].parse().expect("a width");
        let got: f64 = row[e].parse().expect("an energy");
        let want = level(w * 1e-9);
        assert!(
            (got / want - 1.0).abs() < 1e-9,
            "{w} nm: {got:e} against {want:e}"
        );
        assert_eq!(row.last(), Some(&""), "a row that ran says no error");
    }
}

/// **`--list` names every number, the pointer beside it; a file named takes the table; and what
/// cannot be swept is refused with the reason, exit 1, before anything runs.**
#[test]
fn the_command_line_names_its_numbers_and_refuses_what_it_cannot_sweep() {
    let out = pantometry(&["sweep", &scene(), "--list"]);
    assert!(out.status.success());
    let listed = String::from_utf8_lossy(&out.stdout);
    assert!(
        listed
            .lines()
            .any(|l| l == "particle.width_nm\t/domains/0/width_nm\t10\tnm"),
        "{listed}"
    );
    assert!(
        listed
            .lines()
            .any(|l| l.starts_with("duration_s\t/duration_s\t2e-13\t")),
        "{listed}"
    );

    let file = std::env::temp_dir().join(format!("pantometry-sweep-{}.csv", std::process::id()));
    let path = file.to_string_lossy().into_owned();
    let out = pantometry(&["sweep", &scene(), "/domains/0/cells=100:300:3", &path]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        out.stdout.is_empty(),
        "the table went to the file, not stdout"
    );
    let written = std::fs::read_to_string(&file).expect("the file was written");
    let _ = std::fs::remove_file(&file);
    let firsts: Vec<&str> = written
        .lines()
        .skip(1)
        .map(|l| l.split(',').next().unwrap_or(""))
        .collect();
    assert_eq!(firsts, ["1e2", "2e2", "3e2"], "{written}");
    // And each of them ran: the error column, last, is empty — a row that refused would still be a
    // row, and the exit is 0 while any value runs.
    for row in written.lines().skip(1) {
        assert!(row.ends_with(','), "a written row says it refused: {row}");
    }

    for (args, says) in [
        (
            vec!["particle.depth_nm=1,2"],
            "no number called particle.depth_nm",
        ),
        (vec!["particle.cells=150.5"], "whole number"),
        (vec!["particle.width_nm=5,ten"], "not a number"),
        (vec!["particle.width_nm"], "names no values"),
        (vec!["particle.width_nm=0:2:3:log"], "above zero"),
    ] {
        let mut all = vec!["sweep", SCENE];
        let s = scene();
        all[1] = &s;
        all.extend(args.iter().copied());
        let out = pantometry(&all);
        let err = String::from_utf8_lossy(&out.stderr);
        assert_eq!(out.status.code(), Some(1), "{args:?} ran: {err}");
        assert!(err.contains(says), "{args:?} should say {says:?}: {err}");
        assert!(out.stdout.is_empty(), "{args:?} printed a table");
    }

    // **Every value run and every one refused** is a sweep that answered nothing, exit 1 — and the
    // table of why is still printed, each row its reason. The cases above are refused before any
    // run, so none of them reached this.
    let s = scene();
    let out = pantometry(&["sweep", &s, "particle.cells=0,1"]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(1), "{err}");
    assert!(err.contains("no value of particle.cells ran"), "{err}");
    let table = String::from_utf8_lossy(&out.stdout);
    assert_eq!(
        table
            .lines()
            .filter(|l| l.contains("needs at least three interior cells"))
            .count(),
        2,
        "{table}"
    );
}

fn dump(args: &[&str]) -> String {
    let mut all = vec!["--ui-dump", ""];
    let s = scene();
    all[1] = &s;
    all.extend(args.iter().copied());
    let out = pantometry(&all);
    assert!(
        out.status.success(),
        "--ui-dump {args:?}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// One drawn string: left, top, and what it says.
fn drawn(dump: &str) -> Vec<(f32, f32, String)> {
    let mut out = Vec::new();
    for line in dump.lines() {
        if line.len() < 3 || !line.is_char_boundary(2) {
            continue;
        }
        let rest = &line[2..];
        let mut it = rest.split_whitespace();
        let (Some(x), Some(y), Some(w)) = (it.next(), it.next(), it.next()) else {
            continue;
        };
        let (Ok(x), Ok(y), Ok(_)) = (x.parse::<f32>(), y.parse::<f32>(), w.parse::<f32>()) else {
            continue;
        };
        let after = rest
            .split_whitespace()
            .take(3)
            .fold(0usize, |at, f| rest[at..].find(f).unwrap() + at + f.len());
        out.push((x, y, rest[after..].trim_start().to_string()));
    }
    out
}

fn at(dump: &str, text: &str) -> String {
    let (x, y, _) = drawn(dump)
        .into_iter()
        .find(|(_, _, t)| t == text)
        .unwrap_or_else(|| panic!("nothing says {text:?}:\n{dump}"));
    format!("{},{}", x + 4.0, y + 4.0)
}

/// **Run → Sweep a value… opens the window, set up and ready to run.** A value of the scene is
/// chosen — one a domain states — a range is set round it, and the button to run it is there.
#[test]
fn the_run_menu_opens_the_sweep_window() {
    let before = dump(&[]);
    assert!(
        !before.contains("Run sweep"),
        "the window is open before it is asked for"
    );
    let menu = at(&before, "Run");
    let opened = dump(&["--click", &menu]);
    let item = at(&opened, "Sweep a value…");
    let window = dump(&["--click", &menu, "--click", &item]);
    for says in ["Sweep", "Run sweep", "or the values", "particle.cells"] {
        assert!(
            drawn(&window).iter().any(|(_, _, t)| t == says),
            "the window does not say {says:?}:\n{window}"
        );
    }
    // A whole number swept from where it is, by whole numbers: 200 to 204, five runs.
    assert!(
        window.contains("5 run(s), at 200, 201, 202, 203, 204"),
        "{window}"
    );
}

/// **The window's table is the sweep that ran, one row a value, and its chart is that reading
/// against the value.** The position expectation opens the table, being the first reading that
/// moves, and it is the middle of the well at each width: 2.5, 5 and 10 nm, at the four figures a
/// table prints. The chart's horizontal axis is the value in the unit its key names.
#[test]
fn the_windows_table_is_the_sweep_that_ran() {
    let d = dump(&["--sweep", "particle.width_nm=5,10,20"]);
    let rows = drawn(&d);
    let column = rows
        .iter()
        .find(|(_, _, t)| t == "particle.width_nm / nm")
        .map(|(x, y, _)| (*x, *y))
        .unwrap_or_else(|| panic!("no table header:\n{d}"));
    let reading = rows
        .iter()
        .find(|(x, y, t)| {
            (*y - column.1).abs() < 3.0 && *x > column.0 && t.starts_with("particle ")
        })
        .map(|(x, _, t)| (*x, t.clone()))
        .unwrap_or_else(|| panic!("no reading column:\n{d}"));
    assert_eq!(reading.1, "particle <x> / m");
    let mut table = Vec::new();
    for (x, y, t) in &rows {
        if (*x - column.0).abs() < 2.0 && *y > column.1 + 5.0 && t.parse::<f64>().is_ok() {
            let value = rows
                .iter()
                .find(|(x2, y2, _)| (*y2 - y).abs() < 3.0 && (*x2 - reading.0).abs() < 2.0)
                .map(|(_, _, v)| v.clone())
                .unwrap_or_default();
            table.push((t.clone(), value));
        }
    }
    assert_eq!(
        table,
        [
            ("5".to_string(), "2.5000e-9".to_string()),
            ("10".to_string(), "5.0000e-9".to_string()),
            ("20".to_string(), "1.0000e-8".to_string()),
        ],
        "{d}"
    );
    assert!(
        rows.iter()
            .filter(|(_, _, t)| t == "particle.width_nm / nm")
            .count()
            >= 2,
        "the chart's axis is not labelled with the value:\n{d}"
    );
    assert!(d.contains("3 run(s), at 5, 10, 20 nm"), "{d}");

    let out = pantometry(&["--ui-dump", &scene(), "--sweep", "particle.depth_nm=1,2"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("--sweep: the scene has no number"));
}

/// **The window's own Run sweep runs it**, in the background, as the editor does — not only the
/// synchronous `--sweep` door. Run → Sweep a value… → Run sweep, on the range the window set itself
/// round the well's 200 points: five runs, 200 to 204, each a row of the table, and the status says
/// they ran.
#[test]
fn the_windows_run_button_runs_the_sweep() {
    let menu = at(&dump(&[]), "Run");
    let item = at(&dump(&["--click", &menu]), "Sweep a value…");
    let window = dump(&["--click", &menu, "--click", &item]);
    let button = at(&window, "Run sweep");
    let ran = dump(&["--click", &menu, "--click", &item, "--click", &button]);
    // On the status bar, and above the table, which the next thing the editor does cannot take.
    assert_eq!(
        drawn(&ran)
            .iter()
            .filter(|(_, _, t)| t == "swept 5 value(s)")
            .count(),
        2,
        "{ran}"
    );
    let rows = drawn(&ran);
    // The table's header is at the window's left edge, where its other rows' labels are; the value's
    // menu and the chart's axis say `particle.cells` too, elsewhere.
    let left = rows
        .iter()
        .find(|(_, _, t)| t == "or the values")
        .map(|(x, _, _)| *x)
        .unwrap_or_else(|| panic!("no window:\n{ran}"));
    let header = rows
        .iter()
        .find(|(x, _, t)| t == "particle.cells" && (*x - left).abs() < 2.0)
        .map(|(_, y, _)| *y)
        .unwrap_or_else(|| panic!("no table:\n{ran}"));
    let column: Vec<&str> = rows
        .iter()
        .filter(|(x, y, t)| {
            (*x - left).abs() < 2.0 && *y > header + 5.0 && t.parse::<f64>().is_ok()
        })
        .map(|(_, _, t)| t.as_str())
        .collect();
    assert_eq!(column, ["200", "201", "202", "203", "204"], "{ran}");
}

/// **A scene that reports no readings has no table, and the command line says so.** Scene 06's
/// orbits are a picture: a sweep of its sun's mass runs at every value and has nothing to tabulate.
/// It wrote a CSV of one column and exited 0.
#[test]
fn a_sweep_with_nothing_to_tabulate_says_so() {
    let orbits = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/06-orbits.json")
        .to_string_lossy()
        .into_owned();
    let out = pantometry(&["sweep", &orbits, "system.central_kg=1,2"]);
    let err = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(1), "{err}");
    assert!(
        err.contains("2 value(s) of system.central_kg ran") && err.contains("reports no readings"),
        "{err}"
    );
    assert!(out.stdout.is_empty(), "a table of nothing was printed");
}

/// **A reading only some values report goes in a column of its own**, and the readings after it
/// stay in theirs. Scene 25's block reports `generated` only while its source is on, so at 45 W its
/// frame has one more reading than at 0 W; every structure reading after it was blank at 45 W,
/// matched by place. Both rows carry the structure's strain energy now, and the 0 W row has no
/// `generated`.
#[test]
fn a_reading_only_some_values_report_has_its_own_column() {
    let module = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/25-what-140-kelvin-does-to-the-solder.json")
        .to_string_lossy()
        .into_owned();
    let out = pantometry(&["sweep", &module, "module.dissipation[0].watts=0,45"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let csv = String::from_utf8_lossy(&out.stdout);
    let mut lines = csv.lines();
    let header: Vec<&str> = lines.next().expect("a header").split(',').collect();
    let col = |name: &str| {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("no {name}: {header:?}"))
    };
    let (strain, generated) = (col("stress.strain energy [J]"), col("module.generated [J]"));
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split(',').collect()).collect();
    assert_eq!(rows.len(), 2, "{csv}");
    for row in &rows {
        assert!(
            row[strain].parse::<f64>().is_ok(),
            "no strain energy at {}: {csv}",
            row[0]
        );
    }
    assert!(
        rows[0][generated].is_empty(),
        "a source of 0 W generated: {csv}"
    );
    assert!(
        rows[1][generated].parse::<f64>().is_ok_and(|g| g > 0.0),
        "{csv}"
    );
}
