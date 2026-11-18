//! **A sweep runs the scene once per value and reads each run's answer.**
//!
//! [`editor_core::sweep`] is plumbing — a value written into the text, a run, the last frame's
//! readings — and plumbing is held to something that does not come from it. Scene 28 is an electron
//! in the third state of a hard-walled well, and its energy is the discrete Hamiltonian's exact
//! eigenvalue, `(2ħ²/m·dx²)·sin²(nπ/2(N+1))` with `dx = width/(N+1)`. So a sweep of the width has an
//! answer at every value that is written out here from the constants, and a sweep that edited the
//! wrong number, read the wrong reading or kept the wrong frame would miss it.
//!
//! # Not under `wasm32`
//!
//! The scenes are read off a disk, and a `wasm32` target has none.

#![cfg(not(target_family = "wasm"))]

use editor_core::sweep::{self, Point};
use std::sync::atomic::AtomicBool;

fn scene(name: &str) -> String {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes");
    let path = std::fs::read_dir(&dir)
        .expect("the scenes directory reads")
        .map(|e| e.expect("an entry").path())
        .find(|p| {
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with(name))
        })
        .unwrap_or_else(|| panic!("no scene {name}"));
    std::fs::read_to_string(path).expect("the scene reads")
}

fn files() -> editor_core::Beside {
    editor_core::Beside(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes"),
    )
}

/// The energy of the third state of a well `width` metres wide on `cells` interior points: the
/// discrete eigenvalue, from the constants.
fn level(width: f64, cells: f64) -> f64 {
    const HBAR: f64 = 1.054_571_817e-34;
    const ELECTRON_KG: f64 = 9.109_383_701_5e-31;
    let dx = width / (cells + 1.0);
    let theta = 3.0 * std::f64::consts::PI / (2.0 * (cells + 1.0));
    let s = theta.sin();
    2.0 * HBAR * HBAR / (ELECTRON_KG * dx * dx) * s * s
}

fn energy(p: &Point) -> f64 {
    p.readings
        .iter()
        .find(|r| r.label == "<E>")
        .unwrap_or_else(|| panic!("no <E> at {}: {:?}", p.value, p.error))
        .value
}

/// **A sweep of the well's width reads the level each width has.** 5, 10 and 20 nm on scene 28's
/// two hundred points: each run's energy is its own width's eigenvalue to the `1e-9` the scene's own
/// test holds the solver to, which is the level rotating at the operator's rate. So the levels fall
/// as the width squared, four times between each pair, and the 10 nm point is the scene as shipped.
#[test]
fn a_wells_level_falls_as_its_width_squared() {
    let text = scene("28-");
    let width = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.width_nm")
        .expect("the width is a value to sweep");
    assert_eq!((width.now, width.integral, width.unit), (10.0, false, "nm"));
    let stop = AtomicBool::new(false);
    let mut seen = Vec::new();
    let points = sweep::run_all(
        &text,
        &width.pointer,
        &[5.0, 10.0, 20.0],
        &files(),
        &stop,
        |i, p| {
            seen.push((i, p.value));
        },
    );
    assert_eq!(
        seen,
        [(0, 5.0), (1, 10.0), (2, 20.0)],
        "each point as it ran"
    );
    for p in &points {
        assert!(
            p.error.is_none(),
            "{} nm did not run: {:?}",
            p.value,
            p.error
        );
        let want = level(p.value * 1e-9, 200.0);
        let got = energy(p);
        assert!(
            (got / want - 1.0).abs() < 1e-9,
            "{} nm reads {got:e} J and the eigenvalue is {want:e}",
            p.value
        );
    }
    let ratio = energy(&points[0]) / energy(&points[1]);
    assert!(
        (ratio - 4.0).abs() < 1e-8,
        "halving the width multiplied the level by {ratio}"
    );
}

/// **A point at the scene's own value is the scene's own run.** The variant at the value the file
/// already states is the file, byte for byte; and its readings are the run button's, reading for
/// reading, at the seven figures the run file writes — the point keeps every figure, the file six
/// after the first.
#[test]
fn the_scenes_own_value_is_the_scenes_own_run() {
    let text = scene("28-");
    let width = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.width_nm")
        .expect("a width");
    assert_eq!(
        sweep::variant(&text, &width.pointer, width.now).as_deref(),
        Ok(text.as_str()),
        "writing the value the file has changed the file"
    );
    let point = sweep::run_point(&text, &width.pointer, width.now, &files());
    let run = viewer_core::Run::from_json(&editor_core::run(&text, &files()).expect("it runs"))
        .expect("it reads");
    let last = run.frames.last().expect("frames");
    assert_eq!(point.readings.len(), last.readings.len());
    for (a, b) in point.readings.iter().zip(&last.readings) {
        assert_eq!(
            (
                a.domain.as_str(),
                a.label.as_str(),
                format!("{:.6e}", a.value)
            ),
            (
                b.domain.as_str(),
                b.label.as_str(),
                format!("{:.6e}", b.value)
            ),
        );
    }
}

/// **A rate is labelled as a rate, every one of them.** The unit a key names was read off its last
/// word, so the settings form labelled a film coefficient `K` and a stiffness `m`.
///
/// Every key a shipped scene states with `_per_` in it is counted here off the scene's own JSON, and
/// the sweep must list the same number of them — so a rate the listing dropped cannot leave this
/// unchecked — and each must carry the unit written below for its key. The table is the
/// specification: the format's rates, by what they are rates of.
#[test]
fn a_rate_is_labelled_as_a_rate() {
    const RATES: [(&str, &str); 16] = [
        ("convection_w_per_m2_k", "W/m²·K"),
        ("contact_w_per_m2_k", "W/m²·K"),
        ("w_per_m2_k", "W/m²·K"),
        ("w_per_k", "W/K"),
        ("conductivity_w_per_m_k", "W/m·K"),
        ("w_per_m_k", "W/m·K"),
        ("gravity_m_per_s2", "m/s²"),
        ("drive_m_per_s2", "m/s²"),
        ("clearance_l_per_h", "L/h"),
        ("flow_l_per_h", "L/h"),
        ("absorption_per_h", "1/h"),
        ("amplitude_v_per_m", "V/m"),
        ("stiffness_n_per_m", "N/m"),
        ("k0_per_nm", "1/nm"),
        ("lower_m_per_s", "m/s"),
        ("upper_m_per_s", "m/s"),
    ];
    // Every number in the scene whose key says `_per_`, walked here rather than listed by the code
    // under test.
    fn count(v: &serde_json::Value, key: &str) -> usize {
        match v {
            serde_json::Value::Object(m) => m.iter().map(|(k, v)| count(v, k)).sum(),
            serde_json::Value::Array(a) => a.iter().map(|v| count(v, key)).sum(),
            serde_json::Value::Number(_) if key.contains("_per_") || key.starts_with("w_per") => 1,
            _ => 0,
        }
    }
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../pantometry-world/scenes");
    let mut rates = 0;
    for entry in std::fs::read_dir(&dir).expect("reads") {
        let path = entry.expect("an entry").path();
        if path.extension().is_none_or(|x| x != "json") {
            continue;
        }
        let text = std::fs::read_to_string(&path).expect("reads");
        let root: serde_json::Value = serde_json::from_str(&text).expect("parses");
        let mut listed = 0;
        for p in sweep::parameters(&text) {
            let key = p.label.rsplit('.').next().unwrap_or(&p.label);
            let key = key.split('[').next().unwrap_or(key);
            if !(key.contains("_per_") || key.starts_with("w_per")) {
                continue;
            }
            listed += 1;
            let want = RATES
                .iter()
                .find(|(k, _)| *k == key)
                .unwrap_or_else(|| {
                    panic!(
                        "{}: {key} is a rate this table does not know",
                        path.display()
                    )
                })
                .1;
            assert_eq!(p.unit, want, "{}: {}", path.display(), p.label);
        }
        assert_eq!(listed, count(&root, ""), "{}: rates listed", path.display());
        rates += listed;
    }
    assert!(rates > 0, "no rate in any scene: the walk checked nothing");
    assert_eq!(
        sweep::parameters(&scene("28-"))
            .into_iter()
            .find(|p| p.label == "particle.width_nm")
            .map(|p| p.unit),
        Some("nm")
    );
}

/// **A whole number stays whole.** The well's `cells` is written `200`: a sweep of it asks for whole
/// numbers, ten asked between 2 and 5 are the four there are, and a fraction written into it is
/// refused rather than rounded into a different question. The level at 100 and 300 points is that
/// grid's eigenvalue — `θ` moves with `N` here, which a sweep that kept the old grid would miss.
#[test]
fn a_whole_number_stays_whole() {
    let text = scene("28-");
    let cells = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.cells")
        .expect("the grid is a value to sweep");
    assert!(cells.integral, "{cells:?}");
    assert_eq!(
        sweep::spaced(2.0, 5.0, 10, false, true),
        Ok(vec![2.0, 3.0, 4.0, 5.0])
    );
    let refused = sweep::variant(&text, &cells.pointer, 150.5).expect_err("a fraction written");
    assert!(refused.contains("whole number"), "{refused}");
    let points = sweep::run_all(
        &text,
        &cells.pointer,
        &[100.0, 300.0],
        &files(),
        &AtomicBool::new(false),
        |_, _| {},
    );
    assert_eq!(
        points.iter().map(|p| p.value).collect::<Vec<_>>(),
        [100.0, 300.0]
    );
    for p in &points {
        let want = level(10e-9, p.value);
        assert!(
            (energy(p) / want - 1.0).abs() < 1e-9,
            "{} cells: {:e} against {want:e}",
            p.value,
            energy(p)
        );
    }
}

/// **Spacing is what it says.** Linear values step by one amount and reach both ends exactly; a
/// logarithmic sweep steps by one ratio. Nothing to run, a range that is not two numbers and a
/// logarithm of zero are refused.
#[test]
fn spacing_is_what_it_says() {
    assert_eq!(
        sweep::spaced(5.0, 20.0, 4, false, false),
        Ok(vec![5.0, 10.0, 15.0, 20.0])
    );
    let log = sweep::spaced(1.0, 1000.0, 4, true, false).expect("spaced");
    assert_eq!((log[0], log[3]), (1.0, 1000.0), "the ends exactly");
    for pair in log.windows(2) {
        assert!((pair[1] / pair[0] - 10.0).abs() < 1e-12, "{log:?}");
    }
    assert_eq!(sweep::spaced(7.0, 9.0, 1, false, false), Ok(vec![7.0]));
    for (from, to, count, log) in [
        (1.0, 2.0, 0, false),
        (f64::NAN, 2.0, 3, false),
        (0.0, 2.0, 3, true),
        (-1.0, 2.0, 3, true),
    ] {
        assert!(
            sweep::spaced(from, to, count, log, false).is_err(),
            "{from} to {to}, {count}, log {log}"
        );
    }
    assert_eq!(sweep::listed("5, 10,20"), Ok(vec![5.0, 10.0, 20.0]));
    for bad in ["", "5,ten", "5,inf"] {
        assert!(sweep::listed(bad).is_err(), "{bad:?}");
    }
}

/// **A value the scene refuses is an answer, kept beside the value**, and the sweep goes on: a
/// well of negative width does not build, and the 10 nm after it still runs.
#[test]
fn a_refusal_is_kept_beside_its_value() {
    let text = scene("28-");
    let width = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.width_nm")
        .expect("a width");
    let points = sweep::run_all(
        &text,
        &width.pointer,
        &[-1.0, 10.0],
        &files(),
        &AtomicBool::new(false),
        |_, _| {},
    );
    assert_eq!(points.len(), 2);
    assert!(
        points[0].error.is_some() && points[0].readings.is_empty(),
        "{:?}",
        points[0]
    );
    assert!(points[1].error.is_none(), "{:?}", points[1].error);

    // The table says it in its own column, and the row that ran has its numbers.
    let csv = sweep::csv(&width.label, &points);
    let lines: Vec<&str> = csv.lines().collect();
    assert!(
        lines[0].starts_with("particle.width_nm,") && lines[0].ends_with(",error"),
        "{}",
        lines[0]
    );
    assert!(lines[0].contains("particle.<E> [J]"), "{}", lines[0]);
    assert!(
        lines[1].starts_with("-1e0,") && lines[1].len() > "-1e0,,,,".len(),
        "{}",
        lines[1]
    );
    assert!(
        lines[2].starts_with("1e1,") && lines[2].ends_with(','),
        "{}",
        lines[2]
    );
}

/// **A stopped sweep is the prefix it ran**, and says nothing about the values it did not.
#[test]
fn a_stopped_sweep_is_what_it_ran() {
    let text = scene("28-");
    let width = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.width_nm")
        .expect("a width");
    let stop = AtomicBool::new(false);
    let points = sweep::run_all(
        &text,
        &width.pointer,
        &[8.0, 9.0, 10.0],
        &files(),
        &stop,
        |i, _| {
            if i == 0 {
                stop.store(true, std::sync::atomic::Ordering::Relaxed);
            }
        },
    );
    assert_eq!(points.iter().map(|p| p.value).collect::<Vec<_>>(), [8.0]);
}

/// **A column is a reading by its name, wherever a frame lists it.** A frame need not report the
/// same readings at every value — a block's `generated` appears only while its source is on — so a
/// reading that arrives at the second value goes in a column of its own, and every reading after it
/// in that frame still lands in its own column rather than one place over. A name a frame reports
/// twice is two columns, by which of them it is.
#[test]
fn a_column_is_a_reading_by_its_name() {
    let r = |label: &str, value: f64| sweep::Reading {
        domain: "d".into(),
        label: label.into(),
        unit: "K".into(),
        value,
    };
    let points = [
        Point {
            value: 0.0,
            readings: vec![r("peak", 300.0), r("strain", 1.0)],
            error: None,
        },
        Point {
            value: 45.0,
            readings: vec![r("peak", 310.0), r("generated", 9.0), r("strain", 2.0)],
            error: None,
        },
        Point {
            value: 90.0,
            readings: vec![r("peak", 320.0), r("peak", 321.0)],
            error: None,
        },
    ];
    let cols = sweep::columns(&points);
    let names: Vec<(String, usize, usize)> = cols
        .iter()
        .map(|c| (c.label.clone(), c.nth, c.reported))
        .collect();
    assert_eq!(
        names,
        [
            ("peak".to_string(), 0, 3),
            ("strain".to_string(), 0, 2),
            ("generated".to_string(), 0, 1),
            ("peak".to_string(), 1, 1),
        ]
    );
    assert_eq!(
        sweep::at(&points[1], &cols, 1),
        Some(2.0),
        "the strain after generated"
    );
    assert_eq!(
        sweep::at(&points[0], &cols, 2),
        None,
        "no generated with the source off"
    );
    assert_eq!(
        sweep::at(&points[2], &cols, 3),
        Some(321.0),
        "the second peak"
    );
}

/// **What a table opens on is the first reading that moves, measured against its own size.** A
/// column that starts at zero still moves when it rises from there; one point has nothing to move
/// against; and a column of one number does not move at any size.
#[test]
fn the_reading_a_table_opens_on_is_one_that_moves() {
    let point = |value: f64, flat: f64, rising: f64| Point {
        value,
        readings: vec![
            sweep::Reading {
                domain: "d".into(),
                label: "flat".into(),
                unit: String::new(),
                value: flat,
            },
            sweep::Reading {
                domain: "d".into(),
                label: "rising".into(),
                unit: String::new(),
                value: rising,
            },
        ],
        error: None,
    };
    let points = [point(1.0, 7.0, 0.0), point(2.0, 7.0, 3.0)];
    let cols = sweep::columns(&points);
    assert_eq!(sweep::moving(&points, &cols), Some(1));
    assert_eq!(
        sweep::moving(&points[..1], &cols),
        None,
        "one point cannot move"
    );
    let flat = [point(1.0, 7.0, 0.0), point(2.0, 7.0, 0.0)];
    assert_eq!(sweep::moving(&flat, &sweep::columns(&flat)), None);
}

/// **The chart runs in the value's order and breaks where a value has no reading.** Values written
/// `20, 5, 10` are drawn 5, 10, 20, and a value the scene refused is a gap in the line rather than a
/// straight segment across it.
#[test]
fn the_chart_is_in_order_and_broken_at_a_refusal() {
    let at = |value: f64, e: Option<f64>| Point {
        value,
        readings: e
            .map(|v| {
                vec![sweep::Reading {
                    domain: "d".into(),
                    label: "e".into(),
                    unit: "J".into(),
                    value: v,
                }]
            })
            .unwrap_or_default(),
        error: e.is_none().then(|| "refused".to_string()),
    };
    let width = sweep::Parameter {
        pointer: "/domains/0/width_nm".into(),
        label: "w".into(),
        unit: "nm",
        now: 10.0,
        integral: false,
    };
    let points = [
        at(20.0, Some(3.0)),
        at(5.0, Some(1.0)),
        at(10.0, None),
        at(15.0, Some(2.5)),
    ];
    let chart = sweep::chart(&width, &points, 0).expect("a chart");
    let plotted = &chart.axes[0].lines[0].points;
    let xs: Vec<f64> = plotted.iter().map(|p| p[0]).collect();
    assert_eq!(xs, [5.0, 10.0, 15.0, 20.0]);
    assert!(
        plotted[1][1].is_nan(),
        "the refused 10 is joined across: {plotted:?}"
    );
    // Drawn, that is a dot at 5 and a stroke from 15 to 20: nothing through 10.
    let d = chart.draw(600.0, 300.0, None);
    assert_eq!(
        (
            d.discs.len(),
            d.strokes.iter().filter(|s| s.width == 2.0).count()
        ),
        (1, 1)
    );
}

/// **A run that panics is that value's refusal**, and the sweep goes on. A sweep hands a scene values
/// nobody wrote; a panic took its thread down with nothing sent, and the editor waited on it forever.
/// Here the parts a scene names panic when read — scene 29's bracket — and the point says so.
#[test]
fn a_run_that_panics_is_a_refused_value() {
    struct Panics;
    impl pantometry_world::Parts for Panics {
        fn bytes(&self, name: &str) -> Result<Vec<u8>, String> {
            panic!("asked for {name}")
        }
    }
    let text = scene("29-");
    let p = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "bracket.initial_c")
        .expect("the bracket's start");
    let points = sweep::run_all(
        &text,
        &p.pointer,
        &[40.0, 41.0],
        &Panics,
        &AtomicBool::new(false),
        |_, _| {},
    );
    assert_eq!(points.len(), 2, "the panic ended the sweep");
    for pt in &points {
        let why = pt.error.as_deref().unwrap_or("");
        assert!(
            why.starts_with("the run panicked: asked for parts/"),
            "{why}"
        );
    }
}

/// **A count is swept between whole numbers**, and a fractional end is refused rather than rounded
/// into a different sweep: `150.5` alone was refused, and `150.5:150.5:1` ran 151.
#[test]
fn a_counts_ends_are_whole() {
    for (from, to) in [(150.5, 150.5), (99.5, 300.0), (100.0, 300.5)] {
        assert!(
            sweep::spaced(from, to, 3, false, true).is_err(),
            "{from} to {to} was swept as a count"
        );
    }
    assert_eq!(
        sweep::spaced(100.0, 300.0, 3, false, true),
        Ok(vec![100.0, 200.0, 300.0])
    );
}

/// **A logarithmic sweep is the same values on every platform.** Its values become literals in the
/// scene, so a platform's `ln` and `exp` differing by an ulp would be a different run of one sweep.
/// They come from the kernel's, whose bits are pinned; these are the four `1:1000:4:log` gives,
/// measured on 2026-11-18. On this machine the platform's functions give the same first interior
/// value, so this pin earns its keep on CI's other platforms rather than here.
#[test]
fn a_logarithmic_sweep_is_the_same_everywhere() {
    let bits: Vec<u64> = sweep::spaced(1.0, 1000.0, 4, true, false)
        .expect("spaced")
        .iter()
        .map(|v| v.to_bits())
        .collect();
    assert_eq!(
        bits,
        [
            0x3ff0000000000000,
            0x4023ffffffffffff,
            0x4058fffffffffffd,
            0x408f400000000000
        ]
    );
}

/// **The chart is the reading against the value, labelled in the scene's own unit.** Three widths,
/// three points, and the horizontal axis says `particle.width_nm / nm` — the number a person typed,
/// not metres.
#[test]
fn the_chart_is_the_reading_against_the_value() {
    let text = scene("28-");
    let width = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.label == "particle.width_nm")
        .expect("a width");
    let points = sweep::run_all(
        &text,
        &width.pointer,
        &[5.0, 10.0, 20.0],
        &files(),
        &AtomicBool::new(false),
        |_, _| {},
    );
    let k = sweep::columns(&points)
        .iter()
        .position(|c| c.label == "<E>")
        .expect("an energy column");
    let chart = sweep::chart(&width, &points, k).expect("a chart");
    assert_eq!((chart.x0, chart.x1), (5.0, 20.0));
    // What is plotted is the column's: each width against its own energy, not another reading's.
    let plotted = &chart.axes[0].lines[0].points;
    let want: Vec<[f64; 2]> = points.iter().map(|p| [p.value, energy(p)]).collect();
    assert_eq!(plotted, &want);
    let d = chart.draw(600.0, 300.0, None);
    assert!(
        d.texts.iter().any(|t| t.text == "particle.width_nm / nm"),
        "{:?}",
        d.texts.iter().map(|t| &t.text).collect::<Vec<_>>()
    );
}

/// **A point's readings are its run's last frame, on a scene whose readings move.** Scene 28's are
/// stationary — an eigenstate's energy, position and norm do not change — so a point read off the
/// first frame would pass every test above. Scene 15's hot spot cools: its peak at the end is not
/// its peak at the start, the point's is the end's, at the seven figures the run file writes.
#[test]
fn a_point_is_its_runs_last_frame() {
    let text = scene("15-");
    let p = sweep::parameters(&text)
        .into_iter()
        .find(|p| p.pointer.starts_with("/domains/") && !p.integral)
        .expect("a number of the block's");
    let point = sweep::run_point(&text, &p.pointer, p.now, &files());
    assert!(point.error.is_none(), "{:?}", point.error);
    let run = viewer_core::Run::from_json(&editor_core::run(&text, &files()).expect("it runs"))
        .expect("it reads");
    let (first, last) = (
        run.frames.first().expect("frames"),
        run.frames.last().expect("frames"),
    );
    let peak = |f: &viewer_core::Frame| {
        f.readings
            .iter()
            .find(|r| r.label == "peak")
            .expect("a peak")
            .value
    };
    assert!(
        (peak(first) - peak(last)).abs() > 1.0,
        "the peak does not move, so this cannot tell the frames apart: {} and {}",
        peak(first),
        peak(last)
    );
    let mine = point
        .readings
        .iter()
        .find(|r| r.label == "peak")
        .expect("a peak")
        .value;
    assert_eq!(format!("{mine:.6e}"), format!("{:.6e}", peak(last)));
}
