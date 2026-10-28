//! **A run with nothing to place is drawn as its readings over time, and the drawing is the run.**
//!
//! Scene `34`, a tablet through the liver into lean tissue and fat, is the case: a compartment
//! has no position, so there is no panel, and the four concentrations against time are what there
//! is to see. Each check here is against the frames themselves or a closed form, never against a
//! second copy of the mapping:
//!
//! - every point of every line is a frame's reading, bit for bit;
//! - a log axis leaves out what has no logarithm instead of pinning it to the floor;
//! - every placed point is inside its plot;
//! - and **the slope the chart draws is the body's slowest mode.** Four days in, one mode is left,
//!   and on a log axis it is a straight line whose slope is its rate. That rate is an eigenvalue
//!   of the scene's own perfusion matrix, found here by bisection on its characteristic
//!   polynomial: not the run, and not the chart.
//!
//! # Not under `wasm32`
//!
//! It reads the shipped scene off a disk.

#![cfg(not(target_family = "wasm"))]

use editor_core::chart::{chart, Axis, Chart};
use viewer_core::Run;

const SCENE: &str = "../pantometry-world/scenes/34-a-tablet-that-meets-the-liver-first.json";

fn text() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join(SCENE);
    std::fs::read_to_string(path).expect("the shipped scene reads")
}

fn ran() -> Run {
    let json = editor_core::run(&text(), &editor_core::OnDisk).expect("the scene runs");
    Run::from_json(&json).expect("the run parses")
}

fn axis<'a>(c: &'a Chart, unit: &str) -> &'a Axis {
    c.axes
        .iter()
        .find(|a| a.unit == unit)
        .unwrap_or_else(|| panic!("no {unit} axis"))
}

/// **Every point is a frame's reading, to the bit, and every reading is on a line.**
#[test]
fn every_point_is_a_reading() {
    let run = ran();
    for log in [false, true] {
        let c = chart(&run, log).expect("a run with readings has a chart");
        let lines: usize = c.axes.iter().map(|a| a.lines.len()).sum();
        assert_eq!(lines, run.frames[0].readings.len(), "one line per reading");
        for a in &c.axes {
            for l in &a.lines {
                assert_eq!(l.points.len(), run.frames.len(), "{} {}", l.domain, l.label);
                for (p, f) in l.points.iter().zip(&run.frames) {
                    let r = f
                        .readings
                        .iter()
                        .find(|r| r.domain == l.domain && r.label == l.label)
                        .expect("the reading is in the frame");
                    assert_eq!(r.unit, a.unit);
                    assert_eq!(p[0].to_bits(), f.t.to_bits());
                    assert_eq!(p[1].to_bits(), r.value.to_bits());
                }
            }
        }
    }
}

/// **A log axis leaves out zero, and a linear one keeps it.** `cleared` is exactly zero before
/// the first step, which has no logarithm; pinned to the bottom of the axis it would be a cliff
/// the run does not have. Everything placed is inside its plot, on either axis.
#[test]
fn a_log_axis_leaves_out_what_it_cannot_place() {
    let run = ran();
    for log in [false, true] {
        let c = chart(&run, log).expect("a chart");
        let mg = axis(&c, "mg");
        let cleared = mg
            .lines
            .iter()
            .find(|l| l.label == "cleared")
            .expect("a cleared line");
        assert_eq!(
            cleared.points[0][1], 0.0,
            "nothing is cleared before a step"
        );
        assert_eq!(mg.place(0.0, log).is_none(), log);

        let d = c.draw(900.0, 600.0, None);
        let drawn: usize = d
            .strokes
            .iter()
            .filter(|s| s.width == 2.0)
            .map(|s| s.points.len())
            .sum();
        let placeable: usize = c
            .axes
            .iter()
            .flat_map(|a| a.lines.iter().flat_map(|l| l.points.iter()))
            .filter(|p| p[1].is_finite() && (!log || p[1] > 0.0))
            .count();
        assert_eq!(
            drawn, placeable,
            "log {log}: every placeable point, and no other"
        );

        for (a, r) in c.axes.iter().zip(c.plots(900.0, 600.0)) {
            for l in &a.lines {
                for p in &l.points {
                    if let Some(u) = a.place(p[1], log) {
                        assert!(
                            (-1e-12..=1.0 + 1e-12).contains(&u),
                            "{} {} at {} is {u} of its axis",
                            l.domain,
                            l.label,
                            p[1]
                        );
                    }
                }
            }
            assert!(r[0] < r[2] && r[1] < r[3], "a plot with room: {r:?}");
        }
    }
}

/// The scene's numbers, off the file: blood volume and renal clearance, and each tissue's
/// volume, flow, partition and clearance, in litres and hours.
fn body() -> (f64, f64, Vec<[f64; 4]>) {
    let v: serde_json::Value = serde_json::from_str(&text()).expect("the scene parses");
    let d = &v["domains"][0];
    let compartment = |label: &str| {
        d["volumes"]
            .as_array()
            .expect("volumes")
            .iter()
            .find(|c| c["label"] == label)
            .unwrap_or_else(|| panic!("no {label}"))
            .clone()
    };
    let blood = compartment("blood");
    let tissues = d["perfusions"]
        .as_array()
        .expect("perfusions")
        .iter()
        .map(|p| {
            let c = compartment(p["tissue"].as_str().expect("a tissue"));
            [
                c["volume_l"].as_f64().expect("a volume"),
                p["flow_l_per_h"].as_f64().expect("a flow"),
                p["partition"].as_f64().expect("a partition"),
                c["clearance_l_per_h"].as_f64().unwrap_or(0.0),
            ]
        })
        .collect();
    (
        blood["volume_l"].as_f64().expect("a blood volume"),
        blood["clearance_l_per_h"].as_f64().unwrap_or(0.0),
        tissues,
    )
}

/// `det(M − λI)` of the perfused body's matrix in concentrations, per hour, by elimination with
/// partial pivoting.
fn characteristic(lambda: f64) -> f64 {
    let (vb, renal, tissues) = body();
    let n = tissues.len() + 1;
    let mut m = vec![vec![0.0f64; n]; n];
    m[0][0] = -renal / vb;
    for (i, [v, q, k, cl]) in tissues.iter().enumerate() {
        let i = i + 1;
        m[0][0] -= q / vb;
        m[0][i] += q / (k * vb);
        m[i][0] = q / v;
        m[i][i] = -(q + cl) / (k * v);
    }
    for (i, row) in m.iter_mut().enumerate() {
        row[i] -= lambda;
    }
    let mut det = 1.0;
    for c in 0..n {
        let p = (c..n)
            .max_by(|a, b| m[*a][c].abs().total_cmp(&m[*b][c].abs()))
            .expect("a pivot");
        if p != c {
            m.swap(p, c);
            det = -det;
        }
        det *= m[c][c];
        let (above, below) = m.split_at_mut(c + 1);
        let pivot = &above[c];
        for row in below {
            let f = row[c] / pivot[c];
            for (x, p) in row[c..].iter_mut().zip(&pivot[c..]) {
                *x -= f * p;
            }
        }
    }
    det
}

/// **The slope the log chart draws is the body's slowest mode.**
///
/// The slowest eigenvalue is bracketed in `[−0.03, −0.015]` per hour and bisected to the last
/// bit; the next is −0.66, so the bracket holds one root. The stepped model's per-step factor for
/// that mode is exactly `1 + hλ` — an eigenvalue of `I + hM` — so the rate a frame-to-frame ratio
/// shows is `−ln(1 + hλ)/h`, which differs from `−λ` by `hλ/2`, 2.9e-5 of it here.
///
/// The rate is read off the **drawing**: two vertices of the blood's stroke a day apart, the last
/// and the one 24 hours before it, turned back into decades per hour by the plot's own size and
/// its axis's decades. Three days in the next mode is already `e^{−46}` of the slowest, so both
/// are on the one line. Two claims, two tolerances:
///
/// - **The drawing is the readings**: against `ln` of the two readings' ratio over their interval,
///   to 1e-9. What is left is the pixels' arithmetic — subtractions of numbers of a few hundred,
///   each good to about 1e-13.
/// - **The readings are the mode**, to what the run's file can carry. A run reaches the editor as
///   JSON, and its writer keeps seven significant figures, so each reading is good to half a unit
///   in the seventh, 5e-7 of itself; two of them over an interval `Δt` put the rate within
///   `1e-6/Δt`, which over a day is 2.0e-6 of it, and the gut is gone, so that rounding is the
///   whole of the bound. It was first written at 1e-9 across the last two frames, and failed at
///   3.2e-6: the file, not the chart, which is the first claim. A day apart rather than one frame
///   is what lets the bound separate the stepped rate from the continuous one, 2.9e-5 apart,
///   by fourteen times.
#[test]
fn the_slope_the_chart_draws_is_the_slowest_mode() {
    let (mut a, mut b) = (-0.03, -0.015);
    assert!(
        characteristic(a).signum() != characteristic(b).signum(),
        "the bracket does not hold a root"
    );
    for _ in 0..200 {
        let m = 0.5 * (a + b);
        if m == a || m == b {
            break;
        }
        if characteristic(m).signum() == characteristic(a).signum() {
            a = m;
        } else {
            b = m;
        }
    }
    let lambda = 0.5 * (a + b);

    let run = ran();
    let steps_per_hour = {
        let v: serde_json::Value = serde_json::from_str(&text()).expect("the scene parses");
        3600.0 / v["window_s"].as_f64().expect("a window")
    };
    let h = 1.0 / steps_per_hour;
    let want = -(1.0 + h * lambda).ln() / h;

    let c = chart(&run, true).expect("a chart");
    let (w, ht) = (900.0, 600.0);
    let d = c.draw(w, ht, None);
    let i = c
        .axes
        .iter()
        .position(|a| a.unit == "mg/L")
        .expect("a mg/L axis");
    let plot = c.plots(w, ht)[i];
    let blood = c.axes[i]
        .lines
        .iter()
        .position(|l| l.label == "blood")
        .expect("a blood line");
    let colour = editor_core::chart::line_colour(blood);
    let stroke = d
        .strokes
        .iter()
        .find(|s| s.width == 2.0 && s.colour == colour && s.points.len() > 2)
        .expect("the blood's stroke");
    let bl = &c.axes[i].lines[blood].points;
    let day = bl
        .iter()
        .rposition(|p| p[0] <= bl[bl.len() - 1][0] - 86_400.0)
        .expect("a frame a day before the last");
    // The blood is zero before the tablet reaches it, and a log axis leaves those out, so the
    // stroke starts that many readings late.
    let late = bl.iter().take_while(|p| p[1] <= 0.0).count();
    assert_eq!(
        late + stroke.points.len(),
        bl.len(),
        "every positive blood reading is drawn"
    );
    let [p, q] = [stroke.points[day - late], *stroke.points.last().unwrap()];
    let ax = &c.axes[i];
    let decades = ax.hi.log10() - ax.lo.log10();
    let fall = (q[1] - p[1]) / (plot[3] - plot[1]) * decades;
    let hours = (q[0] - p[0]) / (plot[2] - plot[0]) * (c.t1 - c.t0) / 3600.0;
    let drawn = fall * std::f64::consts::LN_10 / hours;

    let (u, v) = (bl[day], bl[bl.len() - 1]);
    let interval = (v[0] - u[0]) / 3600.0;
    let read = -(v[1] / u[1]).ln() / interval;
    let carried = 1e-6 / interval / want;
    println!(
        "  slowest mode {:.12}/h, stepped {want:.12}/h, read {read:.12}/h, drawn {drawn:.12}/h",
        -lambda
    );
    println!(
        "  drawn against stepped {:.2e}, against a bound of {carried:.2e}",
        (drawn / want - 1.0).abs()
    );
    assert!(
        (drawn / read - 1.0).abs() < 1e-9,
        "the chart draws the blood falling at {drawn:.12}/h, and its readings fall at {read:.12}"
    );
    assert!(
        (drawn / want - 1.0).abs() < carried,
        "the chart draws the blood falling at {drawn:.12}/h, and the slowest mode is {want:.12}"
    );
}
