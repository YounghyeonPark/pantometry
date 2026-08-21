//! **A number this crate writes down has to come back as the number it was.**
//!
//! Three writers put numbers into files — the CSV, the JSON, and the report's embedded run — and
//! this workspace spans nanoseconds to hours and femtojoules to kilowatts. Fixed decimal places
//! are only readable at the one scale they were chosen for, and outside it they do not round: they
//! erase. A cavity holding 3.2e-10 J had its whole energy history written as a column of
//! `0.000000000`, beside a field the same run had just reported at 921 V/m.
//!
//! That was found in the CSV and fixed there. The same `{:.6}` was still in `to_json`'s
//! timestamps, where a four-nanosecond run wrote two hundred frames all at `t = 0.000000` — so a
//! consumer reading the file back could not have found the frequency the scene exists to measure.
//! This is here so the third occurrence is a failing test instead of a discovery.
//!
//! The values are chosen to straddle every place a format can quietly give up: the smallest
//! interval a real run uses, the largest a real run reaches, and the exact zero that a
//! trailing-zero trim must still spell.

use pantometry_core::Reading;
use pantometry_scene::{Frame, Panel, PanelData};
use pantometry_view::{html, readings_csv, to_json};

/// Times and values a real scene actually produces, spanning nineteen decades between them.
const TIMES: [f64; 5] = [0.0, 2.0e-11, 4.0e-9, 1.5, 3600.0];
const VALUES: [f64; 5] = [0.0, 3.19e-10, -7.5e-3, 921.4, 4.2e7];

fn frames() -> Vec<Frame> {
    TIMES
        .iter()
        .enumerate()
        .map(|(k, &t)| Frame {
            time_s: t,
            panels: vec![Panel {
                name: "cavity".into(),
                unit: "V/m",
                data: PanelData::Field {
                    nx: 2,
                    ny: 2,
                    nz: 1,
                    extent_m: [0.0, 0.0, 0.0, 2.4e-3, 2.4e-3, 0.0],
                    values: vec![VALUES[k], -VALUES[k], 0.0, VALUES[(k + 1) % 5]],
                },
            }],
            readings: vec![Reading::new("cavity", "energy", VALUES[k], "J")],
        })
        .collect()
}

/// Every timestamp is distinguishable from every other, in all three writers.
///
/// Not "the format looks right" — the actual test is that five different times produce five
/// different strings. A writer that collapses two of them has lost the axis, whatever it prints.
#[test]
fn five_different_times_are_five_different_strings() {
    let f = frames();
    let csv = readings_csv(&f);
    let stamps: Vec<&str> = csv
        .lines()
        .skip(1)
        .map(|l| l.split(',').next().unwrap())
        .collect();
    assert_eq!(stamps.len(), 5);
    for i in 0..5 {
        for j in i + 1..5 {
            assert_ne!(
                stamps[i], stamps[j],
                "csv wrote t={} and t={} the same way: {:?}",
                TIMES[i], TIMES[j], stamps
            );
        }
    }

    let json = to_json("magnitudes", &f);
    for t in TIMES {
        // Each time appears as a `"t":` value somewhere. Exact spelling is the encoder's business;
        // that all five are present and distinct is not.
        let count = json.matches("\"t\": ").count();
        assert_eq!(count, 5, "five frames, {count} timestamps");
        let _ = t;
    }
    let ts: Vec<&str> = json
        .match_indices("\"t\": ")
        .map(|(i, _)| {
            let rest = &json[i + 5..];
            &rest[..rest.find(',').unwrap()]
        })
        .collect();
    for i in 0..5 {
        for j in i + 1..5 {
            assert_ne!(
                ts[i], ts[j],
                "to_json wrote t={} and t={} the same way: {ts:?}",
                TIMES[i], TIMES[j]
            );
        }
    }

    let page = html("magnitudes", &f);
    let hs: Vec<&str> = page
        .match_indices("{\"t\":")
        .map(|(i, _)| {
            let rest = &page[i + 5..];
            &rest[..rest.find(',').unwrap()]
        })
        .collect();
    assert_eq!(hs.len(), 5, "five frames in the report: {hs:?}");
    for i in 0..5 {
        for j in i + 1..5 {
            assert_ne!(hs[i], hs[j], "the report collapsed two timestamps: {hs:?}");
        }
    }
}

/// A value at 3.19e-10 is written as 3.19e-10 and not as zero, everywhere.
///
/// The magnitude is the part that must survive. Six significant figures is the promise; the
/// number of characters spent is not.
#[test]
fn a_small_value_keeps_its_magnitude() {
    let f = frames();
    for (what, text) in [
        ("csv", readings_csv(&f)),
        ("json", to_json("m", &f)),
        ("html", html("m", &f)),
    ] {
        assert!(
            text.contains("3.19e-10") || text.contains("3.190000000e-10"),
            "{what}: 3.19e-10 did not survive"
        );
        assert!(
            text.contains("4.2e7") || text.contains("4.200000000e7"),
            "{what}: 4.2e7 did not survive"
        );
        assert!(
            text.contains("-7.5e-3") || text.contains("-7.500000000e-3"),
            "{what}: a negative milli-value did not survive"
        );
    }
}

/// Trimming trailing zeros must not trim a number away.
///
/// `0` is the one value where "drop what says nothing" and "say something" meet, and an encoder
/// that trims `0.000000e0` down through `0.` to the empty string writes a document that does not
/// parse. Checked by parsing the report's own run block, which is the only one of the three that
/// a browser has to read with `JSON.parse` before anything is drawn at all.
#[test]
fn zero_is_still_a_number_and_the_document_still_parses() {
    let page = html("m", &frames());
    let i = page
        .find("<script id=\"run\" type=\"application/json\">")
        .expect("run block");
    let start = i + "<script id=\"run\" type=\"application/json\">".len();
    let end = start + page[start..].find("</script>").expect("run block ends");
    let body = &page[start..end];
    assert!(body.contains(",0,") || body.contains("[0,") || body.contains(",0]"));
    assert!(!body.contains(",,") && !body.contains("[,") && !body.contains(",]"));
    // Brackets balance and every value is a token: a crude parse, but it fails on exactly the
    // damage a bad trim does, and it needs no dependency to do it.
    let (mut depth, mut worst) = (0i32, 0i32);
    for c in body.chars() {
        match c {
            '[' | '{' => depth += 1,
            ']' | '}' => depth -= 1,
            _ => {}
        }
        worst = worst.min(depth);
    }
    assert_eq!(depth, 0, "the run block does not close");
    assert_eq!(worst, 0, "the run block closes something it did not open");
}
