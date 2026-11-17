//! **A result is a number a closed form already knows.**
//!
//! `editor_core::post` reads values, profiles, cuts and statistics off a run's fields. Every
//! check here is against something that does not come from `post`: a sample the run holds, a
//! linear function the interpolation must reproduce, the error the midpoint and trapezoid rules
//! leave on `x^2` to the last digit, and an insulated block whose mean temperature its initial
//! condition fixes for the whole run.

use editor_core::chart::{self, Across, Line};
use editor_core::post::{self, Grid};
use viewer_core::{Frame, Lattice, Panel, Placed, Run};

/// A grid that is not a cube and not at the origin, so an axis swapped or an offset dropped
/// shows: 3 x 5 x 7 samples over a box from (2, -1, 0.5) mm to (11, 4, 7.5) mm.
const COUNTS: [usize; 3] = [3, 5, 7];
const EXTENT: [f64; 6] = [0.002, -0.001, 0.0005, 0.011, 0.004, 0.0075];

fn field(lattice: Lattice, counts: [usize; 3], extent: [f64; 6], values: Vec<f64>) -> Panel {
    Panel::Field {
        name: "part".to_string(),
        unit: "K".to_string(),
        place: Placed::default(),
        nx: counts[0],
        ny: counts[1],
        nz: counts[2],
        extent_m: Some(extent),
        lattice,
        values,
    }
}

/// Where sample `i` of `n` sits along `[lo, hi]`, written out here rather than asked of `post`.
fn position(lattice: Lattice, lo: f64, hi: f64, i: usize, n: usize) -> f64 {
    let f = match (lattice, n) {
        (_, 1) => 0.5,
        (Lattice::Nodal, _) => i as f64 / (n - 1) as f64,
        (Lattice::Centred, _) => (i as f64 + 0.5) / n as f64,
    };
    lo + f * (hi - lo)
}

/// `f` sampled on the grid, x fastest.
fn sampled(
    lattice: Lattice,
    counts: [usize; 3],
    extent: [f64; 6],
    f: impl Fn([f64; 3]) -> f64,
) -> Vec<f64> {
    let mut out = Vec::new();
    for k in 0..counts[2] {
        for j in 0..counts[1] {
            for i in 0..counts[0] {
                out.push(f([
                    position(lattice, extent[0], extent[3], i, counts[0]),
                    position(lattice, extent[1], extent[4], j, counts[1]),
                    position(lattice, extent[2], extent[5], k, counts[2]),
                ]));
            }
        }
    }
    out
}

/// A field linear in all three, around 300 K, with a different slope on each axis.
fn linear(p: [f64; 3]) -> f64 {
    300.0 + 1.0e3 * p[0] - 2.0e3 * p[1] + 5.0e2 * p[2]
}

const LATTICES: [Lattice; 2] = [Lattice::Nodal, Lattice::Centred];

/// Points spread through the box, inside the outermost samples, by a fixed recipe.
fn interior(lattice: Lattice) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for s in 0..40 {
        let f = |a: usize, m: usize| {
            let n = COUNTS[a];
            let first = position(lattice, EXTENT[a], EXTENT[a + 3], 0, n);
            let last = position(lattice, EXTENT[a], EXTENT[a + 3], n - 1, n);
            first + ((s * m) % 97) as f64 / 96.0 * (last - first)
        };
        out.push([f(0, 13), f(1, 29), f(2, 41)]);
    }
    out
}

#[test]
fn a_probe_at_a_sample_reads_the_sample_bit_for_bit() {
    for lattice in LATTICES {
        // Not linear, so an interpolation that leaned on a neighbour would land somewhere else.
        let values: Vec<f64> = (0..105)
            .map(|n| 290.0 + ((n * 31) % 13) as f64 * 0.37)
            .collect();
        let panel = field(lattice, COUNTS, EXTENT, values.clone());
        let g = Grid::of(&panel).expect("a field");
        for k in 0..7 {
            for j in 0..5 {
                for i in 0..3 {
                    let p = [
                        position(lattice, EXTENT[0], EXTENT[3], i, 3),
                        position(lattice, EXTENT[1], EXTENT[4], j, 5),
                        position(lattice, EXTENT[2], EXTENT[5], k, 7),
                    ];
                    let want = values[i + 3 * (j + 5 * k)];
                    let got = g.at(p).expect("inside");
                    assert_eq!(got.to_bits(), want.to_bits(), "{lattice:?} [{i},{j},{k}]");
                }
            }
        }
    }
}

/// Eight products of a weight at most one and a value near 300, summed: a few units in the last
/// place of 300 each, so 1e-11 K is a bound with two orders to spare and nothing physical in it.
const ROUNDING: f64 = 1e-11;

#[test]
fn a_linear_field_is_exact_between_samples() {
    for lattice in LATTICES {
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        for p in interior(lattice) {
            let got = g.at(p).expect("inside");
            assert!(
                (got - linear(p)).abs() <= ROUNDING,
                "{lattice:?} at {p:?}: {got} against {}",
                linear(p)
            );
        }
    }
}

#[test]
fn past_the_last_cell_centre_the_value_is_that_cells() {
    let panel = field(
        Lattice::Centred,
        COUNTS,
        EXTENT,
        sampled(Lattice::Centred, COUNTS, EXTENT, linear),
    );
    let g = Grid::of(&panel).expect("a field");
    let (y, z) = (0.0012, 0.004);
    // Between either wall and the centre nearest it: that centre's value, not a line run past it.
    let first_x = position(Lattice::Centred, EXTENT[0], EXTENT[3], 0, 3);
    let last_x = position(Lattice::Centred, EXTENT[0], EXTENT[3], 2, 3);
    for (wall, centre) in [(EXTENT[0], first_x), (EXTENT[3], last_x)] {
        for x in [wall, 0.5 * (wall + centre)] {
            let near = g.at([x, y, z]).expect("inside the box");
            let at_centre = g.at([centre, y, z]).expect("inside");
            assert_eq!(near.to_bits(), at_centre.to_bits(), "x = {x}");
            assert!(
                (near - linear([x, y, z])).abs() > 1e-3,
                "it is not extrapolated at x = {x}"
            );
        }
    }
}

#[test]
fn outside_the_box_there_is_no_value() {
    for lattice in LATTICES {
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        let c = g.centre();
        for axis in 0..3 {
            for (edge, step) in [(EXTENT[axis], -1e-9), (EXTENT[axis + 3], 1e-9)] {
                let mut p = c;
                p[axis] = edge;
                assert!(g.at(p).is_some(), "{lattice:?}: on the wall is inside");
                p[axis] = edge + step;
                assert_eq!(g.at(p), None, "{lattice:?}: past the wall along {axis}");
            }
        }
        // A line that leaves the box breaks there rather than reading zero.
        let line = g.line([EXTENT[0], c[1], c[2]], [EXTENT[3] + 0.009, c[1], c[2]], 11);
        assert!(line[0][1].is_finite() && line[10][1].is_nan(), "{line:?}");
    }
}

#[test]
fn a_hole_is_not_averaged_over() {
    let lattice = Lattice::Centred;
    let mut values = sampled(lattice, COUNTS, EXTENT, linear);
    let hole = [1, 2, 3];
    let h = hole[0] + 3 * (hole[1] + 5 * hole[2]);
    let held = values[h];
    values[h] = f64::NAN;
    let panel = field(lattice, COUNTS, EXTENT, values.clone());
    let g = Grid::of(&panel).expect("a field");
    let at = |i: usize, a: usize| position(lattice, EXTENT[a], EXTENT[a + 3], i, COUNTS[a]);

    // Anywhere the hole has a weight: nothing.
    let beside = [
        0.5 * (at(1, 0) + at(2, 0)),
        0.5 * (at(2, 1) + at(3, 1)),
        0.5 * (at(3, 2) + at(4, 2)),
    ];
    assert_eq!(g.at(beside), None, "a point beside the hole");
    // On a neighbour exactly: the neighbour, since the hole's weight is zero. Sample [1, 1, 3]
    // has the hole at [1, 2, 3] as its upper neighbour along y, so the hole is one of the eight
    // corners with a weight of zero — which is the case the skip is for. The first version of
    // this probed [2, 2, 3], the last sample along x, where the hole is not a corner at all.
    let neighbour = [at(1, 0), at(1, 1), at(3, 2)];
    assert_eq!(g.at(neighbour), Some(values[1 + 3 * (1 + 5 * 3)]));
    // The same in a slice: the plane through z sample 2 has the hole's plane as its zero-weight
    // neighbour, and reads its own row whole.
    let s = g.slice(2, at(2, 2)).expect("inside");
    assert!(s.values.iter().all(|v| v.is_finite()), "{:?}", s.values);
    let s = g.slice(2, at(3, 2)).expect("inside");
    assert!(s.value(1, 2).is_nan(), "the plane through the hole has it");

    let whole = Grid::of(&field(lattice, COUNTS, EXTENT, {
        let mut v = values.clone();
        v[h] = held;
        v
    }))
    .map(|g| g.region(None).expect("stats"))
    .expect("a field");
    let s = g.region(None).expect("stats");
    assert_eq!((s.samples, s.holes), (104, 1));
    let cell = (0.009 / 3.0) * (0.005 / 5.0) * (0.007 / 7.0);
    // Two sums of 105 terms differenced: each is good to about 1e-14 of itself, which is 1e-12
    // of one cell's share, so a part in a billion is rounding and the whole cell is not.
    let gone = whole.measure - s.measure;
    assert!(
        (gone / cell - 1.0).abs() <= 1e-9,
        "the hole's volume is not counted: {gone}"
    );
    let lost = whole.integral - s.integral;
    assert!((lost / (held * cell) - 1.0).abs() <= 1e-9, "{lost}");

    let all_holes = field(lattice, COUNTS, EXTENT, vec![f64::NAN; 105]);
    assert_eq!(
        Grid::of(&all_holes).expect("a field").region(None),
        None,
        "no mean of nothing"
    );
}

#[test]
fn a_linear_fields_integral_is_its_centre_value_times_the_volume() {
    let volume = 0.009 * 0.005 * 0.007;
    for lattice in LATTICES {
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        let s = g.region(None).expect("stats");
        let centre = linear(g.centre());
        assert_eq!(s.dims, 3);
        assert!(
            (s.measure / volume - 1.0).abs() <= 1e-14,
            "{lattice:?}: {}",
            s.measure
        );
        assert!(
            (s.mean - centre).abs() <= ROUNDING,
            "{lattice:?}: {} against {centre}",
            s.mean
        );
        assert!(
            (s.integral / (centre * volume) - 1.0).abs() <= 1e-13,
            "{lattice:?}"
        );
        assert_eq!(s.samples, 105);
        // The extremes of a linear field are at opposite corners of the samples.
        let lo = [
            position(lattice, EXTENT[0], EXTENT[3], 0, 3),
            position(lattice, EXTENT[1], EXTENT[4], 4, 5),
            position(lattice, EXTENT[2], EXTENT[5], 0, 7),
        ];
        assert_eq!(s.min_at, lo, "{lattice:?}");
        assert_eq!(s.min.to_bits(), linear(lo).to_bits());
    }
}

/// **The rule, by the error it leaves.** Both rules are exact for a linear field, so the linear
/// test cannot tell them apart. On `x^2` over `[0, L]` each leaves an error known to the last
/// digit: the midpoint rule on `n` cells is `L^3/3 - L^3/(12 n^2)` and the trapezoid rule on `m`
/// intervals `L^3/3 + L^3/(6 m^2)`. Opposite signs, and a factor of two apart.
#[test]
fn each_lattice_is_integrated_by_its_own_rule() {
    let length = 0.02;
    let cube = length * length * length;
    for n in [4usize, 8, 16] {
        // A line: one sample and no thickness across, so the measure is a length.
        let extent = [0.0, 0.0, 0.0, length, 0.0, 0.0];
        for (lattice, want) in [
            (
                Lattice::Centred,
                cube / 3.0 - cube / (12.0 * (n * n) as f64),
            ),
            (
                Lattice::Nodal,
                cube / 3.0 + cube / (6.0 * ((n - 1) * (n - 1)) as f64),
            ),
        ] {
            let counts = [n, 1, 1];
            let panel = field(
                lattice,
                counts,
                extent,
                sampled(lattice, counts, extent, |p| p[0] * p[0]),
            );
            let s = Grid::of(&panel)
                .expect("a field")
                .region(None)
                .expect("stats");
            assert_eq!(s.dims, 1);
            assert!((s.measure - length).abs() <= 1e-17);
            assert!(
                (s.integral / want - 1.0).abs() <= 1e-12,
                "{lattice:?} n={n}: {} against {want}",
                s.integral
            );
            // **And the mean is the integral over the length, not the samples' average.** On
            // samples spaced evenly with a half weight at each end the two differ: at n = 4 the
            // average of x^2 is 0.389 L^2 and the mean by length 0.352 L^2. Every other mean in
            // these tests is of a linear field or of equal weights, where they agree.
            let by_samples = (0..n)
                .map(|i| {
                    let x = position(lattice, 0.0, length, i, n);
                    x * x
                })
                .sum::<f64>()
                / n as f64;
            assert!(
                (s.mean / (want / length) - 1.0).abs() <= 1e-12,
                "{lattice:?} n={n}: mean {} against {}",
                s.mean,
                want / length
            );
            if lattice == Lattice::Nodal {
                assert!(
                    (s.mean / by_samples - 1.0).abs() > 0.01,
                    "the case does not tell the two means apart"
                );
            }
        }
    }
}

/// **A box selects the samples inside it, and an empty box selects nothing.** On the linear field
/// the samples of a box are known by index, so its count, its volume and its mean are closed forms:
/// cell averages give every sample the same share, and the mean of a linear field over a box of
/// whole cells is its value at the middle of those cells.
#[test]
fn a_box_reads_the_samples_inside_it() {
    let lattice = Lattice::Centred;
    let panel = field(
        lattice,
        COUNTS,
        EXTENT,
        sampled(lattice, COUNTS, EXTENT, linear),
    );
    let g = Grid::of(&panel).expect("a field");
    let at = |i: usize, a: usize| position(lattice, EXTENT[a], EXTENT[a + 3], i, COUNTS[a]);
    // x samples 1..=2, y 0..=3, z 2..=4: a 2 x 4 x 3 block of cells, given far corner first.
    let within = [at(2, 0), at(3, 1), at(4, 2), at(1, 0), at(0, 1), at(2, 2)];
    let s = g.region(Some(within)).expect("the box holds cells");
    assert_eq!((s.samples, s.holes), (24, 0));
    let cell = (0.009 / 3.0) * (0.005 / 5.0) * (0.007 / 7.0);
    assert!((s.measure / (24.0 * cell) - 1.0).abs() <= 1e-12);
    let middle = [
        0.5 * (at(1, 0) + at(2, 0)),
        0.5 * (at(0, 1) + at(3, 1)),
        0.5 * (at(2, 2) + at(4, 2)),
    ];
    assert!(
        (s.mean - linear(middle)).abs() <= ROUNDING,
        "{} against {}",
        s.mean,
        linear(middle)
    );
    assert_eq!(s.max_at, [at(2, 0), at(0, 1), at(4, 2)]);
    // A box between samples, and a box outside the domain: no cells, so no statistics.
    let between = 0.5 * (at(0, 0) + at(1, 0));
    let thin = [
        between,
        EXTENT[1],
        EXTENT[2],
        between + 1e-6,
        EXTENT[4],
        EXTENT[5],
    ];
    assert_eq!(g.region(Some(thin)), None);
    let outside = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0];
    assert_eq!(g.region(Some(outside)), None);
}

#[test]
fn a_sheet_is_integrated_over_its_area_and_any_depth_is_on_it() {
    let extent = [0.0, 0.0, 0.003, 0.004, 0.006, 0.003];
    let counts = [4, 6, 1];
    for lattice in LATTICES {
        let panel = field(
            lattice,
            counts,
            extent,
            sampled(lattice, counts, extent, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        let s = g.region(None).expect("stats");
        assert_eq!(s.dims, 2);
        assert!((s.measure / (0.004 * 0.006) - 1.0).abs() <= 1e-14);
        assert!(
            (s.mean - linear(g.centre())).abs() <= ROUNDING,
            "{lattice:?}"
        );
        let mut p = g.centre();
        let on = g.at(p).expect("on the sheet");
        p[2] = 1.0;
        assert_eq!(
            g.at(p),
            Some(on),
            "a sheet has no thickness to be outside of"
        );
    }
}

#[test]
fn a_slice_through_a_row_of_samples_is_that_row() {
    for lattice in LATTICES {
        let values: Vec<f64> = (0..105)
            .map(|n| 290.0 + ((n * 31) % 13) as f64 * 0.37)
            .collect();
        let panel = field(lattice, COUNTS, EXTENT, values.clone());
        let g = Grid::of(&panel).expect("a field");
        for (axis, k) in [(0usize, 1usize), (1, 3), (2, 5)] {
            let at = position(lattice, EXTENT[axis], EXTENT[axis + 3], k, COUNTS[axis]);
            let s = g.slice(axis, at).expect("inside");
            let [ua, va] = s.axes;
            assert!(ua < va && ua != axis && va != axis);
            assert_eq!(s.counts, [COUNTS[ua], COUNTS[va]]);
            for j in 0..s.counts[1] {
                for i in 0..s.counts[0] {
                    let mut idx = [0; 3];
                    idx[ua] = i;
                    idx[va] = j;
                    idx[axis] = k;
                    let want = values[idx[0] + 3 * (idx[1] + 5 * idx[2])];
                    assert_eq!(
                        s.value(i, j).to_bits(),
                        want.to_bits(),
                        "{lattice:?} {axis}"
                    );
                }
            }
        }
        // Between two planes of a linear field, the linear field.
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        let z = 0.5
            * (position(lattice, EXTENT[2], EXTENT[5], 2, 7)
                + position(lattice, EXTENT[2], EXTENT[5], 3, 7));
        let s = g.slice(2, z).expect("inside");
        for j in 0..5 {
            for i in 0..3 {
                let p = [
                    position(lattice, EXTENT[0], EXTENT[3], i, 3),
                    position(lattice, EXTENT[1], EXTENT[4], j, 5),
                    z,
                ];
                assert!((s.value(i, j) - linear(p)).abs() <= ROUNDING);
                assert_eq!(s.point(p[0], p[1]), p);
            }
        }
        assert_eq!(
            g.slice(2, EXTENT[5] + 1e-9),
            None,
            "a plane outside the box"
        );
    }
}

/// **Each pixel's centre on its sample.** For cell averages that is the plane's own extent; for
/// samples on the faces it is half a spacing wider each way, and drawn into the extent itself a
/// nodal field's picture is up to half a pixel out.
#[test]
fn a_slices_image_puts_each_pixel_on_its_sample() {
    for lattice in LATTICES {
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let s = Grid::of(&panel)
            .expect("a field")
            .slice(2, 0.004)
            .expect("inside");
        let e = s.image_extent();
        for (a, axis) in s.axes.iter().enumerate() {
            let n = s.counts[a];
            let pixel = (e[a + 2] - e[a]) / n as f64;
            for i in 0..n {
                let centre = e[a] + (i as f64 + 0.5) * pixel;
                let sample = position(lattice, EXTENT[*axis], EXTENT[*axis + 3], i, n);
                assert!(
                    (centre - sample).abs() <= 1e-15,
                    "{lattice:?} axis {axis} pixel {i}"
                );
            }
        }
    }
    // Rows top first, a hole transparent.
    let mut values = sampled(Lattice::Centred, COUNTS, EXTENT, linear);
    values[2 + 3 * (4 + 5 * 3)] = f64::NAN;
    let panel = field(Lattice::Centred, COUNTS, EXTENT, values);
    let s = Grid::of(&panel)
        .expect("a field")
        .slice(2, position(Lattice::Centred, EXTENT[2], EXTENT[5], 3, 7))
        .expect("inside");
    let colouring = editor_core::Colouring::of("K", &s.values, Some((290.0, 320.0)));
    let rgba = s.rgba(&colouring);
    assert_eq!(rgba.len(), 3 * 5 * 4);
    // Sample (2, 4) is the top row's last pixel.
    assert_eq!(&rgba[8..12], &[0, 0, 0, 0]);
    assert!(rgba.chunks(4).filter(|p| p[3] == 255).count() == 14);
}

#[test]
fn a_line_is_read_at_evenly_spaced_points_from_end_to_end() {
    for lattice in LATTICES {
        let panel = field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, linear),
        );
        let g = Grid::of(&panel).expect("a field");
        let inner = interior(lattice);
        let (a, b) = (inner[3], inner[17]);
        let n = g.points_along(a, b);
        let line = g.line(a, b, n);
        assert_eq!(line.len(), n);
        let length = post::distance(a, b);
        assert_eq!(line[0][0], 0.0);
        assert_eq!(line[n - 1][0], length);
        for (k, [s, v]) in line.iter().enumerate() {
            let f = k as f64 / (n - 1) as f64;
            assert!((s - length * f).abs() <= 1e-18);
            let p: [f64; 3] = std::array::from_fn(|i| a[i] * (1.0 - f) + b[i] * f);
            assert!((v - linear(p)).abs() <= ROUNDING, "{lattice:?} point {k}");
        }
    }
    // A line from one sample to another reads them at its ends, to the bit. It does not tell the
    // two ways of writing the line's points apart: where one misses its end by an ulp, the value
    // read there is the same.
    for lattice in LATTICES {
        let values: Vec<f64> = (0..105)
            .map(|n| 290.0 + ((n * 31) % 13) as f64 * 0.37)
            .collect();
        let panel = field(lattice, COUNTS, EXTENT, values.clone());
        let g = Grid::of(&panel).expect("a field");
        let s = |i: usize, j: usize, k: usize| [g.sample(0, i), g.sample(1, j), g.sample(2, k)];
        for ((a, ia), (b, ib)) in [
            ((s(0, 1, 2), [0, 1, 2]), (s(2, 3, 6), [2, 3, 6])),
            ((s(1, 4, 5), [1, 4, 5]), (s(2, 0, 1), [2, 0, 1])),
        ] {
            let line = g.line(a, b, 7);
            let at = |i: [usize; 3]| values[i[0] + 3 * (i[1] + 5 * i[2])];
            assert_eq!(line[0][1].to_bits(), at(ia).to_bits(), "{lattice:?} first");
            assert_eq!(line[6][1].to_bits(), at(ib).to_bits(), "{lattice:?} last");
        }
    }
    // Four points a cell along the axis crossed most, and both ends.
    let panel = field(Lattice::Centred, COUNTS, EXTENT, vec![300.0; 105]);
    let g = Grid::of(&panel).expect("a field");
    assert_eq!(g.points_along([0.002, 0.0, 0.004], [0.011, 0.0, 0.004]), 13);
    assert_eq!(g.points_along([0.004, 0.0, 0.004], [0.004, 0.0, 0.004]), 2);
}

#[test]
fn a_history_has_one_point_a_frame_and_a_gap_where_the_field_is_not() {
    let lattice = Lattice::Centred;
    let at = |shift: f64| {
        field(
            lattice,
            COUNTS,
            EXTENT,
            sampled(lattice, COUNTS, EXTENT, |p| linear(p) + shift),
        )
    };
    let frame = |t: f64, panels: Vec<Panel>| Frame {
        t,
        panels,
        readings: Vec::new(),
    };
    let run = Run {
        title: "three frames".to_string(),
        frames: vec![
            frame(0.0, vec![at(0.0)]),
            frame(1.0, Vec::new()),
            frame(2.0, vec![at(5.0)]),
        ],
    };
    let p = interior(lattice)[7];
    let h = post::history(&run, "part", p);
    assert_eq!(h.len(), 3);
    assert_eq!(h[0][0], 0.0);
    assert!((h[0][1] - linear(p)).abs() <= ROUNDING);
    assert!(h[1][1].is_nan(), "a frame without the field is a gap");
    assert!((h[2][1] - (linear(p) + 5.0)).abs() <= ROUNDING);
    let r = post::region_history(&run, "part", None);
    assert!(r[1].1.is_none() && r[2].1.is_some());
    assert_eq!(post::fields(&run), ["part"]);
}

#[test]
fn bodies_are_not_a_field_and_say_so() {
    let bodies = Panel::Points {
        name: "crambin".to_string(),
        unit: "m".to_string(),
        place: Placed::default(),
        boxed: false,
        bounds: [0.0; 6],
        positions: vec![0.0; 3],
        values: vec![1.0],
        labels: Vec::new(),
        bonds: Vec::new(),
    };
    let why = Grid::of(&bodies).expect_err("not a field");
    assert!(why.starts_with("crambin is"), "{why}");
    let mut old = field(Lattice::Centred, COUNTS, EXTENT, vec![1.0; 105]);
    if let Panel::Field { extent_m, .. } = &mut old {
        *extent_m = None;
    }
    let why = Grid::of(&old).expect_err("no extent");
    assert!(why.contains("does not say where its samples are"), "{why}");
    let short = field(Lattice::Centred, COUNTS, EXTENT, vec![1.0; 104]);
    assert!(Grid::of(&short).is_err());
    // **A field it cannot read is still a field of the run**, listed so the panel can say why,
    // rather than left out so the panel says the run has none.
    let run = Run {
        title: "a short field".to_string(),
        frames: vec![Frame {
            t: 0.0,
            panels: vec![short],
            readings: Vec::new(),
        }],
    };
    assert_eq!(post::fields(&run), ["part"]);
    let why = post::grid(&run.frames[0], "part").expect_err("unreadable");
    assert!(why.contains("104 values for a 3 x 5 x 7 grid"), "{why}");
}

/// Scene 15's text with its hot cell moved to `at` and, if `cooled` says so, its `x-min` face losing
/// heat at that coefficient to air at the temperature the block starts at.
fn scene_15(at: [usize; 3], cooled: Option<f64>) -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../pantometry-world/scenes/15-a-hot-spot-in-a-block.json");
    let text = std::fs::read_to_string(&path).expect("scene 15 reads");
    let mut scene: serde_json::Value = serde_json::from_str(&text).expect("scene 15 parses");
    let block = &mut scene["domains"][0];
    block["hot_spot"]["at"] = serde_json::json!(at);
    if let Some(h) = cooled {
        block["cooling"] = serde_json::json!([{
            "face": "x-min", "ambient_c": 20.0, "convection_w_per_m2_k": h, "area_cm2": 0.81
        }]);
    }
    scene.to_string()
}

/// Scene `text`'s block, every frame: the statistics of the whole domain.
fn block_over_time(text: &str) -> Vec<(f64, post::Stats)> {
    let json = editor_core::run(text, &editor_core::OnDisk).expect("the scene runs");
    let run = Run::from_json(&json).expect("a run");
    assert_eq!(post::fields(&run), ["block"]);
    let history = post::region_history(&run, "block", None);
    assert_eq!(history.len(), run.frames.len());
    assert!(history.len() > 80, "{} frames", history.len());
    history
        .into_iter()
        .map(|(t, s)| (t, s.expect("every frame has the block")))
        .collect()
}

/// **An insulated block keeps its mean, and a leaking one does not.** Scene 15 is aluminium at
/// 20 C — constant conductivity, heat capacity and density — with one cell of 729 at 60 K above
/// it and no face cooled, so its mean is 293.15 + 60/729 K in every frame, and the hottest cell
/// can only cool and the coldest only warm.
///
/// The run writes each value to seven figures, 1e-4 K at 293 K, so the mean of what it wrote is
/// within 5e-5 K of the true one: that is the tolerance, and it is the file's. Rounding to a fixed
/// grid keeps order, so the maximum principle needs no slack at all.
///
/// **As shipped, the hot cell is in the middle, and that cannot tell insulated from cooled.** Over
/// the run's 6 ms the heat barely reaches the walls: `numerics-reviewer` found the walls could be
/// held at 20 C and the written mean still stay inside 5e-5 K. Measured here, cooling all six
/// faces at 25 W/m^2/K left it unchanged at 9.5e-6 K. So the cell that matters is in a corner,
/// beside three walls. Insulated, its mean held to 3.2e-6 K. With its `x-min` face cooled at
/// 1e4 W/m^2/K it fell 1.6e-3 K, 33 times the tolerance, and the second half of this test asks
/// for that fall.
#[test]
fn an_insulated_block_keeps_its_mean_and_a_leaking_one_does_not() {
    let mean = 293.15 + 60.0 / 729.0;
    let volume = 0.009f64 * 0.009 * 0.009;
    let written = 5e-5;
    for at in [[4, 4, 4], [0, 0, 0]] {
        let history = block_over_time(&scene_15(at, None));
        let mut hottest = f64::INFINITY;
        for (t, s) in &history {
            assert_eq!((s.samples, s.holes, s.dims), (729, 0, 3));
            assert!((s.measure / volume - 1.0).abs() <= 1e-12, "t = {t}");
            assert!(
                (s.mean - mean).abs() <= written,
                "hot cell at {at:?}, t = {t}: mean {} against {mean}",
                s.mean
            );
            // The integral against the initial condition's heat, not against mean times measure,
            // which is how the mean is defined and so cannot disagree with it.
            assert!(
                (s.integral - mean * volume).abs() <= written * volume,
                "hot cell at {at:?}, t = {t}: integral {}",
                s.integral
            );
            assert!(
                s.max <= hottest,
                "hot cell at {at:?}, t = {t}: the hottest cell warmed to {}",
                s.max
            );
            assert!(
                s.min >= 293.15,
                "hot cell at {at:?}, t = {t}: a cell fell to {}, below where all began",
                s.min
            );
            hottest = s.max;
        }
        // At the start the hot cell is the hottest, and it is where the scene put it.
        let first = history[0].1;
        let at_m = at.map(|i| (i as f64 + 0.5) * 1e-3);
        assert_eq!(first.max, 353.15);
        for a in 0..3 {
            assert!(
                (first.max_at[a] - at_m[a]).abs() <= 1e-15,
                "{:?}",
                first.max_at
            );
        }
        // And it has spread: by the end the peak is well below where it started.
        assert!(hottest < 320.0, "hot cell at {at:?}: {hottest}");
    }

    // The same check, on a block that leaks: it has to see the heat leave.
    let leaking = block_over_time(&scene_15([0, 0, 0], Some(1e4)));
    let (t, last) = leaking.last().expect("frames");
    assert!(
        mean - last.mean > 1e-3,
        "t = {t}: a block losing heat through a face kept its mean at {} against {mean}",
        last.mean
    );
}

#[test]
fn a_profile_is_charted_against_distance_on_an_axis_that_fits_it() {
    let line = Line {
        domain: "block".to_string(),
        label: "along x".to_string(),
        points: vec![[0.0, 293.0], [0.004, 353.0], [0.009, f64::NAN]],
    };
    let c = chart::graph(Across::Distance, "K", vec![line]).expect("a chart");
    assert_eq!((c.x0, c.x1), (0.0, 0.004));
    let axis = &c.axes[0];
    assert!(
        (axis.lo - 290.0).abs() < 1e-9 && (axis.hi - 356.0).abs() < 1e-9,
        "{axis:?}"
    );
    let d = c.draw(400.0, 200.0, None);
    assert!(
        d.texts.iter().any(|t| t.text == "distance / mm"),
        "{:?}",
        d.texts
    );

    // A line with no range is drawn in the middle of an axis that has some.
    let flat = Line {
        domain: "block".to_string(),
        label: "mean".to_string(),
        points: vec![[0.0, 293.2], [1.0, 293.2]],
    };
    let c = chart::graph(Across::Time, "K", vec![flat]).expect("a chart");
    assert!(c.axes[0].hi > c.axes[0].lo);
    let middle = c.axes[0].place(293.2, false).expect("placed");
    assert!((middle - 0.5).abs() <= 1e-12, "{middle}");
    assert_eq!(chart::graph(Across::Time, "K", Vec::new()), None);
}

/// **A profile flat to its rounding is drawn flat.** Scene 24's line through its module is
/// 488.5735 K to within 1.7e-13 — the interpolation's rounding of a field flat in the file — and
/// an axis fitted to that range labelled its ticks to twelve decimals, wider than the room the
/// plot leaves for them.
#[test]
fn a_profile_flat_to_its_rounding_is_drawn_flat() {
    let flat = 488.573_45;
    let line = Line {
        domain: "module".to_string(),
        label: "along the line".to_string(),
        points: vec![[0.0, flat], [0.006, flat + 1.7e-13], [0.012, flat]],
    };
    let c = chart::graph(Across::Distance, "K", vec![line]).expect("a chart");
    let axis = &c.axes[0];
    assert!(axis.hi - axis.lo > 0.5, "an axis of the rounding: {axis:?}");
    let d = c.draw(600.0, 200.0, None);
    // The value axis's labels, which are the ones anchored by their right edge to the axis.
    let ticks: Vec<_> = d
        .texts
        .iter()
        .filter(|t| t.anchor == chart::Anchor::Right)
        .collect();
    assert!(!ticks.is_empty(), "{:?}", d.texts);
    for t in ticks {
        assert!(
            t.text.len() <= 12,
            "a tick label of {} characters: {:?}",
            t.text.len(),
            t.text
        );
    }
}

/// **An axis of no length, or of a rounding's length, has one tick; an axis of no number has
/// none.** A run of one frame has a time axis from `t` to `t`, and scene 24's flat profile had a
/// value axis 1.7e-13 wide at 488 K: the first looped for ever looking for a power of ten under
/// zero, and the second counted ticks from 9.77e15, past where adding one changes the count, until
/// it asked for 34 GB.
#[test]
fn an_axis_too_short_to_divide_has_one_tick() {
    let axis = |lo: f64, hi: f64| chart::Axis {
        unit: "K".to_string(),
        lo,
        hi,
        lines: Vec::new(),
    };
    assert_eq!(axis(2.0, 2.0).ticks(false), [2.0]);
    // **First, a rounding-wide axis whose count still counts.** At 1.0 and 1e-14 wide the step is
    // 2e-15 and the count starts at 5e14, under 2^53, so without the guard this returns several
    // ticks rather than looping: it is the case that fails, and fails first, if the guard goes.
    assert_eq!(axis(1.0, 1.0 + 1e-14).ticks(false), [1.0]);
    // Then scene 24's, where without the guard the count starts past 2^53 and never moves.
    assert_eq!(
        axis(488.573_45, 488.573_45 + 1.7e-13).ticks(false),
        [488.573_45]
    );
    assert!(axis(f64::NAN, 1.0).ticks(false).is_empty());
    // And a span that is a span still divides: twelve ulps of 488 K is not, 1e-9 is.
    let fine = axis(488.0, 488.0 + 1e-9).ticks(false);
    assert!((2..=7).contains(&fine.len()), "{fine:?}");
    // The ends of the doubles, which `numerics-reviewer` found loop without allocating: a span
    // whose fifth rounds to nothing, and a span that overflows. One tick each, and a log axis
    // through zero or up to the largest double ends too. **Each on a clock**, because what these
    // guard against is a loop that never returns: without the guard the test fails in ten
    // seconds rather than never.
    let on_a_clock = move |lo: f64, hi: f64, log: bool| -> Vec<f64> {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(axis(lo, hi).ticks(log));
        });
        rx.recv_timeout(std::time::Duration::from_secs(10))
            .unwrap_or_else(|_| panic!("ticks of [{lo:e}, {hi:e}] did not return"))
    };
    assert_eq!(on_a_clock(0.0, 1e-323, false), [0.0]);
    assert_eq!(on_a_clock(-9.35e307, 9.35e307, false), [-9.35e307]);
    assert!(on_a_clock(0.0, 10.0, true).is_empty());
    // Up to the largest double, where `hi * (1 + 1e-9)` is infinite and so is every decade past it.
    let huge = on_a_clock(1e300, f64::MAX, true);
    assert_eq!(huge.len(), 9, "{huge:?}");

    // A run of one frame is charted, and drawn, with one tick on its time axis.
    let run = Run {
        title: "one frame".to_string(),
        frames: vec![Frame {
            t: 0.5,
            panels: Vec::new(),
            readings: vec![viewer_core::Reading {
                domain: "block".to_string(),
                label: "peak".to_string(),
                value: 300.0,
                unit: "K".to_string(),
            }],
        }],
    };
    let c = chart::chart(&run, false).expect("a chart");
    let d = c.draw(400.0, 200.0, Some(0.5));
    // The time axis's labels sit six pixels under the plot; its title, the legend and the unit
    // are elsewhere, and counting by anchor alone found those and could not fail.
    let bottom = c.plots(400.0, 200.0)[0][3];
    let times: Vec<&str> = d
        .texts
        .iter()
        .filter(|t| (t.at[1] - (bottom + 6.0)).abs() < 1e-9)
        .map(|t| t.text.as_str())
        .collect();
    // **Labelled `0.5`, not `0`.** With one tick the step fell back to one, and one asks for no
    // decimals.
    assert_eq!(times, ["0.5"], "{:?}", d.texts);
}

/// **Every tick says a different number, and zero says `0`.** The protein's softest mode is flat
/// at 2.1837e11 Hz, and on an axis stepping 2e8 two decimals of its exponent made four ticks of
/// `2.18e11`; its time axis, stepping 5e-12 s from zero, began `0.000000000000`.
#[test]
fn every_tick_says_a_different_number() {
    let line = Line {
        domain: "crambin".to_string(),
        label: "softest mode".to_string(),
        points: (0..17).map(|k| [k as f64 * 1.25e-12, 2.1837e11]).collect(),
    };
    let c = chart::graph(Across::Time, "Hz", vec![line]).expect("a chart");
    let d = c.draw(900.0, 240.0, None);
    let labels = |anchor: chart::Anchor| -> Vec<String> {
        d.texts
            .iter()
            .filter(|t| t.anchor == anchor)
            .map(|t| t.text.clone())
            .collect()
    };
    let values = labels(chart::Anchor::Right);
    assert!(values.len() >= 3, "{values:?}");
    let mut unique = values.clone();
    unique.dedup();
    assert_eq!(unique, values, "two ticks say the same number");
    let times = labels(chart::Anchor::Top);
    assert_eq!(times.first().map(String::as_str), Some("0"), "{times:?}");
    // **The last tick is labelled inside the axis.** Centred on the plot's right edge, half of
    // it hung past the drawing's margin and was cut by whatever held the drawing; it ends on its
    // tick now, so it is the one label anchored by its right edge below the plot.
    let ends: Vec<String> = d
        .texts
        .iter()
        .filter(|t| t.anchor == chart::Anchor::TopRight && t.text != "time / s")
        .filter(|t| !t.text.starts_with("crambin"))
        .map(|t| t.text.clone())
        .collect();
    assert_eq!(ends, ["2.0e-11"], "{:?}", d.texts);
    let mut unique = times.clone();
    unique.dedup();
    assert_eq!(unique, times, "two times say the same number");
}
