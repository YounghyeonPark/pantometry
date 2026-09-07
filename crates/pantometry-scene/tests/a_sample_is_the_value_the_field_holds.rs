//! **The run file did not carry the state, and nothing said so.**
//!
//! A six-cell block holding `0 0 100 0 0 0 °C` was written out as `0 0 90 0 0 0` — the peak 10%
//! low. The sampler placed its six points at `i/(n−1)` of the extent, which for a cell-based field
//! is the cell *boundaries*, and every consumer in the workspace reads that array: the viewer, the
//! report, glTF, USD and the CSV. `readings` said `peak 100` beside it, from the cells, so the
//! same frame carried two different answers.
//!
//! The fix is not a better place to sample. It is that a field now says where its own values are —
//! [`Lattice`] — and the sampler asks. Both conventions are real here: `Room` and `Hall` take
//! `dx = width/(nx−1)`, so their outermost values sit *on* the wall where a pressure antinode is,
//! and `Bar1D` and `Solid3D` hold cell averages with nothing on the boundary. A single convention
//! would have been wrong for one half or the other.
//!
//! What is checked is the property that matters: **asking a field for its own grid returns its own
//! numbers**, exactly, for either convention.

use pantometry_core::{Lattice, ScalarField};
use pantometry_scene::{sample_field, Extent, PanelData};
use pantometry_units::{LengthVec, Time};

/// A field of `n` known values along x, at whichever lattice it is told to be.
///
/// **Linear between its own positions**, which is what makes this test bite. The first version
/// took the nearest value instead, and a sabotage — the sampler going back to assuming corner to
/// corner for everything — *survived* it: nearest-neighbour over a monotone mapping lands on the
/// same six indices, so the wrong query points gave the right answer. `Solid3D` is trilinear
/// between cell centres and that is where the 10% came from, so the fixture interpolates too.
/// Measured: `0 0 100 0 0 0` at centred positions, queried at the nodal ones, gives 90 at index 2
/// — the exact number the run file used to carry.
struct Row {
    values: Vec<f64>,
    lattice: Lattice,
    length: f64,
}

impl ScalarField for Row {
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        let n = self.values.len();
        let x = p.to_si().x / self.length;
        // Linear between the two of this lattice's own positions that bracket `x`, clamped
        // outside them — the one-dimensional case of what `Solid3D` does between cell centres.
        let mut lo = 0;
        while lo + 1 < n && self.lattice.along(lo + 1, n) <= x {
            lo += 1;
        }
        if lo + 1 >= n {
            return self.values[n - 1];
        }
        let (a, b) = (self.lattice.along(lo, n), self.lattice.along(lo + 1, n));
        if x <= a {
            return self.values[lo];
        }
        let t = (x - a) / (b - a);
        self.values[lo] * (1.0 - t) + self.values[lo + 1] * t
    }

    fn lattice(&self) -> Lattice {
        self.lattice
    }
}

fn taken(values: &[f64], lattice: Lattice) -> Vec<f64> {
    let n = values.len();
    let length = 0.006;
    let field = Row {
        values: values.to_vec(),
        lattice,
        length,
    };
    let extent = Extent::new(
        LengthVec::ZERO,
        LengthVec::from_si(glam::DVec3::new(length, 0.0, 0.0)),
        n,
        1,
        1,
    );
    let panel = sample_field("row", &field, extent, Default::default(), Time::ZERO);
    match panel.data {
        PanelData::Field {
            values, lattice: l, ..
        } => {
            assert_eq!(l, lattice, "the panel did not record the field's lattice");
            values
        }
        _ => panic!("a field is a field"),
    }
}

/// **A cell-based field sampled on its own grid comes back exactly.**
///
/// The regression this file exists for. `0 0 100 0 0 0` used to arrive as `0 0 90 0 0 0`; a
/// blend of any kind fails this, and so does an off-by-one in the placement, because the spike is
/// one cell wide and asymmetric in the row.
#[test]
fn a_centred_field_arrives_intact() {
    let cells = [0.0, 0.0, 100.0, 0.0, 0.0, 0.0];
    assert_eq!(taken(&cells, Lattice::Centred), cells);
}

/// **And a nodal field does too**, which is what the old placement got right.
///
/// `Room`, `Hall` and `Well` are nodal, and their outermost values are *on* the boundary. Changing
/// the sampler globally to cell centres would have moved every one of them by half a cell and lost
/// the wall exactly where an antinode is, which is why this is a property of the field.
#[test]
fn a_nodal_field_arrives_intact() {
    let nodes = [7.0, 0.0, 0.0, 0.0, 0.0, 3.0];
    assert_eq!(taken(&nodes, Lattice::Nodal), nodes);
}

/// **The two conventions are different**, so the declaration is load-bearing rather than
/// decorative.
///
/// A test that passed for both would not notice a field declaring the wrong one. The same values
/// read through the wrong lattice put the spike in the wrong place: measured, `0 0 100 0 0 0` at
/// six cells has its centred position at 0.4167 of the box and its nodal one at 0.4, which is
/// 0.4 mm of a 6 mm row.
#[test]
fn asking_the_wrong_lattice_gives_a_different_answer() {
    let n = 6;
    let (a, b) = (Lattice::Centred.along(2, n), Lattice::Nodal.along(2, n));
    assert!(
        (a - 0.416_666_666_666_666_7).abs() < 1e-12,
        "centred sample 2 of 6 is at {a}"
    );
    assert!((b - 0.4).abs() < 1e-12, "nodal sample 2 of 6 is at {b}");
    // And the ends, which is where the difference is a boundary rather than an offset.
    assert_eq!(Lattice::Nodal.along(0, n), 0.0);
    assert_eq!(Lattice::Nodal.along(n - 1, n), 1.0);
    assert!(Lattice::Centred.along(0, n) > 0.0);
    assert!(Lattice::Centred.along(n - 1, n) < 1.0);
}

/// **One sample along an axis sits at its middle**, whichever the lattice.
///
/// A flat extent asked for one sample has nowhere else to put it, and an extent with real
/// thickness asked for one sample is best represented by its middle rather than by a corner. Both
/// conventions agree here and the code says so once rather than twice.
#[test]
fn a_single_sample_is_taken_at_the_middle() {
    assert_eq!(Lattice::Nodal.along(0, 1), 0.5);
    assert_eq!(Lattice::Centred.along(0, 1), 0.5);
}
