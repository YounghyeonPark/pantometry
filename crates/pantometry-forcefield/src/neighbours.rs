//! A Verlet neighbour list for [`PeriodicForceField`]'s real-space pairs: every pair within the
//! cutoff plus a skin, kept between evaluations and rebuilt only when an atom has moved far enough
//! to bring a pair from outside it to inside the cutoff. **Opt-in**
//! ([`PeriodicForceField::with_neighbour_list`]), and **the same bits** as the cell list it
//! replaces: the same pairs, summed in the same order.
//!
//! # What is kept
//!
//! At a build, every pair `i < j` that is not excluded and whose minimum-image separation is
//! shorter than `R = (r_c + s)(1 + 10⁻⁹)`, found by the cell list's own search at `R`; the
//! positions as they were given; and the box, the cutoff and the excluded pairs they were built
//! for. Each atom's partners are a run of one array. An evaluation then visits only those pairs,
//! which for `s` = 2 Å at `r_c` = 9 Å is `(11/9)³` = 1.83 times the cutoff sphere, against the
//! 125 half-cutoff cells' 3.7.
//!
//! # When it is rebuilt, and why that is enough
//!
//! **Rebuilt when any atom's displacement since the build, by minimum image, is longer than
//! `s/2`**, or when the atom count, the box, the cutoff or the excluded pairs differ from the
//! build's. The last is what makes **no excluded pair ever reach the evaluation**, so that
//! [`PeriodicForceField`]'s pair term does not look each pair up again: that lookup, a binary
//! search of the atom's exclusions, measured 2.5 ms of the 38.7 ms the forces of a 5 089-atom
//! octanol box took with it, and comparing the exclusions is one pass over them. The proof of the
//! first:
//!
//! - Let `u_i = mi(r_i(t) − r_i(0))`, the minimum image of atom `i`'s displacement. Then
//!   `r_i(t) = r_i(0) + u_i + n L` for some whole `n`, whatever path the atom took and however many
//!   times it crossed a face — the minimum image of the displacement is all a periodic system can
//!   see of it.
//! - The minimum-image distance `ρ(x) = min_n |x + n L|` is 1-Lipschitz, being a minimum of
//!   1-Lipschitz functions, and a whole number of box lengths does not change it. So the pair's
//!   distance at `t` is `ρ(x_ij(0) + u_i − u_j) ≥ ρ(x_ij(0)) − |u_i| − |u_j|`.
//! - While no atom has `|u| > s/2`, a pair left out of the list, `ρ(x_ij(0)) ≥ R ≥ r_c + s`, is
//!   still at `ρ ≥ r_c` — and the cutoff takes only `ρ < r_c`. So **every pair inside the cutoff
//!   is in the list**, as long as the rule holds.
//! - The `10⁻⁹` margin on `R` is what makes the floating-point version of that argument true: a
//!   distance and a displacement are each computed to a few ulps of the coordinates, about
//!   `10⁻²⁵` m for positions within a few box lengths of the origin, and the margin is 10⁻¹⁸ m.
//!   It adds pairs; it removes none.
//!
//! Both halves of the rule are needed, and the bound is tight: two atoms moving straight at each
//! other each by `s/2` close a pair at `r_c + s` to `r_c`. A rule at `s` per atom would let that
//! pair reach `r_c − s`. **A frozen atom never moves**: [`crate::MolecularDynamics`] never writes
//! it, so its displacement is exactly zero and a box with its protein frozen rebuilds on the
//! water's motion alone. An atom given in another image — moved by a whole box length, as the
//! force field allows — has not moved either, which is why the displacement is a minimum image.
//!
//! # The order: the cell list's, reproduced
//!
//! The cell list ([`crate::periodic`]'s `for_each_pair`) sums by `i`, then by the 5³ cells about
//! `i`'s in a fixed order, then by `j` within a cell; which cell an atom is in is fixed by its
//! wrapped position at that evaluation, not at the build. So **each evaluation re-derives that
//! order**: it puts every atom in the cutoff's cells as the cell list would, gives each listed
//! partner `j` of `i` the key `(the index of j's cell among the 125 about i's, j)`, and visits
//! `i`'s partners in key order. A partner whose cell is outside the 125 is not visited, as the cell
//! list would not visit it, and is at least the cutoff away. A box fewer than five cells across
//! visits every pair by `j`, and so does the list.
//!
//! The pairs the list visits inside the cutoff are then the cell list's exactly — those the cell
//! list visits and that are in the list, which while the rule holds is all of them — in the cell
//! list's order, with the separation computed by the same function from the same wrapped positions.
//! Every sum that reads them is the same sequence of additions, and gives the same bits. Changing
//! both to one order (by `i`, then `j`) would be simpler and was not done: it would move every bit
//! the default path has, and tests pin them.
//!
//! Re-deriving the order is cheap because it barely changes: the keys are recomputed at every
//! evaluation, and each atom's run is put in order by insertion, which on a run that was in order
//! at the last evaluation costs a comparison a partner; after a build, by a full sort. The keys are
//! distinct, so any sort gives the same order.
//!
//! # What is shared
//!
//! The list lives inside the force field, behind a [`RefCell`] — a cache, not a parameter — so that
//! [`crate::Potential::energy_and_forces`] can keep it through `&self`. Two force fields that differ
//! only in their lists' state compare equal. A [`crate::PeriodicDecoupling`] built from a field with
//! a list keeps a list of its own; its evaluation's pairs go through it, and its
//! [`crate::Alchemical::couplings`], which reads only the group's pairs with the rest by brute
//! force, does not use one. [`PeriodicForceField::with_cell`] starts its list empty.
//!
//! # The cost, measured
//!
//! `tests/a_neighbour_list_is_the_cell_lists_bits.rs`, `the_cost_on_w4s_box_measured`: W4's
//! box, 30 105 atoms in 62.35 × 65.47 × 74.82 Å, the mesh at δ = 10⁻⁶, `r_c` = 9 Å, release, one
//! core, `x86_64-pc-windows-gnu`, the first of two runs. Melted 0.15 ps at 0.5 fs with the protein
//! frozen, then from that one state 50 steps of 2 fs in a 1 ps⁻¹ bath each, every trajectory's
//! end, positions and velocities, the cell list's to the bit:
//!
//! | | pairs kept | ms an evaluation | ms one that builds | ms a step | builds in 50 steps |
//! | --- | --- | --- | --- | --- | --- |
//! | cell list | — | 400.3 | — | 400.5 | — |
//! | skin 1 Å | 6 165 129 | 296.3 | 528.7 | 355.5 | 12 |
//! | skin 1.5 Å | 7 192 760 | 296.3 | 568.8 | 351.6 | 8 |
//! | skin 2 Å | 8 279 403 | 295.7 | 622.0 | 337.7 | 5 |
//! | skin 2.5 Å | 9 430 201 | 315.4 | 666.2 | 347.1 | 4 |
//!
//! The second run, on the same machine, read every row about a fifth faster, the cell list's
//! 329.9 ms a step against 285.9 at 2 Å and 285.5 at 1.5 Å, with the same builds.
//!
//! **An evaluation without a build is 104 ms cheaper, 26%; a step is 63 ms cheaper, 16%, at
//! 2 Å** ([`PeriodicForceField::NEIGHBOUR_SKIN`]), 13% in the second run, because a build costs
//! 190–350 ms — the cell list's search at `r_c + s`, which visits more pairs than the evaluation
//! it saves — and comes every ten steps at 2 Å, every four at 1 Å. Fifty steps are few: the
//! builds are counts of 4–12, and the steps' times two runs each: 1.5 and 2 Å are not told
//! apart (2 Å ahead by 14 ms in one, behind by 0.4 in the other), and 1 and 2.5 Å are slower in
//! both. The pair loop's own iteration, its terms aside, measured once on a lattice of 10 648
//! waters at rest (31 944 atoms, a 1.5 Å skin, a probe not kept): 121 ms by cells, 45 ms by the
//! list, and 271 ms for a build.
//!
//! # Not here
//!
//! **No sharing across fields**: a [`RefCell`] is not `Sync`, and nothing in this crate shares a
//! force field between threads. An evaluation's spans of atoms are visited on several threads,
//! and a build's searched on several, but only the calling thread borrows the list from its cell,
//! and it adds the spans' pairs in order (see [`crate::periodic`]'s module documentation). **No
//! rebuild by the two largest displacements** (`|u|₁ + |u|₂ > s`), which is the tight form of the
//! same bound and would rebuild less often; the rule here is the one whose proof needs nothing
//! about which atoms moved. **No list of the
//! excluded or decoupled pairs**: those are the evaluation's own, by atom.

use crate::ewald::sub;
use crate::ordered;
use crate::periodic::{
    cell_grid, cell_of, separation, CellSearch, Exclusions, PeriodicBox, PeriodicForceField, REACH,
};
use std::cell::RefCell;
use std::fmt;

/// The relative margin on the list's radius: see the module documentation.
pub(crate) const MARGIN: f64 = 1e-9;

/// A key past every cell's: the partner is outside the 125 cells about `i`'s.
const OUTSIDE: u64 = 127;

/// The low half of an entry: the partner's index.
const ATOM: u64 = 0xFFFF_FFFF;

/// What a force field's neighbour list has done ([`PeriodicForceField::neighbour_list`]).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NeighbourListStatus {
    /// The skin, metres.
    pub skin: f64,
    /// How many times the list has been built.
    pub builds: u64,
    /// How many evaluations have gone through it, the builds' included.
    pub evaluations: u64,
    /// How many pairs the last build kept: every pair within the cutoff and the skin, exclusions
    /// aside.
    pub pairs: usize,
}

/// A neighbour list with skin `skin`: see the module documentation.
#[derive(Clone)]
pub(crate) struct NeighbourList {
    skin: f64,
    state: RefCell<State>,
}

#[derive(Clone, Default)]
struct State {
    /// The box edges and cutoff of the build.
    lengths: [f64; 3],
    cutoff: f64,
    /// The pairs the build left out, which every evaluation's must be: what lets the pair term
    /// skip its own check of each pair.
    excluded: Vec<[usize; 2]>,
    /// The positions as given at the build; empty before the first.
    reference: Vec<[f64; 3]>,
    /// Atom `i`'s partners are `entries[start[i]..start[i + 1]]`, each `key << 32 | j`.
    start: Vec<usize>,
    entries: Vec<u64>,
    /// Whether the entries have not been put in the cell list's order since the build.
    unsorted: bool,
    /// Each atom's cell at the last evaluation, when the box is cut into cells.
    home: Vec<[usize; 3]>,
    builds: u64,
    evaluations: u64,
    /// How many threads the last build spawned, which the tests read.
    spawned: usize,
}

impl PartialEq for NeighbourList {
    /// Equal when the skins are: the rest is a cache, and while it is valid it changes no bit.
    fn eq(&self, other: &NeighbourList) -> bool {
        self.skin == other.skin
    }
}

impl fmt::Debug for NeighbourList {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = self.status();
        f.debug_struct("NeighbourList")
            .field("skin", &s.skin)
            .field("builds", &s.builds)
            .field("evaluations", &s.evaluations)
            .field("pairs", &s.pairs)
            .finish()
    }
}

impl NeighbourList {
    /// An empty list with skin `skin`, metres.
    ///
    /// # Panics
    ///
    /// If the skin is not positive and finite.
    pub(crate) fn new(skin: f64) -> NeighbourList {
        assert!(
            skin.is_finite() && skin > 0.0,
            "a neighbour list's skin must be positive and finite: {skin}"
        );
        NeighbourList {
            skin,
            state: RefCell::new(State::default()),
        }
    }

    /// The same skin and nothing built.
    pub(crate) fn emptied(&self) -> NeighbourList {
        NeighbourList::new(self.skin)
    }

    /// How many threads the last build spawned.
    #[cfg(test)]
    pub(crate) fn spawned_at_build(&self) -> usize {
        self.state.borrow().spawned
    }

    pub(crate) fn status(&self) -> NeighbourListStatus {
        let s = self.state.borrow();
        NeighbourListStatus {
            skin: self.skin,
            builds: s.builds,
            evaluations: s.evaluations,
            pairs: s.entries.len(),
        }
    }

    /// Whether the list must be built again before an evaluation at `at` leaving out `ex`: see the
    /// module documentation.
    fn stale(
        &self,
        s: &State,
        cell: &PeriodicBox,
        at: &[[f64; 3]],
        cutoff: f64,
        ex: &Exclusions,
    ) -> bool {
        if s.reference.len() != at.len()
            || s.lengths != cell.lengths()
            || s.cutoff != cutoff
            || s.excluded != ex.pairs()
        {
            return true;
        }
        let limit = 0.25 * self.skin * self.skin;
        at.iter().zip(&s.reference).any(|(p, q)| {
            let u = cell.minimum_image(sub(*p, *q));
            u[0] * u[0] + u[1] * u[1] + u[2] * u[2] > limit
        })
    }

    /// Builds the list by the cell list's search at the larger radius, each span of atoms' partners
    /// found on any of `threads` threads and kept in the spans' order, so the entries are the same
    /// on any number.
    #[allow(clippy::too_many_arguments)]
    fn build(
        &self,
        s: &mut State,
        cell: &PeriodicBox,
        at: &[[f64; 3]],
        cutoff: f64,
        ex: &Exclusions,
        threads: usize,
        pieces: usize,
    ) {
        let n = at.len();
        assert!(
            u64::try_from(n).is_ok_and(|n| n <= ATOM),
            "a neighbour list holds at most 2³² atoms"
        );
        let radius = (cutoff + self.skin) * (1.0 + MARGIN);
        let search = CellSearch::new(cell, at, radius);
        s.entries.clear();
        s.start.clear();
        s.start.push(0);
        let (entries, start) = (&mut s.entries, &mut s.start);
        // Each atom's partners are one run, and a span's are its atoms' runs in order:
        // `ends[k]` is where atom `k` of the span ends among them.
        let spawned = ordered::in_order(
            threads,
            ordered::ranges(n, pieces),
            &mut Vec::new(),
            |atoms, out: &mut Found| {
                for i in atoms {
                    search.atom(i, &mut |i, j, _, _| {
                        if !ex.contains(i, j) {
                            out.partners.push(j as u64);
                        }
                    });
                    out.ends.push(out.partners.len());
                }
            },
            |out| {
                let base = entries.len();
                entries.extend_from_slice(&out.partners);
                start.extend(out.ends.iter().map(|e| base + e));
                out.partners.clear();
                out.ends.clear();
            },
        );
        s.reference.clear();
        s.reference.extend_from_slice(at);
        s.lengths = cell.lengths();
        s.cutoff = cutoff;
        s.excluded.clear();
        s.excluded.extend_from_slice(ex.pairs());
        s.unsorted = true;
        s.builds += 1;
        s.spawned = spawned;
    }

    /// Calls `f(i, j, d, r²)` for exactly the pairs, and in exactly the order, that the cell list
    /// at `cutoff` calls it for, the pairs `exclusions` holds aside — building the list first if it
    /// is stale, which it is whenever `exclusions` differs from its build's. So **no pair of
    /// `exclusions` reaches `f`**, and the caller need not check. The evaluation's own way through
    /// it is [`NeighbourList::for_each_span`]; this is that on one thread, which the tests hold to
    /// the cell list.
    #[cfg(test)]
    pub(crate) fn for_each_pair(
        &self,
        cell: &PeriodicBox,
        at: &[[f64; 3]],
        cutoff: f64,
        exclusions: &Exclusions,
        mut f: impl FnMut(usize, usize, [f64; 3], f64),
    ) {
        self.for_each_span(cell, at, cutoff, exclusions, 1, 1, |spans| {
            for span in spans {
                span.visit(&mut f);
            }
        });
    }

    /// The pairs of `for_each_pair` in spans of atoms: prepares the evaluation — building the list
    /// on `threads` threads if it is stale, and placing every atom in the cutoff's cells — and
    /// hands `visit` the atoms cut into at most `pieces` spans of about equal numbers of entries, in
    /// order. Visiting every span once, in order, is `for_each_pair`; each span
    /// writes only its own atoms' runs and reads only what every span reads, so the spans may be
    /// visited on any threads, and each gives the same pairs and separations.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn for_each_span(
        &self,
        cell: &PeriodicBox,
        at: &[[f64; 3]],
        cutoff: f64,
        exclusions: &Exclusions,
        threads: usize,
        pieces: usize,
        visit: impl FnOnce(Vec<Span<'_>>),
    ) {
        let mut guard = self.state.borrow_mut();
        let s = &mut *guard;
        s.evaluations += 1;
        if self.stale(s, cell, at, cutoff, exclusions) {
            self.build(s, cell, at, cutoff, exclusions, threads, pieces);
        }
        let n = at.len();
        let l = cell.lengths();
        let fresh = s.unsorted;
        let w: Vec<[f64; 3]> = at.iter().map(|p| cell.wrap(*p)).collect();
        let order = cell_grid(l, cutoff).map(|cells| {
            // Per axis, the index among the five cells about `i`'s of the cell `k − cells` from
            // it, or past them: the cell list's `(c_i + cells + dx − reach) % cells = c_j` solved
            // for `dx`.
            let offsets = cells.map(|c| {
                (0..2 * c)
                    .map(|k| {
                        let dx = (k + REACH) % c;
                        if dx <= 2 * REACH {
                            dx as u64
                        } else {
                            OUTSIDE
                        }
                    })
                    .collect::<Vec<u64>>()
            });
            let home: Vec<[usize; 3]> = w.iter().map(|p| cell_of(*p, l, cells)).collect();
            // An atom whose cell is the last evaluation's keeps its keys with every partner whose
            // cell is too.
            let moved: Vec<bool> = if fresh || s.home.len() != n {
                vec![true; n]
            } else {
                home.iter().zip(&s.home).map(|(a, b)| a != b).collect()
            };
            Order {
                cells,
                offsets,
                home,
                moved,
            }
        });
        let view = View {
            w,
            l,
            half: l.map(|x| 0.5 * x),
            rc2: cutoff * cutoff,
            fresh,
            order,
        };
        visit(spans(&s.start, &mut s.entries, pieces, &view));
        s.home = view.order.map_or_else(Vec::new, |o| o.home);
        s.unsorted = false;
    }
}

/// What every span of an evaluation reads: the wrapped positions and, when the box is cut into
/// cells, each atom's cell and the keys' tables.
pub(crate) struct View {
    w: Vec<[f64; 3]>,
    l: [f64; 3],
    half: [f64; 3],
    rc2: f64,
    /// Whether the runs have not been put in order since the build.
    fresh: bool,
    order: Option<Order>,
}

/// The cell list's order, re-derived at an evaluation: see the module documentation.
struct Order {
    cells: [usize; 3],
    offsets: [Vec<u64>; 3],
    home: Vec<[usize; 3]>,
    /// Whether each atom's cell is not the last evaluation's.
    moved: Vec<bool>,
}

impl Order {
    /// The key of a partner in cell `cj` of an atom in cell `ci`: the index of `cj` among the 125
    /// cells about `ci`, in the cell list's order, or [`OUTSIDE`].
    #[inline(always)]
    fn key(&self, ci: [usize; 3], cj: [usize; 3]) -> u64 {
        let (cells, span) = (self.cells, (2 * REACH + 1) as u64);
        let ox = self.offsets[0][cj[0] + cells[0] - ci[0]];
        let oy = self.offsets[1][cj[1] + cells[1] - ci[1]];
        let oz = self.offsets[2][cj[2] + cells[2] - ci[2]];
        if ox == OUTSIDE || oy == OUTSIDE || oz == OUTSIDE {
            OUTSIDE
        } else {
            (ox * span + oy) * span + oz
        }
    }
}

/// A run of atoms and their entries, which a span alone writes ([`NeighbourList::for_each_span`]).
pub(crate) struct Span<'a> {
    atoms: std::ops::Range<usize>,
    /// The whole list's starts, so `entries` begins at `start[atoms.start]`.
    start: &'a [usize],
    entries: &'a mut [u64],
    view: &'a View,
}

/// The atoms cut into at most `pieces` spans of about equal numbers of entries, in order.
fn spans<'a>(
    start: &'a [usize],
    entries: &'a mut [u64],
    pieces: usize,
    view: &'a View,
) -> Vec<Span<'a>> {
    let n = start.len() - 1;
    let per = start[n].div_ceil(pieces.max(1)).max(1);
    let mut out = Vec::new();
    let mut rest = entries;
    let mut a = 0;
    while a < n {
        let mut b = a + 1;
        while b < n && start[b] - start[a] < per {
            b += 1;
        }
        // The atoms after `b` have no entries — the last atom never has a partner `j > i` — so
        // this span takes them rather than leaving them one of their own past `pieces`.
        if start[b] == start[n] {
            b = n;
        }
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(start[b] - start[a]);
        out.push(Span {
            atoms: a..b,
            start,
            entries: head,
            view,
        });
        rest = tail;
        a = b;
    }
    out
}

impl Span<'_> {
    /// Puts each of the span's atoms' runs in the cell list's order and calls `f(i, j, d, r²)` for
    /// its pairs inside the cutoff, in that order.
    pub(crate) fn visit(self, f: &mut impl FnMut(usize, usize, [f64; 3], f64)) {
        let Span {
            atoms,
            start,
            entries,
            view: v,
        } = self;
        let base = start[atoms.start];
        let mut scratch = Vec::new();
        for i in atoms {
            let run = &mut entries[start[i] - base..start[i + 1] - base];
            if let Some(o) = &v.order {
                if run.is_empty() {
                    continue;
                }
                let ci = o.home[i];
                let mut changed = false;
                for e in run.iter_mut() {
                    let j = (*e & ATOM) as usize;
                    if o.moved[i] || o.moved[j] {
                        *e = o.key(ci, o.home[j]) << 32 | j as u64;
                        changed = true;
                    }
                }
                if v.fresh {
                    by_key(run, &mut scratch);
                }
                if changed {
                    insertion_sort(run);
                }
                for &e in run.iter() {
                    if e >> 32 == OUTSIDE {
                        break;
                    }
                    let j = (e & ATOM) as usize;
                    let (d, r2) = separation(v.w[i], v.w[j], v.l, v.half);
                    if r2 < v.rc2 {
                        f(i, j, d, r2);
                    }
                }
            } else {
                // Every pair, by `j`: the build's search at the larger radius is every pair too,
                // by `j`, so the runs are in order, and the sort only says so.
                if v.fresh {
                    insertion_sort(run);
                }
                for &e in run.iter() {
                    let j = (e & ATOM) as usize;
                    let (d, r2) = separation(v.w[i], v.w[j], v.l, v.half);
                    if r2 < v.rc2 {
                        f(i, j, d, r2);
                    }
                }
            }
        }
    }
}

/// A span of atoms' partners found by a build, and where each atom's run ends among them.
#[derive(Default)]
struct Found {
    partners: Vec<u64>,
    ends: Vec<usize>,
}

/// Puts a run in order of its keys alone by counting, through `scratch`: what is left for
/// [`insertion_sort`] is the order by `j` within a key, a few entries each.
fn by_key(run: &mut [u64], scratch: &mut Vec<u64>) {
    let mut count = [0usize; OUTSIDE as usize + 2];
    for &e in run.iter() {
        count[(e >> 32) as usize + 1] += 1;
    }
    for k in 0..=OUTSIDE as usize {
        count[k + 1] += count[k];
    }
    scratch.clear();
    scratch.resize(run.len(), 0);
    for &e in run.iter() {
        let k = (e >> 32) as usize;
        scratch[count[k]] = e;
        count[k] += 1;
    }
    run.copy_from_slice(scratch);
}

/// Sorts a run that is nearly in order already: a comparison an entry when it is in order.
fn insertion_sort(v: &mut [u64]) {
    for k in 1..v.len() {
        let x = v[k];
        if v[k - 1] <= x {
            continue;
        }
        let mut m = k;
        while m > 0 && v[m - 1] > x {
            v[m] = v[m - 1];
            m -= 1;
        }
        v[m] = x;
    }
}

impl PeriodicForceField {
    /// The skin [`PeriodicForceField::with_neighbour_list`] is measured best at on W4's box, 2 Å,
    /// in metres: 338 and 286 ms a step at 2 fs in two runs against the cell list's 401 and 330, a
    /// build every ten steps, and 1.5 Å not told apart from it. See the module documentation of
    /// [`crate::neighbours`].
    pub const NEIGHBOUR_SKIN: f64 = 2e-10;

    /// The same force field with its real-space pairs kept in a Verlet neighbour list of skin
    /// `skin` metres ([`crate::neighbours`]): **the same bits** as without it, energy, forces and
    /// virial, at every evaluation; it is rebuilt whenever an atom has moved more than `skin/2`
    /// since the last build. [`PeriodicForceField::NEIGHBOUR_SKIN`] is the skin measured best on W4's
    /// box. The list may reach past half the box: the pairs inside the cutoff are still the cell
    /// list's.
    ///
    /// # Panics
    ///
    /// If the skin is not positive and finite.
    pub fn with_neighbour_list(mut self, skin: f64) -> PeriodicForceField {
        self.neighbours = Some(NeighbourList::new(skin));
        self
    }

    /// The same force field searching its pairs by the cell list at every evaluation, the default.
    pub fn without_neighbour_list(mut self) -> PeriodicForceField {
        self.neighbours = None;
        self
    }

    /// What the neighbour list has done, or `None` without one.
    pub fn neighbour_list(&self) -> Option<NeighbourListStatus> {
        self.neighbours.as_ref().map(NeighbourList::status)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ccd::ANGSTROM;
    use crate::dynamics::{Bath, Potential};
    use crate::ewald::EwaldParameters;
    use crate::periodic::for_each_pair;
    use crate::water::WaterBox;

    type Call = (usize, usize, [u64; 3], u64);

    /// Every call of the cell list at `at`, excluded pairs aside: `(i, j, d, r²)` as bits.
    fn by_cells(field: &PeriodicForceField, at: &[[f64; 3]]) -> Vec<Call> {
        let mut out = Vec::new();
        let cutoff = field.ewald().parameters().cutoff;
        for_each_pair(&field.cell(), at, cutoff, |i, j, d, r2| {
            if !field.exclusions.contains(i, j) {
                out.push((i, j, d.map(f64::to_bits), r2.to_bits()));
            }
        });
        out
    }

    /// Every call of the list at `at`, the same way.
    fn by_list(field: &PeriodicForceField, list: &NeighbourList, at: &[[f64; 3]]) -> Vec<Call> {
        let mut out = Vec::new();
        list.for_each_pair(
            &field.cell(),
            at,
            field.ewald().parameters().cutoff,
            &field.exclusions,
            |i, j, d, r2| out.push((i, j, d.map(f64::to_bits), r2.to_bits())),
        );
        out
    }

    /// Every pair `i < j`, not excluded, nearer than the cutoff in one of the 27 images about the
    /// rounded minimum image: no cell, no wrapping, no list.
    fn brute(field: &PeriodicForceField, at: &[[f64; 3]]) -> Vec<(usize, usize)> {
        let cell = field.cell();
        let l = cell.lengths();
        let rc = field.ewald().parameters().cutoff;
        let mut out = Vec::new();
        for i in 0..at.len() {
            for j in i + 1..at.len() {
                if field.exclusions.contains(i, j) {
                    continue;
                }
                let d0 = cell.minimum_image(sub(at[i], at[j]));
                let mut nearest = f64::INFINITY;
                for n in 0..27 {
                    let s = [n / 9, (n / 3) % 3, n % 3].map(|k| k as f64 - 1.0);
                    let d = [0, 1, 2].map(|a| d0[a] + s[a] * l[a]);
                    nearest = nearest.min(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                }
                if nearest < rc * rc {
                    out.push((i, j));
                }
            }
        }
        out
    }

    /// A box of rigid TIP3P moving under SETTLE in a 300 K bath, 1 fs a step: the positions at
    /// the start and after each step.
    fn trajectory(
        per_side: usize,
        cutoff: f64,
        steps: usize,
    ) -> (PeriodicForceField, Vec<Vec<[f64; 3]>>) {
        let water = WaterBox::lattice(per_side, 0x7E57);
        let field = water.force_field(cutoff, 1e-5);
        let mut at = water.positions().to_vec();
        let mut md = water
            .dynamics()
            .with_bath(Bath::Langevin {
                temperature: 300.0,
                friction: 5e12,
                seed: 3,
            })
            .thermalised(&at, 300.0, 3);
        let mut out = vec![at.clone()];
        for _ in 0..steps {
            md.step(&field, &mut at, 1e-15);
            out.push(at.clone());
        }
        (field, out)
    }

    /// **The list's calls are the cell list's, to the bit and in order, at every step of a
    /// trajectory that crosses several builds, and their pairs are the brute-force set over 27
    /// images every tenth step.** Three boxes: one six half-cutoff cells across, where the 125
    /// cells about an atom are not the whole box; one fewer than five, where every pair is visited;
    /// and one whose list reaches past half the box.
    #[test]
    fn the_list_is_the_cell_lists_calls_along_a_trajectory() {
        // Waters a side, cutoff and skin (Å), steps, the least number of builds.
        for (per_side, cutoff, skin, steps, least) in [
            (5, 4.6, 1.0, 200, 3),
            (3, 4.6, 1.0, 200, 3),
            (3, 4.6, 2.5, 100, 1),
        ] {
            let (field, path) = trajectory(per_side, cutoff * ANGSTROM, steps);
            let lengths = field.cell().lengths();
            let grid = cell_grid(lengths, cutoff * ANGSTROM);
            assert_eq!(grid.is_some(), per_side == 5, "{per_side}: {grid:?}");
            if let Some(cells) = grid {
                assert!(cells.iter().all(|&c| c > 5), "{cells:?}");
            }
            if skin > 2.0 {
                assert!((cutoff + skin) * ANGSTROM > 0.5 * lengths[0]);
            }
            let list = NeighbourList::new(skin * ANGSTROM);
            let mut inside = 0;
            for (k, at) in path.iter().enumerate() {
                let cells = by_cells(&field, at);
                let listed = by_list(&field, &list, at);
                assert_eq!(listed, cells, "{per_side}, {skin} Å, step {k}");
                if k % 10 == 0 {
                    let mut set: Vec<(usize, usize)> = listed.iter().map(|c| (c.0, c.1)).collect();
                    set.sort_unstable();
                    assert_eq!(set, brute(&field, at), "{per_side}, {skin} Å, step {k}");
                }
                inside += listed.len();
            }
            let s = list.status();
            assert_eq!(s.evaluations, path.len() as u64);
            println!(
                "{per_side}, {skin} Å: {} builds in {}",
                s.builds, s.evaluations
            );
            assert!(
                s.builds >= least && s.builds < s.evaluations / 4,
                "{per_side}, {skin} Å: {} builds",
                s.builds
            );
            assert!(inside > 100 * path.len(), "{inside} pairs inside");
        }
    }

    /// Two waters in a 30 Å cube, `r_c` = 5 Å (twelve cells an axis), their oxygens `r` apart
    /// along x, each water's hydrogens pointing away from the other.
    fn two_waters(r: f64) -> (PeriodicForceField, Vec<[f64; 3]>) {
        let cell = crate::PeriodicBox::cubic(30.0 * ANGSTROM);
        let p = EwaldParameters::for_accuracy(&cell, 5.0 * ANGSTROM, 1e-5);
        let field = PeriodicForceField::tip3p(cell, p, 2);
        let (oh, half) = (crate::water::OH_LENGTH, 0.5 * crate::water::hoh_angle());
        let water = |o: [f64; 3], away: f64| {
            [
                o,
                [o[0] + away * oh * half.cos(), o[1] + oh * half.sin(), o[2]],
                [o[0] + away * oh * half.cos(), o[1] - oh * half.sin(), o[2]],
            ]
        };
        let (a, y) = (10.0 * ANGSTROM, 15.0 * ANGSTROM);
        let mut at = water([a, y, y], -1.0).to_vec();
        at.extend(water([a + r, y, y], 1.0));
        (field, at)
    }

    /// The first water moved by `a` and the second by `b`, along x.
    fn moved(at: &[[f64; 3]], a: f64, b: f64) -> Vec<[f64; 3]> {
        at.iter()
            .enumerate()
            .map(|(i, p)| [p[0] + if i < 3 { a } else { b }, p[1], p[2]])
            .collect()
    }

    /// The listed field's evaluation and the plain one's are the same bits, through `evaluate`
    /// and through the dynamics' `energy_and_forces`, and so are their calls: three evaluations of
    /// the list.
    fn same_bits(plain: &PeriodicForceField, listed: &PeriodicForceField, at: &[[f64; 3]]) {
        assert_eq!(
            evaluation_bits(&plain.evaluate(at)),
            evaluation_bits(&listed.evaluate(at))
        );
        let (mut f1, mut f2) = (vec![[0.0; 3]; at.len()], vec![[0.0; 3]; at.len()]);
        let e1 = plain.energy_and_forces(at, &mut f1);
        let e2 = listed.energy_and_forces(at, &mut f2);
        assert_eq!(e1.to_bits(), e2.to_bits());
        assert_eq!(vector_bits(&f1), vector_bits(&f2));
        // And the calls themselves, `d` to the bit: a separation's sign of zero is absorbed by the
        // sums above (`0.0 + −0.0 = 0.0`), and these waters, all in one plane, have such zeros.
        if let Some(list) = listed.neighbours.as_ref() {
            assert_eq!(by_list(listed, list, at), by_cells(plain, at), "the calls");
        }
    }

    fn vector_bits(v: &[[f64; 3]]) -> Vec<[u64; 3]> {
        v.iter().map(|p| p.map(f64::to_bits)).collect()
    }

    /// An evaluation as bits — every energy term, every force, the virial — so that `−0.0` is not
    /// `0.0`, as `==` would have it.
    fn evaluation_bits(e: &crate::PeriodicEvaluation) -> (Vec<u64>, Vec<[u64; 3]>, [[u64; 3]; 3]) {
        let (t, w) = (&e.energy, &e.energy.ewald);
        let terms = [
            t.bond,
            t.angle,
            t.torsion,
            t.inversion,
            t.van_der_waals,
            t.dispersion_correction,
            t.electrostatic,
            t.total,
            w.real,
            w.reciprocal,
            w.self_energy,
            w.excluded,
            w.background,
            w.total,
        ];
        (
            terms.map(f64::to_bits).to_vec(),
            vector_bits(&e.forces),
            e.virial.map(|r| r.map(f64::to_bits)),
        )
    }

    /// The comparison sees a sign of zero, which `==` does not.
    #[test]
    fn the_bits_comparison_sees_the_sign_of_zero() {
        let (plain, at) = two_waters(6.0 * ANGSTROM);
        let e = plain.evaluate(&at);
        let mut flipped = e.clone();
        flipped.virial[0][0] = -0.0;
        flipped.forces[0][0] = -0.0;
        let mut zero = e.clone();
        zero.virial[0][0] = 0.0;
        zero.forces[0][0] = 0.0;
        assert_eq!(flipped, zero, "== calls them equal");
        assert_ne!(evaluation_bits(&flipped), evaluation_bits(&zero));
    }

    fn list_of(field: &PeriodicForceField) -> &NeighbourList {
        field.neighbours.as_ref().expect("a list")
    }

    fn calls_pair(field: &PeriodicForceField, at: &[[f64; 3]], i: usize, j: usize) -> bool {
        by_list(field, list_of(field), at)
            .iter()
            .any(|c| (c.0, c.1) == (i, j))
    }

    /// **A pair that starts inside the skin and ends inside the cutoff is not missed**: the
    /// oxygens at `r_c + s − ε`, inside the list, each water moved `s/2 − ε/4` towards the other,
    /// so that no atom has moved `s/2` and the list is not rebuilt, and the oxygens end at
    /// `r_c − ε/2`. The pair is in the evaluation, and the list's evaluation is the cell list's.
    #[test]
    fn a_pair_from_the_skin_to_the_cutoff_is_not_missed() {
        let (rc, s, eps) = (5.0 * ANGSTROM, 1.0 * ANGSTROM, 1e-3 * ANGSTROM);
        let (plain, at) = two_waters(rc + s - eps);
        let listed = plain.clone().with_neighbour_list(s);
        same_bits(&plain, &listed, &at);
        assert!(!calls_pair(&listed, &at, 0, 3));
        let step = 0.5 * s - 0.25 * eps;
        let end = moved(&at, step, -step);
        let r = end[3][0] - end[0][0];
        assert!(r < rc && r > rc - eps, "{r:e}");
        same_bits(&plain, &listed, &end);
        assert!(calls_pair(&listed, &end, 0, 3));
        assert_eq!(listed.neighbour_list().unwrap().builds, 1);
        // The pair is there to be missed: its van der Waals is a part of the energy that the cell
        // list's own evaluation would lose with it.
        let vdw = plain.energy(&end).van_der_waals;
        assert!(
            vdw < 0.0 && plain.energy(&at).van_der_waals == 0.0,
            "{vdw:e}"
        );
    }

    /// **Two atoms that each move less than the skin but together more are not missed**: the
    /// oxygens just outside the list, `(r_c + s)(1 + 10⁻⁶)`, each water moved `0.75 s` towards the
    /// other, past `s/2`, so that the list is rebuilt, and the oxygens end `s/2` inside the cutoff.
    /// A rule of `s` an atom would not rebuild, and would miss them.
    #[test]
    fn two_atoms_closing_a_pair_together_rebuild_the_list() {
        let (rc, s) = (5.0 * ANGSTROM, 1.0 * ANGSTROM);
        let (plain, at) = two_waters((rc + s) * (1.0 + 1e-6));
        let listed = plain.clone().with_neighbour_list(s);
        same_bits(&plain, &listed, &at);
        assert!(list_of(&listed)
            .state
            .borrow()
            .entries
            .iter()
            .all(|e| *e & ATOM != 3));
        let end = moved(&at, 0.75 * s, -0.75 * s);
        assert!(end[3][0] - end[0][0] < rc - 0.49 * s);
        same_bits(&plain, &listed, &end);
        assert!(calls_pair(&listed, &end, 0, 3));
        assert_eq!(listed.neighbour_list().unwrap().builds, 2);
    }

    /// **The rule is `s/2`, not less and not more**: one water moved `s/2 + ε` rebuilds the list
    /// and `s/2 − ε` does not, along x and along a diagonal, so that what is compared is a length
    /// and not a component. And **a move by whole box lengths is no move**: one water given three
    /// images away along x and one back along y, the other one back along z, and nothing is
    /// rebuilt, the bits the cell list's at those positions.
    #[test]
    fn the_list_is_rebuilt_past_half_the_skin_and_not_before() {
        let (rc, s) = (5.0 * ANGSTROM, 1.0 * ANGSTROM);
        let eps = 1e-6 * ANGSTROM;
        let (plain, at) = two_waters(rc + 3.0 * s);
        for (sign, builds) in [(-1.0, 1), (1.0, 2)] {
            let listed = plain.clone().with_neighbour_list(s);
            same_bits(&plain, &listed, &at);
            let x = moved(&at, 0.5 * s + sign * eps, 0.0);
            same_bits(&plain, &listed, &x);
            assert_eq!(listed.neighbour_list().unwrap().builds, builds, "x, {sign}");
            let listed = plain.clone().with_neighbour_list(s);
            same_bits(&plain, &listed, &at);
            let c = (0.5 * s + sign * eps) / 3f64.sqrt();
            let diagonal: Vec<[f64; 3]> = at
                .iter()
                .enumerate()
                .map(|(i, p)| if i < 3 { p.map(|v| v + c) } else { *p })
                .collect();
            same_bits(&plain, &listed, &diagonal);
            assert_eq!(
                listed.neighbour_list().unwrap().builds,
                builds,
                "diagonal, {sign}"
            );
        }
        let listed = plain.clone().with_neighbour_list(s);
        same_bits(&plain, &listed, &at);
        let l = plain.cell().lengths();
        let images: Vec<[f64; 3]> = at
            .iter()
            .enumerate()
            .map(|(i, p)| {
                if i < 3 {
                    [p[0] + 3.0 * l[0], p[1] - l[1], p[2]]
                } else {
                    [p[0], p[1], p[2] - l[2]]
                }
            })
            .collect();
        same_bits(&plain, &listed, &images);
        assert_eq!(listed.neighbour_list().unwrap().builds, 1);
    }

    /// A new box or new atoms start the list again; a clone keeps a list of its own; and two
    /// fields that differ only in their lists' state are equal.
    #[test]
    fn a_new_box_or_new_atoms_start_the_list_again() {
        let (plain, at) = two_waters(6.0 * ANGSTROM);
        let listed = plain.clone().with_neighbour_list(ANGSTROM);
        same_bits(&plain, &listed, &at);
        let status = listed.neighbour_list().unwrap();
        assert_eq!((status.builds, status.evaluations), (1, 3));
        assert!(status.pairs > 0);
        let bigger = listed.with_cell(crate::PeriodicBox::cubic(31.0 * ANGSTROM));
        assert_eq!(bigger.neighbour_list().unwrap().builds, 0);
        let ion = crate::solvated::Ion::chloride();
        let more = listed.clone().with_ions(&[ion]);
        assert_eq!(more.neighbour_list().unwrap().builds, 0);
        let copy = listed.clone();
        copy.evaluate(&at);
        assert_eq!(copy.neighbour_list().unwrap().evaluations, 4);
        assert_eq!(listed.neighbour_list().unwrap().evaluations, 3);
        assert_eq!(copy, listed, "the state is a cache");
        assert_ne!(copy, plain);
        assert_eq!(listed.clone().without_neighbour_list(), plain);
        assert!(plain.neighbour_list().is_none());
    }

    /// **No pair the list is told to leave out reaches the evaluation, whatever it was built
    /// with**, which is what lets the force field's pair term skip its own lookup: two waters'
    /// six intramolecular pairs, the only pairs inside the cutoff, listed with nothing excluded;
    /// then asked to leave out one of them, the list is built again and does; then another one
    /// instead, as many exclusions as before, so that comparing their number would not do; then
    /// asked with nothing excluded again, built again, and gives them back. The same exclusions
    /// twice build nothing.
    #[test]
    fn other_exclusions_build_the_list_again() {
        let (plain, at) = two_waters(6.0 * ANGSTROM);
        let (cell, rc) = (plain.cell(), plain.ewald().parameters().cutoff);
        let list = NeighbourList::new(ANGSTROM);
        let calls = |ex: &Exclusions| {
            let mut v = Vec::new();
            list.for_each_pair(&cell, &at, rc, ex, |i, j, _, _| v.push((i, j)));
            v.sort_unstable();
            v
        };
        let none = Exclusions::from_pairs(at.len(), &[]);
        let one = Exclusions::from_pairs(at.len(), &[[0, 1]]);
        let all = calls(&none);
        assert_eq!(all, [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)]);
        assert_eq!(list.status().builds, 1);
        let other = Exclusions::from_pairs(at.len(), &[[0, 2]]);
        let fewer = calls(&one);
        assert_eq!(fewer, all[1..]);
        assert_eq!(list.status().builds, 2);
        let instead = calls(&other);
        assert_eq!(instead, [all[0], all[2], all[3], all[4], all[5]]);
        assert_eq!(list.status().builds, 3);
        assert_eq!(calls(&none), all);
        assert_eq!(list.status().builds, 4);
        assert_eq!(calls(&none), all);
        assert_eq!(list.status().builds, 4);
    }

    /// **A cutoff far below an atom's spacing is searched, not aborted**: 0.01 Å in a 9.3 Å box
    /// asked for more than 1 800³ cells uncapped, and now gets 128³, every pair found against the
    /// minimum image of every pair — none at 0.01 Å, and at 1.6 and 1.9 Å the O···H pairs of
    /// neighbouring waters. And W4's box keeps its cells, so its order.
    #[test]
    fn a_tiny_cutoff_is_capped_not_aborted() {
        let a = ANGSTROM;
        assert_eq!(
            cell_grid([62.35 * a, 65.47 * a, 74.82 * a], 9.0 * a),
            Some([13, 14, 16])
        );
        let water = WaterBox::lattice(3, 0x7E57);
        let cell = water.cell();
        assert_eq!(cell_grid(cell.lengths(), 0.01 * a), Some([128; 3]));
        assert!((2.0 * cell.lengths()[0] / (0.01 * a)).floor() > 1800.0);
        let field = water.force_field(4.6 * a, 1e-5);
        let at = water.positions();
        for cutoff in [0.01, 1.6, 1.9] {
            let rc = cutoff * a;
            let mut found = Vec::new();
            for_each_pair(&cell, at, rc, |i, j, _, _| {
                if !field.exclusions.contains(i, j) {
                    found.push((i, j));
                }
            });
            found.sort_unstable();
            let mut want = Vec::new();
            for i in 0..at.len() {
                for j in i + 1..at.len() {
                    let d = cell.minimum_image(sub(at[i], at[j]));
                    if !field.exclusions.contains(i, j)
                        && d[0] * d[0] + d[1] * d[1] + d[2] * d[2] < rc * rc
                    {
                        want.push((i, j));
                    }
                }
            }
            assert_eq!(found, want, "{cutoff} Å");
            assert_eq!(want.is_empty(), cutoff < 1.0, "{cutoff} Å: {}", want.len());
        }
    }

    /// **A capped grid finds every pair, where there are pairs to find**: 2 500 points in a
    /// 100 Å cube at `r_c` = 1.2 Å, 166³ = 4.57 × 10⁶ cells uncapped and 128³ capped, every pair
    /// the cell list calls the brute-force set over 27 images — not empty — and a list's calls the
    /// cell list's to the bit, through the capped grid's keys, before and after every point moves.
    #[test]
    fn a_capped_grid_finds_every_pair() {
        let a = ANGSTROM;
        let cell = crate::PeriodicBox::cubic(100.0 * a);
        let rc = 1.2 * a;
        assert!((2.0 * 100.0 / 1.2f64).floor() == 166.0);
        assert_eq!(cell_grid(cell.lengths(), rc), Some([128; 3]));
        let n = 2500;
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut at: Vec<[f64; 3]> = (0..n).map(|_| [0; 3].map(|_| 100.0 * a * next())).collect();
        let none = Exclusions::from_pairs(n, &[]);
        let list = NeighbourList::new(0.3 * a);
        for round in 0..2 {
            let mut cells = Vec::new();
            for_each_pair(&cell, &at, rc, |i, j, d, r2| {
                cells.push((i, j, d.map(f64::to_bits), r2.to_bits()))
            });
            let mut listed = Vec::new();
            list.for_each_pair(&cell, &at, rc, &none, |i, j, d, r2| {
                listed.push((i, j, d.map(f64::to_bits), r2.to_bits()))
            });
            assert_eq!(listed, cells, "round {round}");
            let mut found: Vec<(usize, usize)> = cells.iter().map(|c| (c.0, c.1)).collect();
            found.sort_unstable();
            let l = cell.lengths();
            let mut want = Vec::new();
            for i in 0..n {
                for j in i + 1..n {
                    let d0 = cell.minimum_image(sub(at[i], at[j]));
                    // Every other image is at least L − |d0_x| ≥ L/2 + r_c away on that axis.
                    let inner = (0..3).all(|x| d0[x].abs() <= 0.5 * l[x] - rc);
                    if inner && d0[0] * d0[0] + d0[1] * d0[1] + d0[2] * d0[2] >= rc * rc {
                        continue;
                    }
                    let nearest = (0..27)
                        .map(|m| {
                            let s = [m / 9, (m / 3) % 3, m % 3].map(|k| k as f64 - 1.0);
                            let d = [0, 1, 2].map(|x| d0[x] + s[x] * l[x]);
                            d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
                        })
                        .fold(f64::INFINITY, f64::min);
                    if nearest < rc * rc {
                        want.push((i, j));
                    }
                }
            }
            assert!(want.len() > 10, "round {round}: {} pairs", want.len());
            assert_eq!(found, want, "round {round}");
            // Every point moved by up to 0.05 Å an axis, at most 0.087 Å, inside half the skin
            // (0.15 Å): no build, and some points in new cells.
            let before: Vec<[usize; 3]> = at
                .iter()
                .map(|p| cell_of(cell.wrap(*p), l, [128; 3]))
                .collect();
            for p in at.iter_mut() {
                for x in p.iter_mut() {
                    *x += 0.1 * a * (next() - 0.5);
                }
            }
            let crossed = (0..n)
                .filter(|&i| cell_of(cell.wrap(at[i]), l, [128; 3]) != before[i])
                .count();
            assert!(crossed > 100, "{crossed} points changed cell");
        }
        assert_eq!(list.status().builds, 1);
    }

    /// **A build on threads is the build on one, and any cut of the atoms visits the same pairs**:
    /// the entries and their starts, and every call `(i, j, d, r²)` in order, the same bits on 2,
    /// 3, 7 and 16 threads cut into 2 to 10 000 spans as on one thread in one span — the cuts past
    /// the atoms giving spans of one atom or a few. Each build is asserted to have spawned more
    /// than one thread, but on WebAssembly. Three boxes: 27 waters at 3.5 Å, visited by cells and
    /// built over every pair, the build's radius of 5.5 Å being too wide for five cells across; at
    /// 4.6 Å, both over every pair; and 64 waters at 3.5 Å with a 0.6 Å skin, both by cells.
    #[test]
    fn a_build_on_threads_is_the_build_on_one() {
        let wasm = cfg!(target_family = "wasm");
        for (per_side, cutoff, skin, visit_by_cells, build_by_cells) in [
            (3, 3.5, 2.0, true, false),
            (3, 4.6, 2.0, false, false),
            (4, 3.5, 0.6, true, true),
        ] {
            let (cutoff, skin) = (cutoff * ANGSTROM, skin * ANGSTROM);
            let water = WaterBox::lattice(per_side, 0x3C0);
            let field = water.force_field(cutoff, 1e-8);
            let (cell, at, ex) = (field.cell(), water.positions(), field.skipped());
            let radius = (cutoff + skin) * (1.0 + MARGIN);
            assert_eq!(cell_grid(cell.lengths(), cutoff).is_some(), visit_by_cells);
            assert_eq!(cell_grid(cell.lengths(), radius).is_some(), build_by_cells);
            type Calls = Vec<(usize, usize, [u64; 3], u64)>;
            let run = |threads: usize, pieces: usize| -> (Calls, Vec<usize>, Vec<u64>, usize) {
                let list = NeighbourList::new(skin);
                let mut calls = Vec::new();
                let mut spans = 0;
                list.for_each_span(&cell, at, cutoff, ex, threads, pieces, |cut| {
                    spans = cut.len();
                    for span in cut {
                        span.visit(&mut |i, j, d, r2| {
                            calls.push((i, j, d.map(f64::to_bits), r2.to_bits()));
                        });
                    }
                });
                let s = list.state.borrow();
                assert_eq!(s.builds, 1);
                if threads > 1 && !wasm {
                    assert!(s.spawned > 1, "{threads} threads built on {}", s.spawned);
                } else {
                    assert_eq!(s.spawned, 0);
                }
                (calls, s.start.clone(), s.entries.clone(), spans)
            };
            let (calls, start, entries, spans) = run(1, 1);
            assert_eq!(spans, 1);
            assert!(calls.len() > 500, "{} pairs", calls.len());
            for (threads, pieces) in [(2, 2), (3, 24), (7, 56), (16, 81), (16, 10_000)] {
                let label = format!("{per_side}³ waters, {threads} threads, {pieces} pieces");
                let (c, st, en, sp) = run(threads, pieces);
                assert!(sp > 1 && sp <= pieces, "{label}: {sp} spans");
                assert_eq!(st, start, "{label}");
                assert_eq!(en, entries, "{label}");
                assert_eq!(c, calls, "{label}");
            }
        }
    }

    /// **Cut into spans, the list visits what it visits in one at every step of a trajectory**,
    /// across its builds and the steps between them, when atoms change cell and are keyed and
    /// sorted again span by span: three lists of their own, visited in one span, in 24 on three
    /// threads' builds and in 81 on sixteen's, every call the same bits in the same order at each
    /// of 120 steps of 125 waters at 4.6 Å, six half-cutoff cells across.
    #[test]
    fn a_list_cut_into_spans_visits_what_it_visits_whole_along_a_trajectory() {
        let (field, path) = trajectory(5, 4.6 * ANGSTROM, 120);
        let (cell, cutoff, ex) = (
            field.cell(),
            field.ewald().parameters().cutoff,
            field.skipped(),
        );
        assert!(cell_grid(cell.lengths(), cutoff).is_some());
        let skin = 1.0 * ANGSTROM;
        let lists = [(1, 1), (3, 24), (16, 81)].map(|(t, p)| (t, p, NeighbourList::new(skin)));
        let mut changed = 0;
        let mut previous: Option<Vec<[usize; 3]>> = None;
        for (k, at) in path.iter().enumerate() {
            let builds = lists[0].2.status().builds;
            let mut seen: Vec<Vec<Call>> = Vec::new();
            for (threads, pieces, list) in &lists {
                let mut calls = Vec::new();
                list.for_each_span(&cell, at, cutoff, ex, *threads, *pieces, |cut| {
                    for span in cut {
                        span.visit(&mut |i, j, d, r2| {
                            calls.push((i, j, d.map(f64::to_bits), r2.to_bits()));
                        });
                    }
                });
                seen.push(calls);
            }
            assert_eq!(seen[1], seen[0], "step {k}, 24 spans");
            assert_eq!(seen[2], seen[0], "step {k}, 81 spans");
            // Atoms that changed cell at a step that built nothing: keyed and sorted again in
            // place, span by span.
            let home = lists[0].2.state.borrow().home.clone();
            if let (Some(before), true) = (&previous, lists[0].2.status().builds == builds) {
                changed += home.iter().zip(before).filter(|(a, b)| a != b).count();
            }
            previous = Some(home);
        }
        let s = lists[0].2.status();
        assert!(s.builds >= 2, "{} builds", s.builds);
        assert!(s.builds < s.evaluations / 4, "{} builds", s.builds);
        assert!(changed > 100, "{changed} changes of cell");
        for (_, _, list) in &lists {
            assert_eq!(list.status().builds, s.builds);
        }
    }
}
