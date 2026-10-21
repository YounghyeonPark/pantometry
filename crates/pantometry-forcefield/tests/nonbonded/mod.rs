//! **What the non-bonded overrides are checked against, written out**: a molecule's 1-2, 1-3 and
//! 1-4 pairs by breadth-first search over its bonds — not the force field's walk over paths — and
//! its non-bonded energy and forces in vacuum as a plain sum over every pair, by either mixing
//! rule, with the 1-4 pairs scaled. And a set of invented van der Waals values to override with,
//! which are no force field's.
//!
//! Used by `the_nonbonded_overrides_against_closed_forms.rs` and, with GAFF's values, by
//! `liquid_octanol_with_gaff.rs`.

#![allow(dead_code)]
// The pair sums are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::uff::{UffType, KCAL_PER_MOL};
use pantometry_forcefield::{
    Component, Element, EwaldParameters, ForceField, Mixing, PeriodicBox, PeriodicForceField,
    Potential,
};
use std::f64::consts::PI;

pub const ANGSTROM: f64 = 1e-10;

/// **Invented van der Waals, `[x, D]` in UFF's form**: round numbers chosen to be unlike UFF's and
/// unlike each other, so that the two mixing rules part (`x` from 1.2 to 3.9 Å) — and no force
/// field's. A hydrogen on an oxygen or a nitrogen gets a small one, as an alcohol's does in many
/// force fields; every other hydrogen another.
pub fn invented(c: &Component) -> Vec<[f64; 2]> {
    let mut neighbour = vec![None; c.atoms().len()];
    for b in c.bonds() {
        let [i, j] = b.atoms;
        neighbour[i] = Some(c.atoms()[j].element);
        neighbour[j] = Some(c.atoms()[i].element);
    }
    c.atoms()
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let (x, d) = match a.element {
                Element::C => (3.9, 0.090),
                Element::O => (3.3, 0.190),
                Element::N => (3.6, 0.160),
                Element::H => match neighbour[i] {
                    Some(Element::O) | Some(Element::N) => (1.2, 0.010),
                    _ => (2.7, 0.020),
                },
                other => panic!("no invented value for {other:?}"),
            };
            [x * ANGSTROM, d * KCAL_PER_MOL]
        })
        .collect()
}

/// UFF's own `[x, D]` for each type: what a field without overrides carries.
pub fn uff_values(types: &[UffType]) -> Vec<[f64; 2]> {
    types
        .iter()
        .map(|t| {
            let p = t.parameters();
            [p.vdw_distance, p.vdw_energy]
        })
        .collect()
}

/// The bonds of `c` as index pairs.
pub fn bonds(c: &Component) -> Vec<[usize; 2]> {
    c.bonds().iter().map(|b| b.atoms).collect()
}

/// **The shortest path through the bonds between every two atoms**, by breadth-first search from
/// each atom, `None` for atoms in different molecules: `[i][j]`.
pub fn graph_distances(n: usize, bonds: &[[usize; 2]]) -> Vec<Vec<Option<usize>>> {
    let mut adjacent = vec![Vec::new(); n];
    for &[i, j] in bonds {
        adjacent[i].push(j);
        adjacent[j].push(i);
    }
    (0..n)
        .map(|start| {
            let mut d = vec![None; n];
            d[start] = Some(0);
            let mut queue = std::collections::VecDeque::from([start]);
            while let Some(a) = queue.pop_front() {
                let da = d[a].expect("queued atoms are reached");
                for &b in &adjacent[a] {
                    if d[b].is_none() {
                        d[b] = Some(da + 1);
                        queue.push_back(b);
                    }
                }
            }
            d
        })
        .collect()
}

/// Every pair `[i, j]`, `i < j`, whose shortest path is exactly three bonds, sorted.
pub fn one_four_by_search(n: usize, bonds: &[[usize; 2]]) -> Vec<[usize; 2]> {
    let d = graph_distances(n, bonds);
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            if d[i][j] == Some(3) {
                out.push([i, j]);
            }
        }
    }
    out
}

/// One pair's van der Waals by the mixing rule written out: `D_ij [(x_ij/r)¹² − 2 (x_ij/r)⁶]` and
/// its `dE/dr`, with `x_ij` the geometric or arithmetic mean and `D_ij` the geometric.
pub fn pair_vdw(a: [f64; 2], b: [f64; 2], mixing: Mixing, r: f64) -> (f64, f64) {
    let x = match mixing {
        Mixing::Geometric => (a[0] * b[0]).sqrt(),
        Mixing::LorentzBerthelot => 0.5 * (a[0] + b[0]),
    };
    let d = (a[1] * b[1]).sqrt();
    let s = x / r;
    let s2 = s * s;
    let s6 = s2 * s2 * s2;
    let s12 = s6 * s6;
    (d * (s12 - 2.0 * s6), 12.0 * d * (s6 - s12) / r)
}

/// `(2π/V) Σ_i Σ_j D_ij [x_ij¹²/(9 r_c⁹) − 2 x_ij⁶/(3 r_c³)]` and A&T's `P_tail`, written out over
/// every ordered pair of `a` with `b`, by `mixing`, and the sums of the two parts' magnitudes for
/// the rounding.
pub fn tail_by_pairs(
    a: &[[f64; 2]],
    b: &[[f64; 2]],
    mixing: Mixing,
    rc: f64,
    v: f64,
) -> ([f64; 2], [f64; 2]) {
    let (rc3, mut t12, mut t6) = (rc * rc * rc, 0.0, 0.0);
    for p in a {
        for q in b {
            let x = match mixing {
                Mixing::Geometric => (p[0] * q[0]).sqrt(),
                Mixing::LorentzBerthelot => 0.5 * (p[0] + q[0]),
            };
            let d = (p[1] * q[1]).sqrt();
            let x2 = x * x;
            let x6 = x2 * x2 * x2;
            t12 += d * x6 * x6;
            t6 += d * x6;
        }
    }
    let rc9 = rc3 * rc3 * rc3;
    let (e12, e6) = (
        2.0 * PI / v * t12 / (9.0 * rc9),
        2.0 * PI / v * 2.0 * t6 / (3.0 * rc3),
    );
    let (p12, p6) = (
        2.0 * PI / (3.0 * v * v) * 4.0 * t12 / (3.0 * rc9),
        2.0 * PI / (3.0 * v * v) * 4.0 * t6 / rc3,
    );
    ([e12 - e6, p12 - p6], [e12 + e6, p12 + p6])
}

/// A molecule's non-bonded terms in vacuum, by part, and the forces they give.
pub struct Nonbonded {
    /// Van der Waals, joules.
    pub van_der_waals: f64,
    /// Coulomb, joules.
    pub coulomb: f64,
    /// The sum of every term's magnitude, for the rounding of the two above.
    pub magnitude: f64,
    /// How many pairs.
    pub pairs: usize,
    /// The force on each atom, newtons.
    pub forces: Vec<[f64; 3]>,
}

/// **Molecules in vacuum with another force field's non-bonded terms on UFF's bonded ones**: the
/// bonded terms are UFF's, from [`PeriodicForceField::bonded`], and the non-bonded every pair not
/// 1-2 or 1-3 by [`graph_distances`], with no cutoff and no image, its van der Waals by
/// [`pair_vdw`] and its Coulomb `k_e q_i q_j / r`, a 1-4 pair's scaled by `scale`.
#[derive(Clone)]
pub struct Vacuum {
    pub bonded: ForceField,
    pub vdw: Vec<[f64; 2]>,
    pub charges: Vec<f64>,
    pub mixing: Mixing,
    /// Every pair that is neither 1-2 nor 1-3, with its `[s_v, s_c]`.
    pub pairs: Vec<([usize; 2], [f64; 2])>,
}

impl Vacuum {
    pub fn new(
        component: &Component,
        types: &[UffType],
        vdw: Vec<[f64; 2]>,
        charges: Vec<f64>,
        mixing: Mixing,
        scale: [f64; 2],
    ) -> Vacuum {
        // Any box serves: only its bonded terms are taken.
        let cell = PeriodicBox::cubic(100.0 * ANGSTROM);
        let p = EwaldParameters::for_accuracy(&cell, 10.0 * ANGSTROM, 1e-3);
        let bonded = PeriodicForceField::new(component, types, cell, p)
            .expect("supported")
            .bonded()
            .expect("a solute")
            .clone();
        let n = component.atoms().len();
        let d = graph_distances(n, &bonds(component));
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                match d[i][j] {
                    Some(1) | Some(2) => {}
                    Some(3) => pairs.push(([i, j], scale)),
                    _ => pairs.push(([i, j], [1.0, 1.0])),
                }
            }
        }
        Vacuum {
            bonded,
            vdw,
            charges,
            mixing,
            pairs,
        }
    }

    /// The non-bonded terms at `at`.
    pub fn nonbonded(&self, at: &[[f64; 3]]) -> Nonbonded {
        let mut out = Nonbonded {
            van_der_waals: 0.0,
            coulomb: 0.0,
            magnitude: 0.0,
            pairs: self.pairs.len(),
            forces: vec![[0.0; 3]; at.len()],
        };
        for &([i, j], [sv, sc]) in &self.pairs {
            let d = [0, 1, 2].map(|k| at[i][k] - at[j][k]);
            let r = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let (v, dv) = pair_vdw(self.vdw[i], self.vdw[j], self.mixing, r);
            let c = sc * COULOMB * self.charges[i] * self.charges[j] / r;
            out.van_der_waals += sv * v;
            out.coulomb += c;
            out.magnitude += (sv * v).abs() + c.abs();
            // dE/dr of the pair, then −dE/dr along d on i and the opposite on j.
            let de_dr = sv * dv - c / r;
            for k in 0..3 {
                let f = -de_dr * d[k] / r;
                out.forces[i][k] += f;
                out.forces[j][k] -= f;
            }
        }
        out
    }
}

impl Potential for Vacuum {
    fn energy_and_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let b = self.bonded.evaluate(at);
        let nb = self.nonbonded(at);
        for i in 0..at.len() {
            for k in 0..3 {
                forces[i][k] = b.forces[i][k] + nb.forces[i][k];
            }
        }
        b.energy.total + nb.van_der_waals + nb.coulomb
    }
}
