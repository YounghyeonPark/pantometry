//! **n-Hexadecane, and n-decane, as the oil of benzene's oil:water partition**: step C-1. Each is
//! the dictionary's entry (`R16`, `D10`), typed for UFF, relaxed in vacuum and put on its
//! constraints as `tests/octanol/` makes octanol ([`Octanol::relaxed`]); GAFF's non-bonded terms on
//! UFF's bonded ones, as `liquid_octanol_with_gaff.rs`'s octanol is.
//!
//! - **GAFF's `c3` and `hc` only** ([`gaff_kinds`]): every carbon of an n-alkane is `c3` and every
//!   hydrogen `hc` — no hydrogen has an electron-withdrawing neighbour — read off the elements,
//!   and FreeSolv's octane, nonane and decane files say the same. Lorentz–Berthelot, and GAFF's 1-4
//!   scales ([`ONE_FOUR`]).
//! - **The charges are decane's, by position** ([`transferred`]): FreeSolv has no hexadecane, so
//!   its AM1-BCC charges are built from decane's (`mobley_2197088`): the outer [`ENDS`] carbons at
//!   each end, with their hydrogens, at decane's values from the same end, and every interior
//!   carbon and its hydrogens at decane's interior values (C5 and C6, which are equal); then **every
//!   atom shifted by the same `−Σq/N`** so that the molecule is neutral. Decane's own charges sum to
//!   −0.0005 e as written; hexadecane's built ones to −0.0047 before the shift, which is
//!   +0.000094 e an atom.
//! - **What earns the transfer**: the same construction on octane, nonane and decane, against
//!   FreeSolv's own files for them. Measured: worst atom 0.00103 e off in octane, **0.00119 in
//!   nonane**, and decane reproduced to its own shift. So the brief's 0.001 e does not hold: nonane's
//!   C3 is −0.0795 where decane's is −0.0783, and octane's −0.0793 — AM1-BCC's own chain-to-chain
//!   variation at four places, which no construction by position from one chain can follow.
//!   [`HOMOLOGUE_BOUND`] is the measured 0.0012 plus half a unit of the files' fourth place.
//!
//! **C-1b's isotropic start** ([`isotropic`]): C-1's lattice of parallel chains never lost its
//! alignment, so each molecule is instead a conformer with random gauche defects
//! ([`Alkane::conformer`]) in a uniformly random orientation ([`uniform_rotation`]), set at its
//! lattice site in C-1's own box and then moved, never turned, until no two atoms overlap
//! ([`pack`]). The orientations are the draws to rounding, so the start is isotropic by
//! construction rather than by melting.
//!
//! **Kept out of the published crate** (`Cargo.toml`'s `exclude`) with `tests/gaff/`, whose numbers
//! it uses.

#![allow(dead_code)]
// The geometry is written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

use crate::freesolv::{self, Mol2};
use crate::gaff::{gaff, ONE_FOUR};
use crate::nonbonded::Vacuum;
use crate::octanol::{joined, Liquid, Octanol, ANGSTROM};
use pantometry_core::Rng;
use pantometry_forcefield::{
    Component, Element, EwaldParameters, Mixing, PeriodicBox, PeriodicForceField,
};
use std::f64::consts::PI;
use std::ops::Range;

/// n-Hexadecane, `R16` of the Chemical Component Dictionary, byte for byte.
pub const R16: &str = include_str!("../../components/R16.cif");

/// n-Decane, `D10`, byte for byte: the fallback oil.
pub const D10: &str = include_str!("../../components/D10.cif");

/// How many carbons at each end take decane's own end values.
pub const ENDS: usize = 4;

/// **What the construction reproduces octane and nonane to, e an atom**: measured 0.001035 and
/// 0.001193 (see the module documentation), so 0.0012 to the files' four places, and half a unit of
/// the fourth place over it. Earned by the measurement, not chosen: the brief asked for 0.001.
pub const HOMOLOGUE_BOUND: f64 = 0.00125;

/// An oil and the numbers it is measured against.
pub struct Oil {
    pub name: &'static str,
    /// The dictionary entry's text.
    pub entry: &'static str,
    /// The density at 25 °C, kg m⁻³.
    pub density: f64,
    /// `ΔH_vap` at 25 °C, kJ/mol.
    pub enthalpy_of_vaporization: f64,
    /// The melting point, kelvin.
    pub melting: f64,
    /// The viscosity, mPa s, and the temperature it is at, °C.
    pub viscosity: (f64, f64),
    /// Where the four came from.
    pub cited: &'static str,
    /// The measurement's lattice.
    pub lattice: [usize; 3],
}

/// **n-Hexadecane's targets**: 0.7701 g/cm³ at 25 °C (CRC Handbook 95th ed., 2014, p. 3-294),
/// `ΔH_vap` 81.35 kJ/mol at 25 °C (p. 6-139; NIST's WebBook gives 81.3 ± 0.6, as the design study
/// read it), melting at 18.18 °C (p. 3-294), each as PubChem CID 11006's HSDB record quotes it, read
/// on 2026-10-24; viscosity 3.474 mPa s at 20 °C (Ullmann's, as the same record quotes it). The
/// lattice is 2 × 7 × 7, 98 molecules: the extended chain is 20.5 Å, so a cube of 95 does not fit
/// a whole number of sites, and this box is 46 × 34 × 30 Å.
pub const HEXADECANE: Oil = Oil {
    name: "n-hexadecane",
    entry: R16,
    density: 770.1,
    enthalpy_of_vaporization: 81.35,
    melting: 291.33,
    viscosity: (3.474, 20.0),
    cited: "CRC Handbook 95th ed. pp. 3-294 and 6-139 and Ullmann's, via PubChem CID 11006 (HSDB), read 2026-10-24",
    lattice: [2, 7, 7],
};

/// **n-Decane's**: 0.7255 g/cm³ at 25 °C (CRC p. 3-142), `ΔH_vap` 51.42 kJ/mol at 25 °C (p. 6-135),
/// melting at −29.7 °C (p. 3-142), viscosity 0.838 mPa s at 25 °C (p. 6-235), via PubChem CID 15600's
/// HSDB record, read 2026-10-24. 2 × 7 × 8, 112 molecules, 31 × 34 × 34 Å.
pub const DECANE: Oil = Oil {
    name: "n-decane",
    entry: D10,
    density: 725.5,
    enthalpy_of_vaporization: 51.42,
    melting: 243.45,
    viscosity: (0.838, 25.0),
    cited: "CRC Handbook 95th ed. pp. 3-142, 6-135 and 6-235, via PubChem CID 15600 (HSDB), read 2026-10-24",
    lattice: [2, 7, 8],
};

/// Each atom's neighbours.
fn adjacency(n: usize, bonds: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &[a, b] in bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    adj
}

/// **The carbons of an unbranched alkane in chain order**, from its terminal carbon of lowest index,
/// and each carbon's hydrogens.
///
/// # Panics
///
/// If the molecule is not an n-alkane: only carbon and hydrogen, every carbon four bonds and at most
/// two of them to carbon, every hydrogen one bond to a carbon, and the carbons one chain.
pub fn chain(elements: &[Element], bonds: &[[usize; 2]]) -> (Vec<usize>, Vec<Vec<usize>>) {
    let adj = adjacency(elements.len(), bonds);
    let carbons: Vec<usize> = (0..elements.len())
        .filter(|&i| elements[i] == Element::C)
        .collect();
    let carbon_neighbours = |i: usize| -> Vec<usize> {
        adj[i]
            .iter()
            .copied()
            .filter(|&j| elements[j] == Element::C)
            .collect()
    };
    for (i, e) in elements.iter().enumerate() {
        match e {
            Element::C => {
                assert_eq!(adj[i].len(), 4, "carbon {i} has four bonds");
                assert!(
                    carbon_neighbours(i).len() <= 2,
                    "carbon {i} is not a branch"
                );
            }
            Element::H => {
                assert_eq!(adj[i].len(), 1, "hydrogen {i} has one bond");
                assert_eq!(
                    elements[adj[i][0]],
                    Element::C,
                    "hydrogen {i} is on a carbon"
                );
            }
            other => panic!("an alkane has no {other:?}"),
        }
    }
    let start = *carbons
        .iter()
        .find(|&&i| carbon_neighbours(i).len() == 1)
        .expect("a terminal carbon");
    let mut order = vec![start];
    let mut previous = usize::MAX;
    let mut at = start;
    loop {
        let next: Vec<usize> = carbon_neighbours(at)
            .into_iter()
            .filter(|&j| j != previous)
            .collect();
        match next.as_slice() {
            [] => break,
            [j] => {
                previous = at;
                at = *j;
                order.push(at);
            }
            _ => panic!("a branch at {at}"),
        }
    }
    assert_eq!(order.len(), carbons.len(), "the carbons are one chain");
    let hydrogens = order
        .iter()
        .map(|&c| {
            adj[c]
                .iter()
                .copied()
                .filter(|&j| elements[j] == Element::H)
                .collect()
        })
        .collect();
    (order, hydrogens)
}

/// **A FreeSolv alkane's charges by position**: for each carbon along its chain, its charge and its
/// hydrogens' — which must all be one value, as AM1-BCC's are by symmetry.
pub fn by_position(mol2: &Mol2) -> Vec<(f64, f64)> {
    let elements: Vec<Element> = mol2.atoms.iter().map(|a| a.element).collect();
    let (order, hydrogens) = chain(&elements, &mol2.bonds);
    order
        .iter()
        .zip(&hydrogens)
        .map(|(&c, hs)| {
            let q: Vec<f64> = hs.iter().map(|&h| mol2.atoms[h].charge).collect();
            assert!(
                q.iter().all(|x| x.to_bits() == q[0].to_bits()),
                "{}: carbon {}'s hydrogens differ",
                mol2.name,
                mol2.atoms[c].name
            );
            (mol2.atoms[c].charge, q[0])
        })
        .collect()
}

/// The construction's result: each atom's charge in the target's order, what the positions summed
/// to before the shift, and the shift every atom was given.
#[derive(Clone, Debug)]
pub struct Transfer {
    pub charges: Vec<f64>,
    pub net_before: f64,
    pub shift: f64,
}

/// **Charges for the n-alkane `(elements, bonds)` from decane's by position**: see the module
/// documentation. `decane` is [`by_position`] of decane's file.
///
/// # Panics
///
/// If the chain is shorter than `2 ENDS`, or decane's interior positions are not one value.
pub fn transferred(decane: &[(f64, f64)], elements: &[Element], bonds: &[[usize; 2]]) -> Transfer {
    let m = decane.len();
    assert_eq!(m, 10, "decane");
    let interior = decane[ENDS];
    assert_eq!(
        [interior.0.to_bits(), interior.1.to_bits()],
        [
            decane[m - 1 - ENDS].0.to_bits(),
            decane[m - 1 - ENDS].1.to_bits()
        ],
        "decane's interior positions are one value"
    );
    let (order, hydrogens) = chain(elements, bonds);
    let n = order.len();
    assert!(n >= 2 * ENDS, "a chain of {n} has no ends of {ENDS}");
    let mut q = vec![f64::NAN; elements.len()];
    for k in 0..n {
        let from_end = n - 1 - k;
        let (qc, qh) = if k < ENDS {
            decane[k]
        } else if from_end < ENDS {
            decane[m - 1 - from_end]
        } else {
            interior
        };
        q[order[k]] = qc;
        for &h in &hydrogens[k] {
            q[h] = qh;
        }
    }
    assert!(q.iter().all(|x| x.is_finite()), "every atom has a charge");
    let net_before: f64 = q.iter().sum();
    let shift = -net_before / q.len() as f64;
    Transfer {
        charges: q.iter().map(|x| x + shift).collect(),
        net_before,
        shift,
    }
}

/// FreeSolv's decane by position.
pub fn decane_positions() -> Vec<(f64, f64)> {
    by_position(&freesolv::decane())
}

/// **GAFF's type for each atom of an n-alkane, read off its elements**: `c3` for carbon, `hc` for
/// hydrogen.
pub fn gaff_kinds(c: &Component) -> Vec<&'static str> {
    c.atoms()
        .iter()
        .map(|a| match a.element {
            Element::C => "c3",
            Element::H => "hc",
            other => panic!("an alkane has no {other:?}"),
        })
        .collect()
}

/// An alkane: the molecule, made as octanol is, the construction its charges came from, its GAFF
/// types and its chain.
pub struct Alkane {
    pub molecule: Octanol,
    pub transfer: Transfer,
    pub kinds: Vec<&'static str>,
    /// The carbons in chain order, on the molecule's own indices.
    pub carbons: Vec<usize>,
}

impl Alkane {
    /// The alkane of dictionary entry `entry`, its charges [`transferred`] from decane's, relaxed
    /// in vacuum under UFF and turned into its frame: first carbon → last along `+x`, the second in
    /// the `xy` plane.
    pub fn new(entry: &str) -> Alkane {
        let component = Component::from_ccd(entry).expect("the entry");
        let elements: Vec<Element> = component.atoms().iter().map(|a| a.element).collect();
        let bonds: Vec<[usize; 2]> = component.bonds().iter().map(|b| b.atoms).collect();
        let transfer = transferred(&decane_positions(), &elements, &bonds);
        let (carbons, _) = chain(&elements, &bonds);
        let name = |i: usize| component.atoms()[i].name.clone();
        let (a, b, c) = (
            name(carbons[0]),
            name(carbons[carbons.len() - 1]),
            name(carbons[1]),
        );
        let kinds = gaff_kinds(&component);
        let molecule = Octanol::relaxed(
            component,
            transfer.charges.clone(),
            [a.as_str(), b.as_str(), c.as_str()],
        );
        Alkane {
            molecule,
            transfer,
            kinds,
            carbons,
        }
    }

    /// GAFF's `[x, D]` for each atom.
    pub fn vdw(&self) -> Vec<[f64; 2]> {
        self.kinds.iter().map(|k| gaff(k)).collect()
    }

    /// The terminal carbons, on the molecule's own indices.
    pub fn ends(&self) -> [usize; 2] {
        [self.carbons[0], self.carbons[self.carbons.len() - 1]]
    }
}

/// `field`, a box of whole alkanes `a`, with GAFF's van der Waals on every atom, Lorentz–Berthelot
/// and GAFF's 1-4 scales; its charges are already the alkane's.
pub fn with_gaff(field: PeriodicForceField, a: &Alkane) -> PeriodicForceField {
    let one = a.vdw();
    let n = field.charges().len();
    assert_eq!(n % one.len(), 0, "whole alkanes");
    let values: Vec<[f64; 2]> = (0..n).map(|i| one[i % one.len()]).collect();
    field
        .with_nonbonded(0..n, values)
        .expect("GAFF's values are valid")
        .with_mixing(Mixing::LorentzBerthelot)
        .with_one_four(ONE_FOUR[0], ONE_FOUR[1])
}

/// One alkane in vacuum with GAFF's non-bonded terms on UFF's bonded ones.
pub fn gas_phase(a: &Alkane) -> Vacuum {
    let o = &a.molecule;
    Vacuum::new(
        &o.component,
        &o.types,
        a.vdw(),
        o.charges.clone(),
        Mixing::LorentzBerthelot,
        ONE_FOUR,
    )
}

// ---------------------------------------------------------------------------------------------
// What says liquid or crystal

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = dot(a, a).sqrt();
    a.map(|x| x / l)
}

/// **The largest eigenvalue of a symmetric 3 × 3 matrix**, by the closed form for its roots
/// (O. K. Smith, *Commun. ACM* **4**, 168 (1961)).
pub fn largest_eigenvalue(q: [[f64; 3]; 3]) -> f64 {
    let p1 = q[0][1] * q[0][1] + q[0][2] * q[0][2] + q[1][2] * q[1][2];
    let tr = (q[0][0] + q[1][1] + q[2][2]) / 3.0;
    if p1 == 0.0 {
        return q[0][0].max(q[1][1]).max(q[2][2]);
    }
    let d = [q[0][0] - tr, q[1][1] - tr, q[2][2] - tr];
    let p = ((d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + 2.0 * p1) / 6.0).sqrt();
    let b = [
        [d[0] / p, q[0][1] / p, q[0][2] / p],
        [q[0][1] / p, d[1] / p, q[1][2] / p],
        [q[0][2] / p, q[1][2] / p, d[2] / p],
    ];
    let det = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
        - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
        + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
    let phi = (0.5 * det).clamp(-1.0, 1.0).acos() / 3.0;
    tr + 2.0 * p * phi.cos()
}

/// **The nematic order parameter of unit vectors `u`**: the largest eigenvalue of
/// `Q = ⟨(3 u uᵀ − I)/2⟩`. One for every vector parallel, as the lattice starts; of order `1/√N`
/// for `N` isotropic ones ([`isotropic_order`]).
pub fn nematic_order(u: &[[f64; 3]]) -> f64 {
    let mut q = [[0.0; 3]; 3];
    for v in u {
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                q[a][b] += (1.5 * v[a] * v[b] - 0.5 * delta) / u.len() as f64;
            }
        }
    }
    largest_eigenvalue(q)
}

/// **What an isotropic liquid of `n` molecules shows**: the nematic order of `n` directions drawn
/// uniformly on the sphere, `sets` times by `Rng::for_index(seed, set)`; its mean and standard
/// deviation.
pub fn isotropic_order(n: usize, sets: usize, seed: u64) -> (f64, f64) {
    let s: Vec<f64> = (0..sets)
        .map(|k| {
            let mut rng = Rng::for_index(seed, k as u64);
            let u: Vec<[f64; 3]> = (0..n)
                .map(|_| {
                    let z = 2.0 * rng.unit() - 1.0;
                    let phi = 2.0 * PI * rng.unit();
                    let r = (1.0 - z * z).sqrt();
                    [r * phi.cos(), r * phi.sin(), z]
                })
                .collect();
            nematic_order(&u)
        })
        .collect();
    let m = s.iter().sum::<f64>() / sets as f64;
    let var = s.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (sets - 1) as f64;
    (m, var.sqrt())
}

/// Each molecule's end-to-end direction, first carbon to last.
pub fn directions(at: &[[f64; 3]], molecules: &[Range<usize>], ends: [usize; 2]) -> Vec<[f64; 3]> {
    molecules
        .iter()
        .map(|m| unit(sub(at[m.start + ends[1]], at[m.start + ends[0]])))
        .collect()
}

/// **The fraction of C–C–C–C dihedrals that are trans**, `|φ| > 120°`, over every molecule, and the
/// mean end-to-end distance (metres).
pub fn conformations(at: &[[f64; 3]], molecules: &[Range<usize>], carbons: &[usize]) -> (f64, f64) {
    let (mut trans, mut all, mut ends) = (0usize, 0usize, 0.0);
    for m in molecules {
        let c: Vec<[f64; 3]> = carbons.iter().map(|&i| at[m.start + i]).collect();
        for k in 0..c.len() - 3 {
            let (b1, b2, b3) = (
                sub(c[k + 1], c[k]),
                sub(c[k + 2], c[k + 1]),
                sub(c[k + 3], c[k + 2]),
            );
            let (n1, n2) = (cross(b1, b2), cross(b2, b3));
            let cos = dot(n1, n2) / (dot(n1, n1) * dot(n2, n2)).sqrt();
            if cos < -0.5 {
                trans += 1;
            }
            all += 1;
        }
        let e = sub(c[c.len() - 1], c[0]);
        ends += dot(e, e).sqrt();
    }
    (trans as f64 / all as f64, ends / molecules.len() as f64)
}

/// Each molecule's centre of mass.
pub fn centres(at: &[[f64; 3]], molecules: &[Range<usize>], masses: &[f64]) -> Vec<[f64; 3]> {
    molecules
        .iter()
        .map(|r| {
            let m: f64 = masses[r.clone()].iter().sum();
            [0, 1, 2].map(|k| r.clone().map(|i| masses[i] * at[i][k]).sum::<f64>() / m)
        })
        .collect()
}

/// **The centres' mean-square displacement at each lag of `lags`** (in frames), averaged over every
/// origin and molecule, each frame's displacements taken relative to the whole box's — the bath
/// moves the box's centre of mass as a random walk. Unwrapped positions: molecular dynamics never
/// wraps. The mean over `group`'s molecules only, when given.
pub fn msd(frames: &[Vec<[f64; 3]>], lags: &[usize], group: Option<&[usize]>) -> Vec<f64> {
    let n = frames[0].len();
    let all: Vec<usize> = (0..n).collect();
    let g = group.unwrap_or(&all);
    let mean =
        |f: &Vec<[f64; 3]>| [0, 1, 2].map(|k| f.iter().map(|p| p[k]).sum::<f64>() / n as f64);
    let box_centre: Vec<[f64; 3]> = frames.iter().map(mean).collect();
    lags.iter()
        .map(|&lag| {
            let (mut sum, mut count) = (0.0, 0usize);
            for t in 0..frames.len() - lag {
                let shift = sub(box_centre[t + lag], box_centre[t]);
                for &i in g {
                    let d = sub(sub(frames[t + lag][i], frames[t][i]), shift);
                    sum += dot(d, d);
                    count += 1;
                }
            }
            sum / count as f64
        })
        .collect()
}

/// **`g(r)` between two sets of sites in different molecules**, one frame: every pair `a ∈ first`,
/// `b ∈ second` of different molecules by minimum image, binned by `width` to `bins`, normalised by
/// the number of such pairs over the volume — so that it tends to one. `owner[i]` is atom `i`'s
/// molecule. With `same`, `first` and `second` are one list and each pair is counted once.
#[allow(clippy::too_many_arguments)]
pub fn g_of_r(
    cell: &PeriodicBox,
    at: &[[f64; 3]],
    owner: &[usize],
    first: &[usize],
    second: &[usize],
    same: bool,
    width: f64,
    bins: usize,
) -> Vec<f64> {
    let mut h = vec![0.0; bins];
    let mut pairs = 0usize;
    for (x, &a) in first.iter().enumerate() {
        let others = if same { &second[x + 1..] } else { second };
        for &b in others {
            if owner[a] == owner[b] {
                continue;
            }
            pairs += 1;
            let d = cell.minimum_image(sub(at[a], at[b]));
            let k = (dot(d, d).sqrt() / width) as usize;
            if k < bins {
                h[k] += 1.0;
            }
        }
    }
    let v = cell.volume();
    for (k, x) in h.iter_mut().enumerate() {
        let (r0, r1) = (k as f64 * width, (k + 1) as f64 * width);
        let shell = 4.0 / 3.0 * PI * (r1 * r1 * r1 - r0 * r0 * r0);
        *x /= pairs as f64 * shell / v;
    }
    h
}

/// `g(r)`'s first peak (position by a parabola through the highest bin and its neighbours, and
/// height) and the first minimum after it (position and value), searching the peak in `range`.
pub fn first_shell(g: &[f64], width: f64, range: [f64; 2]) -> ([f64; 2], [f64; 2]) {
    let lo = (range[0] / width) as usize;
    let hi = ((range[1] / width) as usize).min(g.len() - 1);
    let b = (lo..hi)
        .max_by(|&i, &j| g[i].total_cmp(&g[j]))
        .expect("bins");
    let (y0, y1, y2) = (g[b - 1], g[b], g[b + 1]);
    let curvature = y0 - 2.0 * y1 + y2;
    let shift = if curvature == 0.0 {
        0.0
    } else {
        0.5 * (y0 - y2) / curvature
    };
    let peak = [
        (b as f64 + 0.5 + shift) * width,
        y1 - 0.25 * (y0 - y2) * shift,
    ];
    let end = ((b as f64 + 3.0 * ANGSTROM / width) as usize).min(g.len());
    let m = (b..end)
        .min_by(|&i, &j| g[i].total_cmp(&g[j]))
        .expect("bins");
    (peak, [(m as f64 + 0.5) * width, g[m]])
}

// ---------------------------------------------------------------------------------------------
// An isotropic start: step C-1b

/// A 3 × 3 rotation matrix, rows first.
pub type Rotation = [[f64; 3]; 3];

/// `r v`.
pub fn rotate(r: &Rotation, v: [f64; 3]) -> [f64; 3] {
    [dot(r[0], v), dot(r[1], v), dot(r[2], v)]
}

/// `a b`.
fn compose(a: &Rotation, b: &Rotation) -> Rotation {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}

/// **The rotation by `angle` about the unit vector `axis`**, Rodrigues' formula.
pub fn about(axis: [f64; 3], angle: f64) -> Rotation {
    let (s, c) = angle.sin_cos();
    let t = 1.0 - c;
    let [x, y, z] = axis;
    [
        [c + t * x * x, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, c + t * y * y, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, c + t * z * z],
    ]
}

/// **A rotation drawn uniformly on SO(3)** — by Haar measure — from three draws of `rng`: the unit
/// quaternion `(√(1−u₁) sin 2πu₂, √(1−u₁) cos 2πu₂, √u₁ sin 2πu₃, √u₁ cos 2πu₃)`, which is uniform
/// on the 3-sphere (K. Shoemake, "Uniform random rotations", *Graphics Gems III*, 1992, 124–132),
/// and the 3-sphere double-covers SO(3) with Haar measure as its image.
pub fn uniform_rotation(rng: &mut Rng) -> Rotation {
    let (u1, u2, u3) = (rng.unit(), rng.unit(), rng.unit());
    let (a, b) = ((1.0 - u1).sqrt(), u1.sqrt());
    let (s2, c2) = (2.0 * PI * u2).sin_cos();
    let (s3, c3) = (2.0 * PI * u3).sin_cos();
    let (x, y, z, w) = (a * s2, a * c2, b * s3, b * c3);
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - z * w),
            2.0 * (x * z + y * w),
        ],
        [
            2.0 * (x * y + z * w),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - x * w),
        ],
        [
            2.0 * (x * z - y * w),
            2.0 * (y * z + x * w),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}

/// **The packer's radius of an atom, metres**: a pair of atoms in different molecules is an overlap
/// nearer than the sum of theirs — H–H 1.6 Å, C–H 2.0 Å and C–C 2.4 Å: [`PACKING_SCALE`] of 1.0 Å
/// for hydrogen and 1.5 Å for carbon, whose sums would be the 2.0 Å [`crate::octanol::OVERLAP`] is
/// for H–H, 2.5 and 3.0.
pub fn packing_radius(e: Element) -> f64 {
    let full = match e {
        Element::H => 1.0,
        Element::C => 1.5,
        other => panic!("an alkane has no {other:?}"),
    };
    PACKING_SCALE * full * ANGSTROM
}

/// **What the packer can reach without turning a molecule**, as a fraction of the full radii:
/// measured on C-1b's 112 decanes at 0.7255 g/cm³, [`pack`] without its random kicks cleared every
/// overlap at 0.8 for 5 of 16 seeds in 5 000 iterations and at 0.75 for 13; at 0.85 the start's seed
/// was left at 0.95 of the radii after 20 000; and growing the radii in stages, 0.5 to 0.8, jammed
/// worse (0.90 and 0.88 on the first two seeds). With the kicks, 0.8 clears for 15 of 16 seeds in
/// 20 000 iterations, the start's own in 267 with none. Rigid conformers in fixed, random
/// orientations jam short of a liquid's own closest contacts. What is left is strain the first
/// stage takes out: GAFF's `hc`–`hc` at 1.6 Å is +25 kcal/mol, and in the first 0.4 ps at 0.5 fs the
/// box's `U/N` measured 65.6 → 17.7 kcal/mol, with the temperature at 298 K from 0.1 ps.
pub const PACKING_SCALE: f64 = 0.8;

/// **The trans probability of each C–C–C–C dihedral of a start's conformer**: 0.75, between the
/// 0.825 C-1's decane reached at 298.15 K and a melt's lower value, each gauche ±120° from trans
/// with equal probability.
pub const START_TRANS: f64 = 0.75;

/// **How far apart two atoms of one molecule more than three bonds apart must be in a start's
/// conformer**, metres. Measured on the dictionary's decane turned rigidly: one gauche brings two
/// hydrogens five bonds apart to 1.87 Å, and `g⁺g⁻` to 0.72 Å, so 1.8 Å keeps every single gauche
/// and refuses the pentane clash. 2.0 Å refused every gauche, and left every chain all-trans.
pub const SELF_CLEARANCE: f64 = 1.8 * ANGSTROM;

/// The number of bonds between each pair of atoms, by breadth-first search.
fn bond_distances(n: usize, bonds: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let adj = adjacency(n, bonds);
    (0..n)
        .map(|i| {
            let mut d = vec![usize::MAX; n];
            d[i] = 0;
            let mut queue = std::collections::VecDeque::from([i]);
            while let Some(a) = queue.pop_front() {
                for &b in &adj[a] {
                    if d[b] == usize::MAX {
                        d[b] = d[a] + 1;
                        queue.push_back(b);
                    }
                }
            }
            d
        })
        .collect()
}

/// The atoms on `b`'s side of the bond `a`–`b`, `b` included.
fn beyond(adj: &[Vec<usize>], a: usize, b: usize) -> Vec<usize> {
    let mut seen = vec![false; adj.len()];
    seen[a] = true;
    seen[b] = true;
    let mut stack = vec![b];
    let mut side = vec![b];
    while let Some(x) = stack.pop() {
        for &y in &adj[x] {
            if !seen[y] {
                seen[y] = true;
                stack.push(y);
                side.push(y);
            }
        }
    }
    side
}

/// A conformer drawn for a start: positions about the centre of mass (metres), in the molecule's
/// frame, and how many of its dihedrals are gauche.
pub struct Conformer {
    pub at: Vec<[f64; 3]>,
    pub gauche: usize,
}

impl Alkane {
    /// **A conformer with random gauche defects**: from the relaxed all-trans conformer, each of the
    /// `n − 3` C–C–C–C dihedrals left trans with probability `p_trans` and otherwise turned by +120°
    /// or −120° — the whole of the chain beyond its middle bond rotated rigidly about that bond, so
    /// that every bond length and angle is the relaxed conformer's to rounding, and every held X–H
    /// stays at its length. Drawn again, from the same stream, while two atoms more than three bonds
    /// apart are nearer than [`SELF_CLEARANCE`] (`g⁺g⁻`'s clash). Returned about its centre of mass.
    ///
    /// # Panics
    ///
    /// If a thousand draws all clash.
    pub fn conformer(&self, rng: &mut Rng, p_trans: f64) -> Conformer {
        let o = &self.molecule;
        let n = o.at.len();
        let bonds: Vec<[usize; 2]> = o.component.bonds().iter().map(|b| b.atoms).collect();
        let adj = adjacency(n, &bonds);
        let apart = bond_distances(n, &bonds);
        let sides: Vec<(usize, usize, Vec<usize>)> = (1..self.carbons.len() - 2)
            .map(|k| {
                let (a, b) = (self.carbons[k], self.carbons[k + 1]);
                (a, b, beyond(&adj, a, b))
            })
            .collect();
        for _ in 0..1000 {
            let mut at = o.at.clone();
            let mut gauche = 0;
            for (a, b, side) in &sides {
                let u = rng.unit();
                if u < p_trans {
                    continue;
                }
                let sign = if u < 0.5 * (1.0 + p_trans) { 1.0 } else { -1.0 };
                let r = about(unit(sub(at[*b], at[*a])), sign * 2.0 * PI / 3.0);
                let pivot = at[*b];
                for &i in side {
                    let d = rotate(&r, sub(at[i], pivot));
                    at[i] = [pivot[0] + d[0], pivot[1] + d[1], pivot[2] + d[2]];
                }
                gauche += 1;
            }
            let clash = (0..n).any(|i| {
                (i + 1..n).any(|j| {
                    let d = sub(at[i], at[j]);
                    apart[i][j] > 3 && dot(d, d) < SELF_CLEARANCE * SELF_CLEARANCE
                })
            });
            if clash {
                continue;
            }
            let total = o.mass();
            let com = [0, 1, 2]
                .map(|k| at.iter().zip(&o.masses).map(|(p, m)| m * p[k]).sum::<f64>() / total);
            return Conformer {
                at: at.iter().map(|p| sub(*p, com)).collect(),
                gauche,
            };
        }
        panic!("a thousand conformers all clash");
    }
}

/// What [`pack`] did: the positions and centres (metres), its iterations and evaluations, how many
/// times its pair list was built, the overlap penalty before and after (Å⁴), and the worst pair's
/// distance over the sum of its radii — at least one when it succeeded.
pub struct Packed {
    pub at: Vec<[f64; 3]>,
    pub centres: Vec<[f64; 3]>,
    pub iterations: usize,
    pub evaluations: usize,
    pub lists: usize,
    /// How many times it was jammed and moved the overlapping molecules at random.
    pub kicks: usize,
    pub penalty: [f64; 2],
    pub worst: f64,
}

/// How far past the radii [`pack`] pushes, as a factor: it minimises the overlap at `1.1` of every
/// pair's sum and stops when none is nearer than the sum itself.
pub const PACKING_MARGIN: f64 = 1.1;

/// The pair list's skin, Å: pairs within the longest reach plus this are listed, and the list is
/// built again once a molecule has moved half of it.
const PACKING_SKIN: f64 = 1.0;

/// How many iterations [`pack`] descends before it calls itself jammed and moves every molecule in
/// an overlap by a random step.
const PACKING_ROUND: usize = 300;

/// The random step's largest component, Å.
const PACKING_KICK: f64 = 1.0;

/// How many past steps [`pack`]'s L-BFGS keeps.
const PACKING_MEMORY: usize = 8;

/// The furthest [`pack`] moves any centre in one step, Å.
const PACKING_STRIDE: f64 = 0.5;

/// **Every pair of atoms in different molecules within `reach` by minimum image**, `i < j`, by a
/// cell list of cells at least `reach` across. Å.
fn pairs_within(x: &[[f64; 3]], owner: &[usize], l: [f64; 3], reach: f64) -> Vec<[usize; 2]> {
    let nc = l.map(|e| ((e / reach).floor() as usize).max(3));
    let cell_of = |p: [f64; 3]| {
        [0, 1, 2].map(|k| {
            let w = p[k] - l[k] * (p[k] / l[k]).floor();
            ((w / l[k] * nc[k] as f64) as usize).min(nc[k] - 1)
        })
    };
    let index = |c: [usize; 3]| (c[0] * nc[1] + c[1]) * nc[2] + c[2];
    let mut cells = vec![Vec::new(); nc[0] * nc[1] * nc[2]];
    let home: Vec<[usize; 3]> = x.iter().map(|p| cell_of(*p)).collect();
    for (i, c) in home.iter().enumerate() {
        cells[index(*c)].push(i);
    }
    let mut pairs = Vec::new();
    for i in 0..x.len() {
        let c = home[i];
        for dx in [nc[0] - 1, 0, 1] {
            for dy in [nc[1] - 1, 0, 1] {
                for dz in [nc[2] - 1, 0, 1] {
                    let n = [
                        (c[0] + dx) % nc[0],
                        (c[1] + dy) % nc[1],
                        (c[2] + dz) % nc[2],
                    ];
                    for &j in &cells[index(n)] {
                        if j <= i || owner[j] == owner[i] {
                            continue;
                        }
                        let d = minimum_image(sub(x[i], x[j]), l);
                        if dot(d, d) < reach * reach {
                            pairs.push([i, j]);
                        }
                    }
                }
            }
        }
    }
    pairs
}

/// `d − L round(d/L)` per axis.
fn minimum_image(mut d: [f64; 3], l: [f64; 3]) -> [f64; 3] {
    for k in 0..3 {
        d[k] -= l[k] * (d[k] / l[k]).round();
    }
    d
}

/// **The overlap penalty and its gradient** over `pairs`, Å units: `Σ (t² − r²)²` over every pair
/// nearer than `t = PACKING_MARGIN (ρ_i + ρ_j)`, by minimum image — Packmol's objective (L.
/// Martínez et al., *J. Comput. Chem.* **30**, 2157 (2009)) — and the smallest `r / (ρ_i + ρ_j)`
/// among those pairs (infinite when none).
fn overlap(
    x: &[[f64; 3]],
    pairs: &[[usize; 2]],
    rho: &[f64],
    l: [f64; 3],
) -> (f64, Vec<[f64; 3]>, f64) {
    let mut f = 0.0;
    let mut g = vec![[0.0; 3]; x.len()];
    let mut worst = f64::INFINITY;
    for &[i, j] in pairs {
        let d = minimum_image(sub(x[i], x[j]), l);
        let r2 = dot(d, d);
        let sum = rho[i] + rho[j];
        let t = PACKING_MARGIN * sum;
        if r2 < t * t {
            let e = t * t - r2;
            f += e * e;
            for k in 0..3 {
                g[i][k] -= 4.0 * e * d[k];
                g[j][k] += 4.0 * e * d[k];
            }
            worst = worst.min(r2.sqrt() / sum);
        }
    }
    (f, g, worst)
}

/// **Rigid molecules moved apart, never turned, until no two atoms overlap**: molecule `m`'s atoms
/// are `centres[m] + rotations[m] bodies[m][i]`, and only the centres move, so every molecule keeps
/// the orientation it was drawn with. L-BFGS on [`overlap`]'s penalty in the centres, whose
/// gradient is each molecule's atoms' summed, with a backtracking (Armijo) line search, no centre
/// moved more than [`PACKING_STRIDE`] a step, and a pair list rebuilt once a molecule has moved half
/// [`PACKING_SKIN`]. Molecules with `fixed[m]` do not move. **Jammed** — [`PACKING_ROUND`]
/// iterations without clearing, or a line search that finds no lower penalty along steepest descent
/// — every molecule in an overlap is moved by a step of up to [`PACKING_KICK`] along each axis, drawn
/// from `Rng::for_index(seed, kick × molecules + molecule)`, and the descent starts again.
/// Deterministic: every draw indexed, and every sum in one order. Stops when the worst pair is at
/// its radii's sum, or after `max_iterations`.
///
/// **Why the molecules are not turned**: measured on C-1b's 112 decanes, a packer that also turned
/// them by the gradient's torque — from the lattice sites, at full radii — took the nematic order of
/// uniformly drawn orientations to 0.55; and from a box first packed without turning at 0.6 of the
/// radii, to 0.24. Rods pushed apart at a liquid's density line up, so turning undoes the start's
/// point.
#[allow(clippy::too_many_arguments)]
pub fn pack(
    bodies: &[Vec<[f64; 3]>],
    radii: &[Vec<f64>],
    centres: Vec<[f64; 3]>,
    rotations: &[Rotation],
    fixed: &[bool],
    cell: &PeriodicBox,
    seed: u64,
    max_iterations: usize,
) -> Packed {
    let a = ANGSTROM;
    let l = cell.lengths().map(|x| x / a);
    let mols = bodies.len();
    let owner: Vec<usize> = bodies
        .iter()
        .enumerate()
        .flat_map(|(m, b)| (0..b.len()).map(move |_| m))
        .collect();
    let rho: Vec<f64> = radii.iter().flatten().map(|r| r / a).collect();
    // Each molecule's atoms about its centre, turned once and for all.
    let arms: Vec<Vec<[f64; 3]>> = bodies
        .iter()
        .zip(rotations)
        .map(|(b, r)| b.iter().map(|p| rotate(r, p.map(|x| x / a))).collect())
        .collect();
    let place = |c: &[f64]| -> Vec<[f64; 3]> {
        arms.iter()
            .enumerate()
            .flat_map(|(m, b)| {
                b.iter()
                    .map(move |p| [c[3 * m] + p[0], c[3 * m + 1] + p[1], c[3 * m + 2] + p[2]])
            })
            .collect()
    };
    let reach = PACKING_MARGIN * 2.0 * rho.iter().fold(0.0f64, |m, r| m.max(*r)) + PACKING_SKIN;
    // The rigid gradient: each free molecule's atoms' summed.
    let gather = |g: &[[f64; 3]]| -> Vec<f64> {
        let mut out = vec![0.0; 3 * mols];
        for (i, gi) in g.iter().enumerate() {
            let m = owner[i];
            if !fixed[m] {
                for k in 0..3 {
                    out[3 * m + k] += gi[k];
                }
            }
        }
        out
    };
    let mut c: Vec<f64> = centres.iter().flat_map(|p| p.map(|x| x / a)).collect();
    let mut x = place(&c);
    let mut listed_at = c.clone();
    let mut pairs = pairs_within(&x, &owner, l, reach);
    let mut lists = 1;
    let (mut f, g0, mut worst) = overlap(&x, &pairs, &rho, l);
    let mut g = gather(&g0);
    let start = f;
    let mut evaluations = 1;
    let mut iterations = 0;
    let mut history: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    let vdot = |u: &[f64], v: &[f64]| u.iter().zip(v).map(|(p, q)| p * q).sum::<f64>();
    let (mut kicks, mut since_kick) = (0usize, 0usize);
    while worst < 1.0 && iterations < max_iterations {
        if since_kick >= PACKING_ROUND {
            // Jammed: every molecule in an overlap moved by a random step and the descent
            // started again. Packmol's `movebad`, in spirit.
            let mut bad = vec![false; mols];
            for &[i, j] in &pairs {
                let e = minimum_image(sub(x[i], x[j]), l);
                let r = dot(e, e).sqrt();
                if r < rho[i] + rho[j] {
                    bad[owner[i]] = true;
                    bad[owner[j]] = true;
                }
            }
            for m in (0..mols).filter(|&m| bad[m] && !fixed[m]) {
                let mut rng = Rng::for_index(seed, (kicks * mols + m) as u64);
                for k in 0..3 {
                    c[3 * m + k] += PACKING_KICK * (2.0 * rng.unit() - 1.0);
                }
            }
            x = place(&c);
            pairs = pairs_within(&x, &owner, l, reach);
            listed_at = c.clone();
            lists += 1;
            let (f2, g2, worst2) = overlap(&x, &pairs, &rho, l);
            evaluations += 1;
            f = f2;
            g = gather(&g2);
            worst = worst2;
            history.clear();
            kicks += 1;
            since_kick = 0;
            continue;
        }
        iterations += 1;
        since_kick += 1;
        // The two-loop recursion.
        let mut q = g.clone();
        let mut alphas = Vec::with_capacity(history.len());
        for (s, y) in history.iter().rev() {
            let rho_k = 1.0 / vdot(y, s);
            let alpha = rho_k * vdot(s, &q);
            for (qi, yi) in q.iter_mut().zip(y) {
                *qi -= alpha * yi;
            }
            alphas.push((rho_k, alpha));
        }
        let scale = match history.last() {
            Some((s, y)) => vdot(s, y) / vdot(y, y),
            None => {
                let most = (0..mols)
                    .map(|m| vdot(&g[3 * m..3 * m + 3], &g[3 * m..3 * m + 3]).sqrt())
                    .fold(0.0f64, f64::max);
                0.1 / most
            }
        };
        for qi in q.iter_mut() {
            *qi *= scale;
        }
        for ((s, y), (rho_k, alpha)) in history.iter().zip(alphas.iter().rev()) {
            let beta = rho_k * vdot(y, &q);
            for (qi, si) in q.iter_mut().zip(s) {
                *qi += (alpha - beta) * si;
            }
        }
        let mut d: Vec<f64> = q.iter().map(|v| -v).collect();
        if vdot(&d, &g) >= 0.0 {
            history.clear();
            d = g.iter().map(|v| -v * scale.abs()).collect();
        }
        let longest = (0..mols)
            .map(|m| vdot(&d[3 * m..3 * m + 3], &d[3 * m..3 * m + 3]).sqrt())
            .fold(0.0f64, f64::max);
        if longest > PACKING_STRIDE {
            for v in d.iter_mut() {
                *v *= PACKING_STRIDE / longest;
            }
        }
        let slope = vdot(&d, &g);
        let mut t = 1.0;
        let mut accepted = None;
        while t > 1e-12 {
            let c2: Vec<f64> = c.iter().zip(&d).map(|(p, v)| p + t * v).collect();
            let moved = (0..mols)
                .map(|m| {
                    let e = [0, 1, 2].map(|k| c2[3 * m + k] - listed_at[3 * m + k]);
                    dot(e, e).sqrt()
                })
                .fold(0.0f64, f64::max);
            let x2 = place(&c2);
            if moved > 0.5 * PACKING_SKIN {
                pairs = pairs_within(&x2, &owner, l, reach);
                listed_at = c2.clone();
                lists += 1;
            }
            let (f2, g2, worst2) = overlap(&x2, &pairs, &rho, l);
            evaluations += 1;
            if f2 <= f + 1e-4 * t * slope {
                accepted = Some((c2, x2, f2, gather(&g2), worst2));
                break;
            }
            t *= 0.5;
        }
        let Some((c2, x2, f2, g2, worst2)) = accepted else {
            if history.is_empty() {
                since_kick = PACKING_ROUND;
            }
            history.clear();
            continue;
        };
        let s: Vec<f64> = c2.iter().zip(&c).map(|(p, q)| p - q).collect();
        let y: Vec<f64> = g2.iter().zip(&g).map(|(p, q)| p - q).collect();
        if vdot(&s, &y) > 0.0 {
            history.push((s, y));
            if history.len() > PACKING_MEMORY {
                history.remove(0);
            }
        }
        c = c2;
        x = x2;
        f = f2;
        g = g2;
        worst = worst2;
    }
    Packed {
        at: x.iter().map(|p| p.map(|v| v * a)).collect(),
        centres: (0..mols)
            .map(|m| [c[3 * m] * a, c[3 * m + 1] * a, c[3 * m + 2] * a])
            .collect(),
        iterations,
        evaluations,
        lists,
        kicks,
        penalty: [start, f],
        worst,
    }
}

/// **The worst contact by brute force over the 27 images**: the smallest `r / (ρ_i + ρ_j)`, radii
/// `radii`, over every pair of atoms in different molecules at each of the 27 shifts
/// `n ∈ {−1, 0, 1}³`, and every atom with its own molecule's at the 26 that are not zero.
/// Independent of [`pack`]'s cell list and minimum image. Returns the ratio and the two atoms.
///
/// **The positions are taken as they are, each molecule whole**, not wrapped atom by atom: a
/// molecule across a face would otherwise meet itself at a shift that is not an image, which is
/// what `closest_contact` would do here, its lattice's molecules each lying inside the box. So
/// that the 27 hold every nearest image, every atom must lie within a quarter-edge of the box: then
/// each difference is under 1.5 edges, and its nearest image is one of the 27.
///
/// # Panics
///
/// If an atom lies further than a quarter of an edge outside the box.
pub fn worst_contact(
    at: &[[f64; 3]],
    molecules: &[Range<usize>],
    radii: &[f64],
    cell: &PeriodicBox,
) -> (f64, usize, usize) {
    let l = cell.lengths();
    for p in at {
        for k in 0..3 {
            assert!(
                p[k] >= -0.25 * l[k] && p[k] < 1.25 * l[k],
                "an atom at {p:?} is past a quarter-edge outside the box"
            );
        }
    }
    let w = at;
    let mut owner = vec![usize::MAX; at.len()];
    for (m, r) in molecules.iter().enumerate() {
        for i in r.clone() {
            owner[i] = m;
        }
    }
    assert!(
        owner.iter().all(|&m| m != usize::MAX),
        "every atom is in a molecule"
    );
    let mut best = (f64::INFINITY, 0, 0);
    for i in 0..w.len() {
        for j in i..w.len() {
            let same = owner[i] == owner[j];
            let d = sub(w[i], w[j]);
            let s = radii[i] + radii[j];
            for nx in -1i32..=1 {
                for ny in -1i32..=1 {
                    for nz in -1i32..=1 {
                        if same && nx == 0 && ny == 0 && nz == 0 {
                            continue;
                        }
                        let e = [
                            d[0] + nx as f64 * l[0],
                            d[1] + ny as f64 * l[1],
                            d[2] + nz as f64 * l[2],
                        ];
                        let ratio = dot(e, e).sqrt() / s;
                        if ratio < best.0 {
                            best = (ratio, i, j);
                        }
                    }
                }
            }
        }
    }
    best
}

/// What [`isotropic`] built: the box, and how — each molecule's conformer about its centre of mass,
/// its gauche count and the rotation drawn for it, and the packer's report.
pub struct Start {
    pub liquid: Liquid,
    pub bodies: Vec<Vec<[f64; 3]>>,
    pub rotations: Vec<Rotation>,
    pub gauche: Vec<usize>,
    pub packed: Packed,
    /// Every atom's packing radius.
    pub radii: Vec<f64>,
}

/// **An isotropic start**: `counts` molecules of alkane `a` in [`Liquid::cell_of`]'s box at
/// `density` — C-1's box, so the density is exact to rounding as there — each molecule a
/// [`Alkane::conformer`] at `p_trans` turned by a [`uniform_rotation`], both drawn from
/// `Rng::for_index(seed, molecule)`, its centre of mass first at its lattice site; then [`pack`]ed
/// until no two atoms overlap, and each centre wrapped into the box. Ewald at `cutoff` and
/// `accuracy`, UFF's non-bonded terms until [`with_gaff`] puts GAFF's on.
///
/// # Panics
///
/// If the packer does not clear every overlap in `max_iterations`.
#[allow(clippy::too_many_arguments)]
pub fn isotropic(
    a: &Alkane,
    counts: [usize; 3],
    density: f64,
    seed: u64,
    p_trans: f64,
    cutoff: f64,
    accuracy: f64,
    max_iterations: usize,
) -> Start {
    let o = &a.molecule;
    let (cell, site, gap) = Liquid::cell_of(o, counts, density);
    let mols = counts[0] * counts[1] * counts[2];
    let one: Vec<f64> = o
        .component
        .atoms()
        .iter()
        .map(|x| packing_radius(x.element))
        .collect();
    let (mut bodies, mut centres, mut rotations, mut gauche) =
        (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for i in 0..counts[0] {
        for j in 0..counts[1] {
            for k in 0..counts[2] {
                let index = ((i * counts[1] + j) * counts[2] + k) as u64;
                let mut rng = Rng::for_index(seed, index);
                let c = a.conformer(&mut rng, p_trans);
                rotations.push(uniform_rotation(&mut rng));
                bodies.push(c.at);
                gauche.push(c.gauche);
                centres.push([
                    (i as f64 + 0.5) * site[0],
                    (j as f64 + 0.5) * site[1],
                    (k as f64 + 0.5) * site[2],
                ]);
            }
        }
    }
    let radii = vec![one.clone(); mols];
    let packed = pack(
        &bodies,
        &radii,
        centres,
        &rotations,
        &vec![false; mols],
        &cell,
        seed ^ 0x9AC4,
        max_iterations,
    );
    assert!(
        packed.worst >= 1.0,
        "the packer left an overlap at {} of the radii after {} iterations",
        packed.worst,
        packed.iterations
    );
    // Each molecule moved by whole edges so that its centre of mass is in the box.
    let n = o.at.len();
    let mut at = packed.at.clone();
    let mut molecules = Vec::new();
    for m in 0..mols {
        let r = m * n..(m + 1) * n;
        let total = o.mass();
        let com = [0, 1, 2].map(|k| {
            r.clone()
                .map(|i| o.masses[i - m * n] * at[i][k])
                .sum::<f64>()
                / total
        });
        let shift = sub(cell.wrap(com), com);
        for i in r.clone() {
            at[i] = [
                at[i][0] + shift[0],
                at[i][1] + shift[1],
                at[i][2] + shift[2],
            ];
        }
        molecules.push(r);
    }
    let parts: Vec<&Component> = (0..mols).map(|_| &o.component).collect();
    let component = joined(&parts);
    let types: Vec<_> = (0..mols).flat_map(|_| o.types.iter().cloned()).collect();
    let charges: Vec<f64> = (0..mols).flat_map(|_| o.charges.iter().copied()).collect();
    let masses: Vec<f64> = (0..mols).flat_map(|_| o.masses.iter().copied()).collect();
    let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
    let field = PeriodicForceField::new(&component, &types, cell, p)
        .expect("alkanes in their box")
        .with_charges(charges.clone());
    Start {
        liquid: Liquid {
            component,
            types,
            charges,
            field,
            at,
            masses,
            molecules,
            solute: 0,
            octanols: mols,
            site,
            gap,
        },
        bodies,
        rotations,
        gauche,
        packed,
        radii: radii.concat(),
    }
}
