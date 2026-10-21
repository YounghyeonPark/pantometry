//! **1-octanol, and a box of it**, for step 4c of the explicit-solvent track: benzene's partition
//! between octanol and water, log P. The molecule is the dictionary's `OC9`, typed for UFF, given
//! charges, relaxed in vacuum and put on its constraints; the box is a lattice of that one
//! conformer at the liquid's density, with no two atoms of different molecules (or of a molecule
//! and its own image) closer than the lattice's gap.
//!
//! - **One component for the whole liquid** ([`joined`]): the crate's periodic force field takes a
//!   solute as one [`Component`] and water as rigid TIP3P, and has no third kind of molecule, so a
//!   box of 192 octanols is written as one dictionary entry of 192 disconnected molecules, each
//!   atom named `<name>_<k>`, and read back by [`Component::from_ccd`]. Every term of the field then
//!   applies as it stands: each molecule is made whole through its own bonds, its 1-2 and 1-3 pairs
//!   are excluded, and every other pair, inside a molecule or between two, is a cut-off UFF pair
//!   with Ewald's electrostatics.
//! - **The lattice** ([`Liquid::lattice`]): the relaxed conformer in its own frame — centre of mass
//!   at the origin, C1 → O along `+x`, the carbon plane `xy` — on an orthorhombic lattice of
//!   `counts` sites, each site's cell its span along each axis plus one common gap `g`, solved so
//!   that the cell's volume is the molecule's share at the density ([`cell_for`]). Each site turns
//!   its molecule end for end, about its long axis, both or neither, by two draws of
//!   `Rng::for_index(seed, site)`; every one of the four is a proper rotation that keeps the
//!   carbon plane `xy`, so the gaps stay `g` to within the conformer's asymmetry.
//! - **The density is exact to rounding**: two edges of the cell are the spans plus `g`, and the
//!   third is the share's volume divided by their product.

#![allow(dead_code)]
// The geometry is written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

use pantometry_core::Rng;
use pantometry_forcefield::uff::{self, UffType};
use pantometry_forcefield::{
    qeq, Component, Element, EwaldParameters, ForceField, MolecularDynamics, PeriodicBox,
    PeriodicForceField, Shake,
};
use std::ops::Range;

pub const ANGSTROM: f64 = 1e-10;
pub const FS: f64 = 1e-15;

/// Octan-1-ol, `OC9` of the Chemical Component Dictionary, byte for byte.
pub const OC9: &str = include_str!("../../components/OC9.cif");

/// **1-octanol's density at 25 °C, kg m⁻³**: 0.8262 g/cm³, the CRC Handbook's (95th ed., 2014,
/// p. 3-426) as PubChem's HSDB record for CID 957 quotes it, read there on 2026-10-21. The
/// commonly quoted 0.827 g/cm³ is 0.1% above it, and is nearer the 20 °C value; PubChem's JECFA
/// record gives a specific-gravity range of 0.822–0.830.
pub const DENSITY: f64 = 826.2;

/// The atoms nearer than this, metres, are an overlap: an H–H pair of UFF's there is
/// `0.044 [(2.886/2)¹² − 2 (2.886/2)⁶]` = 2.8 kcal/mol, and the lattice's gap is 2.3 Å.
pub const OVERLAP: f64 = 2.0e-10;

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

/// The index of the atom called `name`.
pub fn index(c: &Component, name: &str) -> usize {
    c.atoms()
        .iter()
        .position(|a| a.name == name)
        .unwrap_or_else(|| panic!("no atom {name}"))
}

/// 1-octanol alone: its component, types, charges, vacuum force field and a relaxed geometry on
/// its constraints, in its own frame.
pub struct Octanol {
    pub component: Component,
    pub types: Vec<UffType>,
    pub charges: Vec<f64>,
    pub vacuum: ForceField,
    /// Relaxed in vacuum, every X–H put on UFF's natural length along itself, then turned into
    /// the molecule's frame: centre of mass at the origin, C1 (`CAA`) → O (`OAB`) along `+x`,
    /// C1, C2 and O in the `xy` plane with C2 at `+y`.
    pub at: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
}

impl Octanol {
    /// 1-octanol with `charges`, or QEq's at the dictionary's ideal geometry when `None`.
    pub fn new(charges: Option<Vec<f64>>) -> Octanol {
        let component = Component::from_ccd(OC9).expect("OC9");
        let types = uff::assign(&component);
        let ideal: Vec<[f64; 3]> = component.atoms().iter().map(|a| a.at).collect();
        let charges =
            charges.unwrap_or_else(|| qeq::charges(&component, &ideal).expect("QEq").charges);
        let vacuum = ForceField::new(&component, &types)
            .expect("UFF")
            .with_charges(charges.clone());
        let mut at = ideal;
        // To 10⁻⁴ kcal mol⁻¹ Å⁻¹ on every atom.
        vacuum.minimise(&mut at, 5000, 1e-4 * uff::KCAL_PER_MOL / ANGSTROM);
        let masses: Vec<f64> = component.atoms().iter().map(|a| a.element.mass()).collect();
        let shake = Shake::to_hydrogen(&component, &vacuum);
        let before = at.clone();
        shake.constrain_positions(&masses, &vec![false; at.len()], &before, &mut at);
        // The frame. A rotation and a translation keep every length to rounding.
        let total: f64 = masses.iter().sum();
        let com =
            [0, 1, 2].map(|k| at.iter().zip(&masses).map(|(p, m)| m * p[k]).sum::<f64>() / total);
        let (c1, c2, o) = (
            index(&component, "CAA"),
            index(&component, "CAC"),
            index(&component, "OAB"),
        );
        let u = unit(sub(at[o], at[c1]));
        let w = unit(cross(u, sub(at[c2], at[c1])));
        let v = cross(w, u);
        let at = at
            .iter()
            .map(|p| {
                let d = sub(*p, com);
                [dot(d, u), dot(d, v), dot(d, w)]
            })
            .collect();
        Octanol {
            component,
            types,
            charges,
            vacuum,
            at,
            masses,
        }
    }

    /// The eighteen X–H constraints, on the molecule's own indices.
    pub fn shake(&self) -> Shake {
        Shake::to_hydrogen(&self.component, &self.vacuum)
    }

    /// The molecule's mass, kilograms.
    pub fn mass(&self) -> f64 {
        self.masses.iter().sum()
    }

    /// The conformer's extent along each axis of its frame, `[lowest, highest]`, metres.
    pub fn extent(&self) -> [[f64; 2]; 3] {
        [0, 1, 2].map(|k| {
            self.at
                .iter()
                .fold([f64::INFINITY, f64::NEG_INFINITY], |e, p| {
                    [e[0].min(p[k]), e[1].max(p[k])]
                })
        })
    }
}

/// **The lattice's cell for a molecule of spans `span` (metres) and a share `volume` (m³)**: the
/// three edges and the common gap `g`. `g` is the root of `Π (span_k + g) = volume` by bisection,
/// and the third edge is then `volume / (a b)`, so the product is the volume to rounding whatever
/// the bisection left.
///
/// # Panics
///
/// If the molecule's spans alone exceed the volume.
pub fn cell_for(span: [f64; 3], volume: f64) -> ([f64; 3], f64) {
    let f = |g: f64| (span[0] + g) * (span[1] + g) * (span[2] + g);
    assert!(f(0.0) < volume, "the molecule does not fit its share");
    let (mut lo, mut hi) = (0.0, volume.cbrt());
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < volume {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let (a, b) = (span[0] + lo, span[1] + lo);
    ([a, b, volume / (a * b)], lo)
}

/// **One component of several molecules**: each part's atoms and bonds in turn, every atom named
/// `<name>_<k>` for the part's place `k`, read back by [`Component::from_ccd`]. The coordinates
/// written are each part's own and mean nothing: the box keeps its positions apart.
pub fn joined(parts: &[&Component]) -> Component {
    let yn = |b: bool| if b { "Y" } else { "N" };
    let mut s = String::from(
        "data_LIQ\n_chem_comp.id LIQ\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (k, c) in parts.iter().enumerate() {
        for a in c.atoms() {
            s += &format!(
                "LIQ {}_{k} {} {} {} {} {} {}\n",
                a.name,
                a.element.symbol(),
                a.charge,
                yn(a.aromatic),
                a.at[0] / ANGSTROM,
                a.at[1] / ANGSTROM,
                a.at[2] / ANGSTROM
            );
        }
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (k, c) in parts.iter().enumerate() {
        for b in c.bonds() {
            let order = match b.order.count() {
                1 => "SING",
                2 => "DOUB",
                _ => "TRIP",
            };
            s += &format!(
                "{}_{k} {}_{k} {order} {}\n",
                c.atoms()[b.atoms[0]].name,
                c.atoms()[b.atoms[1]].name,
                yn(b.aromatic)
            );
        }
    }
    Component::from_ccd(&s).expect("the joined entry parses")
}

/// A molecule placed in the box before the octanols: its component, types, charges and
/// positions (metres, already where it goes).
pub struct Solute<'a> {
    pub component: &'a Component,
    pub types: &'a [UffType],
    pub charges: &'a [f64],
    pub at: Vec<[f64; 3]>,
    /// Every lattice octanol with an atom nearer than this to a solute atom, by minimum image, is
    /// left out.
    pub clearance: f64,
}

/// A box of 1-octanol, and perhaps one solute first: see the module documentation.
pub struct Liquid {
    pub component: Component,
    pub types: Vec<UffType>,
    pub charges: Vec<f64>,
    pub field: PeriodicForceField,
    pub at: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
    /// Each molecule's atoms: the solute's first, if any, then each octanol's.
    pub molecules: Vec<Range<usize>>,
    /// How many atoms the solute has, zero for none.
    pub solute: usize,
    /// How many octanols.
    pub octanols: usize,
    /// The lattice's cell, metres, and its gap.
    pub site: [f64; 3],
    pub gap: f64,
}

impl Liquid {
    /// `counts[0] × counts[1] × counts[2]` octanols at `density` (kg m⁻³), turned by `seed`;
    /// Ewald at `cutoff` and `accuracy`.
    pub fn lattice(
        octanol: &Octanol,
        counts: [usize; 3],
        density: f64,
        seed: u64,
        cutoff: f64,
        accuracy: f64,
    ) -> Liquid {
        Liquid::build(None, octanol, counts, density, seed, cutoff, accuracy)
    }

    /// The same lattice with `solute` first, less every octanol it overlaps.
    #[allow(clippy::too_many_arguments)]
    pub fn around(
        solute: Solute,
        octanol: &Octanol,
        counts: [usize; 3],
        density: f64,
        seed: u64,
        cutoff: f64,
        accuracy: f64,
    ) -> Liquid {
        Liquid::build(
            Some(solute),
            octanol,
            counts,
            density,
            seed,
            cutoff,
            accuracy,
        )
    }

    /// The lattice's box for `counts` sites at `density`, and the site's cell and gap.
    pub fn cell_of(
        octanol: &Octanol,
        counts: [usize; 3],
        density: f64,
    ) -> (PeriodicBox, [f64; 3], f64) {
        let e = octanol.extent();
        let span = [0, 1, 2].map(|k| e[k][1] - e[k][0]);
        let (site, gap) = cell_for(span, octanol.mass() / density);
        let cell = PeriodicBox::new([0, 1, 2].map(|k| counts[k] as f64 * site[k]));
        (cell, site, gap)
    }

    fn build(
        solute: Option<Solute>,
        octanol: &Octanol,
        counts: [usize; 3],
        density: f64,
        seed: u64,
        cutoff: f64,
        accuracy: f64,
    ) -> Liquid {
        let (cell, site, gap) = Liquid::cell_of(octanol, counts, density);
        let e = octanol.extent();
        // The conformer's own middle at the site's centre, so that every gap is `g`.
        let middle = [0, 1, 2].map(|k| 0.5 * (e[k][0] + e[k][1]));
        let mut at: Vec<[f64; 3]> = Vec::new();
        let mut types = Vec::new();
        let mut charges = Vec::new();
        let mut masses = Vec::new();
        let mut molecules = Vec::new();
        let mut parts: Vec<&Component> = Vec::new();
        let n_solute = solute.as_ref().map_or(0, |s| s.at.len());
        if let Some(s) = &solute {
            parts.push(s.component);
            types.extend_from_slice(s.types);
            charges.extend_from_slice(s.charges);
            masses.extend(s.component.atoms().iter().map(|a| a.element.mass()));
            at.extend_from_slice(&s.at);
            molecules.push(0..s.at.len());
        }
        let mut octanols = 0;
        for i in 0..counts[0] {
            for j in 0..counts[1] {
                for k in 0..counts[2] {
                    let index = ((i * counts[1] + j) * counts[2] + k) as u64;
                    let mut rng = Rng::for_index(seed, index);
                    let end_for_end = rng.unit() < 0.5;
                    let about_its_axis = rng.unit() < 0.5;
                    let centre = [
                        (i as f64 + 0.5) * site[0],
                        (j as f64 + 0.5) * site[1],
                        (k as f64 + 0.5) * site[2],
                    ];
                    let placed: Vec<[f64; 3]> = octanol
                        .at
                        .iter()
                        .map(|p| {
                            let mut q = sub(*p, middle);
                            if end_for_end {
                                q = [-q[0], -q[1], q[2]];
                            }
                            if about_its_axis {
                                q = [q[0], -q[1], -q[2]];
                            }
                            [0, 1, 2].map(|a| centre[a] + q[a])
                        })
                        .collect();
                    if let Some(s) = &solute {
                        let near = placed.iter().any(|p| {
                            s.at.iter().any(|q| {
                                let d = cell.minimum_image(sub(*p, *q));
                                dot(d, d) < s.clearance * s.clearance
                            })
                        });
                        if near {
                            continue;
                        }
                    }
                    let start = at.len();
                    at.extend(placed);
                    molecules.push(start..at.len());
                    types.extend_from_slice(&octanol.types);
                    charges.extend_from_slice(&octanol.charges);
                    masses.extend_from_slice(&octanol.masses);
                    parts.push(&octanol.component);
                    octanols += 1;
                }
            }
        }
        let component = joined(&parts);
        let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
        let field = PeriodicForceField::new(&component, &types, cell, p)
            .expect("octanol in its box")
            .with_charges(charges.clone());
        Liquid {
            component,
            types,
            charges,
            field,
            at,
            masses,
            molecules,
            solute: n_solute,
            octanols,
            site,
            gap,
        }
    }

    pub fn cell(&self) -> PeriodicBox {
        self.field.cell()
    }

    /// The mass in the box over its volume, kg m⁻³.
    pub fn density(&self) -> f64 {
        self.masses.iter().sum::<f64>() / self.cell().volume()
    }

    /// Every bond to a hydrogen of every molecule held at UFF's natural length.
    pub fn shake(&self) -> Shake {
        Shake::to_hydrogen(
            &self.component,
            self.field.bonded().expect("the octanols are the solute"),
        )
    }

    /// Every X–H held, at rest, in no bath.
    pub fn dynamics(&self) -> MolecularDynamics {
        MolecularDynamics::new(self.masses.clone()).with_bond_constraints(self.shake())
    }
}

/// **The closest contact by brute force over the 27 images**: every pair of atoms in different
/// molecules at each of the 27 shifts `n ∈ {−1, 0, 1}³` of the box, and every atom with every atom
/// of its own molecule — itself included — at the 26 that are not zero, the positions first
/// wrapped into the box so that the 27 hold every nearest image. Returns the distance (metres) and
/// the two atoms.
pub fn closest_contact(
    at: &[[f64; 3]],
    molecules: &[Range<usize>],
    cell: &PeriodicBox,
) -> (f64, usize, usize) {
    let l = cell.lengths();
    let w: Vec<[f64; 3]> = at.iter().map(|p| cell.wrap(*p)).collect();
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
                        let r2 = dot(e, e);
                        if r2 < best.0 {
                            best = (r2, i, j);
                        }
                    }
                }
            }
        }
    }
    (best.0.sqrt(), best.1, best.2)
}

/// How many bonds of `c` have a hydrogen at one end, counted from the elements.
pub fn bonds_to_hydrogen(c: &Component) -> Vec<[usize; 2]> {
    c.bonds()
        .iter()
        .filter(|b| b.atoms.iter().any(|&i| c.atoms()[i].element == Element::H))
        .map(|b| b.atoms)
        .collect()
}
