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
//! - **4c-1's measurement of the liquid** ([`measure_liquid`]): the melt ([`PLAN`]), the pressure
//!   by the molecular virial, `ΔH_vap` against a gas sampled under a potential it is given, and
//!   `g_OO` — shared by `liquid_octanol_in_a_box.rs`'s UFF models and
//!   `liquid_octanol_with_gaff.rs`'s GAFF one, so that each is the same protocol to the step.

#![allow(dead_code)]
// The geometry is written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

use pantometry_core::Rng;
use pantometry_forcefield::uff::{self, UffType, KCAL_PER_MOL};
use pantometry_forcefield::{
    qeq, Bath, Component, Element, Estimate, EwaldParameters, ForceField, MolecularDynamics,
    PeriodicBox, PeriodicForceField, PmeParameters, Potential, Shake,
};
use pantometry_units::BOLTZMANN;
use std::f64::consts::PI;
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

/// The bath's temperature, kelvin: 25 °C.
pub const KELVIN: f64 = 298.15;

/// The measurement's lattice: 192 molecules, 41.6 × 36.4 × 33.6 Å.
pub const LIQUID: [usize; 3] = [3, 8, 8];

/// One kJ/mol in joules per molecule: the crate's `KCAL_PER_MOL` over the thermochemical 4.184.
pub const KJ_PER_MOL: f64 = KCAL_PER_MOL / 4.184;

/// Joules per molecule in kcal/mol.
pub fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

/// `x` to the whole power `n` by repeated multiplication: `powi` is not used for any input or
/// reference, its precision being unspecified and platform-dependent.
pub fn pow(x: f64, n: u32) -> f64 {
    (0..n).fold(1.0, |p, _| p * x)
}

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

fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
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
        let ideal: Vec<[f64; 3]> = component.atoms().iter().map(|a| a.at).collect();
        let charges =
            charges.unwrap_or_else(|| qeq::charges(&component, &ideal).expect("QEq").charges);
        Octanol::relaxed(component, charges, ["CAA", "OAB", "CAC"])
    }

    /// **Any molecule made the same way**, `component` with `charges`: relaxed in vacuum under
    /// UFF from the dictionary's ideal geometry, every X–H put on its length, and turned into its
    /// frame by the atoms `frame = [a, b, c]` — `a → b` along `+x`, `c` in the `xy` plane at `+y`.
    /// Octanol's is `["CAA", "OAB", "CAC"]`; the struct keeps octanol's name because octanol was
    /// the first, and `tests/hexadecane/` builds its alkanes with it.
    pub fn relaxed(component: Component, charges: Vec<f64>, frame: [&str; 3]) -> Octanol {
        let types = uff::assign(&component);
        let ideal: Vec<[f64; 3]> = component.atoms().iter().map(|a| a.at).collect();
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
            index(&component, frame[0]),
            index(&component, frame[2]),
            index(&component, frame[1]),
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

// ---------------------------------------------------------------------------------------------
// 4c-1's measurement of the liquid

/// `steps` of constrained BAOAB at `dt` in a 5 ps⁻¹ bath at 298.15 K, from rest-thermalised
/// velocities, on `potential`; the dynamics and the positions.
pub fn bathed(
    template: &MolecularDynamics,
    potential: &impl Potential,
    start: &[[f64; 3]],
    dt: f64,
    steps: usize,
    seed: u64,
) -> (MolecularDynamics, Vec<[f64; 3]>) {
    let mut at = start.to_vec();
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 5e12,
            seed,
        })
        .thermalised(&at, KELVIN, seed);
    md.run(potential, &mut at, dt, steps);
    (md, at)
}

/// The pressure of the untruncated fluid at `at`, pascals, by the molecular virial: each
/// molecule's forces taken at its centre of mass, so that the constraint forces, internal to it,
/// drop out; the kinetic part the bath's `N_mol k_B T` exactly; and A&T's tail in place of the
/// virial's. W2's, with flexible molecules: the molecular virial is the atomic one less
/// `Σ_i (r_i − R_mol)·F_i`, on whole molecules, which molecular dynamics keeps them.
pub fn pressure(field: &PeriodicForceField, l: &Liquid, at: &[[f64; 3]]) -> f64 {
    let ev = field.evaluate(at);
    let mut w = ev.virial[0][0] + ev.virial[1][1] + ev.virial[2][2];
    for r in &l.molecules {
        let m: f64 = l.masses[r.clone()].iter().sum();
        let com = [0, 1, 2].map(|k| r.clone().map(|i| l.masses[i] * at[i][k]).sum::<f64>() / m);
        for i in r.clone() {
            w -= dot(sub(at[i], com), ev.forces[i]);
        }
    }
    let v = field.cell().volume();
    let n = l.molecules.len() as f64;
    (n * BOLTZMANN.to_si() * KELVIN + w / 3.0) / v - field.dispersion_correction() / v
        + field.pressure_correction()
}

/// The O–O histogram of one frame, normalised to `g(r)`, bins of `width` to `bins`.
pub fn g_oo(
    cell: &PeriodicBox,
    at: &[[f64; 3]],
    oxygens: &[usize],
    width: f64,
    bins: usize,
) -> Vec<f64> {
    let n = oxygens.len();
    let mut h = vec![0.0; bins];
    for a in 0..n {
        for b in a + 1..n {
            let d = cell.minimum_image(sub(at[oxygens[a]], at[oxygens[b]]));
            let k = (len(d) / width) as usize;
            if k < bins {
                h[k] += 1.0;
            }
        }
    }
    let v = cell.volume();
    for (k, x) in h.iter_mut().enumerate() {
        let (r0, r1) = (k as f64 * width, (k + 1) as f64 * width);
        let shell = 4.0 / 3.0 * PI * (pow(r1, 3) - pow(r0, 3));
        *x /= (n * (n - 1)) as f64 / 2.0 * shell / v;
    }
    h
}

/// The measurement's field: PME at δ = 10⁻⁵ by `PmeParameters::for_accuracy`, and the neighbour
/// list at its measured skin.
pub fn production_field(l: &Liquid) -> PeriodicForceField {
    production(l.field.clone(), l)
}

/// `field`, on [`Liquid`] `l`'s atoms, made the measurement's: the mesh and the list as
/// [`production_field`] gives them.
pub fn production(field: PeriodicForceField, l: &Liquid) -> PeriodicForceField {
    let p = field.ewald().parameters();
    let mesh = PmeParameters::for_accuracy(&l.cell(), &p, 1e-5, l.at.len());
    field
        .with_mesh(mesh)
        .with_neighbour_list(PeriodicForceField::NEIGHBOUR_SKIN)
}

/// One record of the melt: the time (ps), `U/N` (kcal/mol), the pressure (bar), the temperature (K).
pub type Row = (f64, f64, f64, f64);

/// One stage of the melt: its length, step, bath and temperature.
pub struct Stage {
    pub name: &'static str,
    pub picoseconds: f64,
    pub dt: f64,
    pub friction: f64,
    pub temperature: f64,
}

/// **The melt and the measurement**: 0.4 ps at 0.5 fs in a 50 ps⁻¹ bath at 298.15 K, which takes
/// out the lattice's strain; then 10 ps at 500 K and 2 fs to lose the lattice's order; then 14 ps
/// at 298.15 K to settle, then 26 ps sampled in a 1 ps⁻¹ bath: 25 000 steps of 2 fs, chosen to fit
/// an hour at the 115 ms a step [`the_cost_measured`] first read. It read 46.5 ms the second time,
/// and the runs took 48–64 ms a step (26 minutes) and 48–116 (44 minutes): the same code, the
/// machine's speed changing under it. **They are short for this liquid**, whose viscosity is 7.3 mPa s at 25 °C (CRC Handbook p. 6-234, as PubChem quotes it),
/// eight times water's, and the melt's blocks are printed so that what has not settled shows.
pub const PLAN: [Stage; 4] = [
    Stage {
        name: "strain out",
        picoseconds: 0.4,
        dt: 0.5 * FS,
        friction: 50e12,
        temperature: KELVIN,
    },
    Stage {
        name: "hot",
        picoseconds: 10.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: 500.0,
    },
    Stage {
        name: "settle",
        picoseconds: 14.0,
        dt: 2.0 * FS,
        friction: 5e12,
        temperature: KELVIN,
    },
    Stage {
        name: "sample",
        picoseconds: 26.0,
        dt: 2.0 * FS,
        friction: 1e12,
        temperature: KELVIN,
    },
];

/// `U` of one molecule in vacuum at 298.15 K under `potential`, joules: eight runs of `ps_each`
/// after 10 ps each, every X–H held, a 5 ps⁻¹ bath; the mean over the runs and its standard error
/// from their spread.
pub fn gas(o: &Octanol, potential: &impl Potential, ps_each: f64) -> (f64, f64) {
    let means: Vec<f64> = (0..8u64)
        .map(|s| {
            let template =
                MolecularDynamics::new(o.masses.clone()).with_bond_constraints(o.shake());
            let (mut md, mut at) = bathed(&template, potential, &o.at, 2.0 * FS, 5000, 0x6A5 + s);
            let steps = (ps_each * 1e-12 / (2.0 * FS)).round() as usize;
            let mut sum = 0.0;
            let mut count = 0usize;
            for k in 0..steps {
                md.step(potential, &mut at, 2.0 * FS);
                if k % 50 == 0 {
                    sum += md.potential_energy().expect("stepped");
                    count += 1;
                }
            }
            sum / count as f64
        })
        .collect();
    let m = means.iter().sum::<f64>() / 8.0;
    let var = means.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / 7.0;
    (m, (var / 8.0).sqrt())
}

/// **4c-1's measurement of liquid octanol `o`**, the field's non-bonded terms given by `rules` and
/// the gas sampled under `gas_potential`: [`LIQUID`]'s 192 molecules at [`DENSITY`], `r_c` = 9 Å,
/// PME at δ = 10⁻⁵ with the neighbour list, every X–H held, 2 fs; melted by [`PLAN`]. Reported:
/// every 0.1 ps of the melt, `U/N`, the pressure and the temperature in 2.5 ps blocks; over the
/// sample, `U/N`, `P` with their autocorrelation-aware errors and the sample's halves,
/// `ΔH_vap = U_gas − U_liq/N + RT`, and `g_OO`'s first peak and minimum. Nothing about the model is
/// asserted: only that the box has no overlap, every held bond stays held, and the run is finite.
pub fn measure_liquid(
    name: &str,
    o: &Octanol,
    rules: &dyn Fn(PeriodicForceField) -> PeriodicForceField,
    gas_potential: &impl Potential,
) {
    let started = std::time::Instant::now();
    let l = Liquid::lattice(o, LIQUID, DENSITY, 0x0C71, 9.0 * ANGSTROM, 1e-5);
    let (r, _, _) = closest_contact(&l.at, &l.molecules, &l.cell());
    println!(
        "\n{name}: {} octanols, {} atoms, box {:.3} × {:.3} × {:.3} Å at {:.6} g/cm³; closest contact {:.3} Å",
        l.octanols,
        l.at.len(),
        l.cell().lengths()[0] / ANGSTROM,
        l.cell().lengths()[1] / ANGSTROM,
        l.cell().lengths()[2] / ANGSTROM,
        l.density() / 1e3,
        r / ANGSTROM
    );
    assert!(r >= OVERLAP);
    let field = production(rules(l.field.clone()), &l);
    let n = l.octanols as f64;
    let oxygens: Vec<usize> = l
        .molecules
        .iter()
        .map(|m| m.start + index(&o.component, "OAB"))
        .collect();
    let mut at = l.at.clone();
    let mut md = l.dynamics();
    let mut series: Vec<Row> = Vec::new();
    let mut sample_u = Vec::new();
    let mut sample_p = Vec::new();
    let width = 0.02 * ANGSTROM;
    let bins = (l.cell().half_shortest_edge() / width) as usize;
    let mut g = vec![0.0; bins];
    let mut frames = 0usize;
    let mut time = 0.0;
    for (k, stage) in PLAN.iter().enumerate() {
        md = md.with_bath(Bath::Langevin {
            temperature: stage.temperature,
            friction: stage.friction,
            seed: 0x0C10 + k as u64,
        });
        if k == 0 {
            md = md.thermalised(&at, stage.temperature, 0x0C10);
        }
        md.forget_forces();
        let steps = (stage.picoseconds * 1e-12 / stage.dt).round() as usize;
        let every = ((0.1e-12 / stage.dt).round() as usize).max(1);
        let t0 = std::time::Instant::now();
        for s in 0..steps {
            md.step(&field, &mut at, stage.dt);
            if (s + 1) % every == 0 {
                time += every as f64 * stage.dt;
                let u = kcal(md.potential_energy().expect("stepped")) / n;
                let p = pressure(&field, &l, &at) / 1e5;
                series.push((time * 1e12, u, p, md.half_step_temperature()));
                if stage.name == "sample" {
                    sample_u.push(u);
                    sample_p.push(p);
                    if (s + 1) % (5 * every) == 0 {
                        for (x, y) in g
                            .iter_mut()
                            .zip(g_oo(&l.cell(), &at, &oxygens, width, bins))
                        {
                            *x += y;
                        }
                        frames += 1;
                    }
                }
            }
        }
        let elapsed = t0.elapsed().as_secs_f64();
        println!(
            "stage {} ({}): {} ps at {} fs, {} K, γ {} /ps: {:.0} s, {:.1} ms a step",
            k,
            stage.name,
            stage.picoseconds,
            stage.dt / FS,
            stage.temperature,
            stage.friction / 1e12,
            elapsed,
            1e3 * elapsed / steps as f64
        );
    }
    assert!(series.iter().all(|s| s.1.is_finite() && s.2.is_finite()));
    let shake = l.shake();
    let worst = shake
        .bonds()
        .iter()
        .zip(shake.lengths())
        .map(|(&[i, j], d)| (len(sub(at[i], at[j])) - d).abs() / d)
        .fold(0.0f64, f64::max);
    println!("the worst held bond is {worst:.2e} of its length off it");
    assert!(worst < 1e-10);

    println!("\nthe melt, in blocks of 2.5 ps:\n| t (ps) | U/N (kcal/mol) | P (bar) | T (K) |\n| --- | --- | --- | --- |");
    for block in series.chunks(25) {
        let m = |f: &dyn Fn(&Row) -> f64| block.iter().map(f).sum::<f64>() / block.len() as f64;
        println!(
            "| {:.1}–{:.1} | {:.3} | {:.0} | {:.1} |",
            block[0].0 - 0.1,
            block[block.len() - 1].0,
            m(&|s| s.1),
            m(&|s| s.2),
            m(&|s| s.3)
        );
    }
    let u = Estimate::of(&sample_u);
    let p = Estimate::of(&sample_p);
    // The sample's halves, as a check on drift the autocorrelation cannot see.
    let half = sample_u.len() / 2;
    let (u1, u2) = (
        Estimate::of(&sample_u[..half]),
        Estimate::of(&sample_u[half..]),
    );
    let (p1, p2) = (
        Estimate::of(&sample_p[..half]),
        Estimate::of(&sample_p[half..]),
    );
    println!(
        "\nU/N = {:.3} ± {:.3} kcal/mol (τ {:.1} samples of 0.1 ps); halves {:.3} ± {:.3}, {:.3} ± {:.3}",
        u.mean, u.error, u.tau, u1.mean, u1.error, u2.mean, u2.error
    );
    println!(
        "P = {:.0} ± {:.0} bar (τ {:.1}); halves {:.0} ± {:.0}, {:.0} ± {:.0}; experiment at this density: 1 bar",
        p.mean, p.error, p.tau, p1.mean, p1.error, p2.mean, p2.error
    );
    let (ug, ug_err) = gas(o, gas_potential, 200.0);
    let rt = BOLTZMANN.to_si() * KELVIN;
    let u_liq = u.mean * KCAL_PER_MOL;
    let dh = ug - u_liq + rt;
    let dh_err = (ug_err * ug_err + (u.error * KCAL_PER_MOL) * (u.error * KCAL_PER_MOL)).sqrt();
    println!(
        "U_gas = {:.3} ± {:.3} kcal/mol (eight runs of 200 ps); ΔH_vap = U_gas − U/N + RT = {:.2} \
         ± {:.2} kcal/mol = {:.1} ± {:.1} kJ/mol",
        kcal(ug),
        kcal(ug_err),
        kcal(dh),
        kcal(dh_err),
        dh / KJ_PER_MOL,
        dh_err / KJ_PER_MOL
    );
    println!("  experiment: 70.98 kJ/mol at 25 °C (CRC Handbook 95th ed. p. 6-142, as PubChem CID 957 quotes it, read)");
    let mean_g: Vec<f64> = g.iter().map(|x| x / frames as f64).collect();
    let lo = (2.3 * ANGSTROM / width) as usize;
    let hi = (3.4 * ANGSTROM / width) as usize;
    let b = (lo..hi)
        .max_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    let (y0, y1, y2) = (mean_g[b - 1], mean_g[b], mean_g[b + 1]);
    let shift = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2);
    let r_peak = (b as f64 + 0.5 + shift) * width;
    let min_hi = (4.5 * ANGSTROM / width) as usize;
    let m = (b..min_hi)
        .min_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    let second = (m..((6.5 * ANGSTROM / width) as usize).min(bins))
        .max_by(|&i, &j| mean_g[i].total_cmp(&mean_g[j]))
        .expect("bins");
    // The coordination number to the first minimum: the oxygens a hydroxyl's neighbours are.
    let rho = n / l.cell().volume();
    let coordination: f64 = (0..m)
        .map(|k| {
            let (r0, r1) = (k as f64 * width, (k + 1) as f64 * width);
            mean_g[k] * rho * 4.0 / 3.0 * PI * (pow(r1, 3) - pow(r0, 3))
        })
        .sum();
    println!(
        "g_OO over {frames} frames: first peak at {:.3} Å, height {:.2}; first minimum at {:.2} Å, \
         {:.3}; {:.2} oxygens inside it; the next maximum at {:.2} Å, {:.2}",
        r_peak / ANGSTROM,
        y1 - 0.25 * (y0 - y2) * shift,
        (m as f64 + 0.5) * width / ANGSTROM,
        mean_g[m],
        coordination,
        (second as f64 + 0.5) * width / ANGSTROM,
        mean_g[second]
    );
    println!("the whole run {:.0} s", started.elapsed().as_secs_f64());
}
