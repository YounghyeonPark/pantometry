//! What the constrained-bond and hydration tests share: benzene from the dictionary, typed, its
//! QEq charges, relaxed in vacuum and put on its constraints; and benzene in a box of TIP3P water
//! with every water rigid and every C–H bond held.

#![allow(dead_code)]

use pantometry_forcefield::uff::{self, UffType};
use pantometry_forcefield::water::{self, Settle, WaterBox};
use pantometry_forcefield::{
    qeq, AtLambda, Bath, Component, EwaldParameters, ForceField, Lambda, MolecularDynamics,
    PeriodicBox, PeriodicDecoupling, PeriodicForceField, Shake,
};

pub const ANGSTROM: f64 = 1e-10;
pub const FS: f64 = 1e-15;

/// Benzene, `BNZ` of the Chemical Component Dictionary, byte for byte.
pub const BNZ: &str = include_str!("../../components/BNZ.cif");

/// Benzene alone: its component, types, charges, vacuum force field and a relaxed geometry on
/// its constraints.
pub struct Benzene {
    pub component: Component,
    pub types: Vec<UffType>,
    pub charges: Vec<f64>,
    pub vacuum: ForceField,
    /// Relaxed in vacuum, then every C–H set to UFF's natural length along itself.
    pub at: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
}

impl Benzene {
    /// Benzene with `charges`, or QEq's at the dictionary's ideal geometry when `None`.
    pub fn new(charges: Option<Vec<f64>>) -> Benzene {
        let component = Component::from_ccd(BNZ).expect("BNZ");
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
        // Before and after the same: each bond moved along itself onto its length.
        let before = at.clone();
        shake.constrain_positions(&masses, &vec![false; at.len()], &before, &mut at);
        Benzene {
            component,
            types,
            charges,
            vacuum,
            at,
            masses,
        }
    }

    /// The six C–H constraints, on benzene's own atom indices.
    pub fn shake(&self) -> Shake {
        Shake::to_hydrogen(&self.component, &self.vacuum)
    }

    /// The positions moved so that their centroid is at `centre`.
    pub fn centred_at(&self, centre: [f64; 3]) -> Vec<[f64; 3]> {
        let n = self.at.len() as f64;
        let c = [0, 1, 2].map(|k| self.at.iter().map(|p| p[k]).sum::<f64>() / n);
        self.at
            .iter()
            .map(|p| [0, 1, 2].map(|k| p[k] - c[k] + centre[k]))
            .collect()
    }
}

/// Benzene in a box of rigid TIP3P water: the force field, the positions (benzene's twelve first),
/// the masses and the dynamics template that holds every water by SETTLE and every C–H by SHAKE.
pub struct Solvated {
    pub field: PeriodicForceField,
    pub at: Vec<[f64; 3]>,
    pub masses: Vec<f64>,
    pub waters: usize,
    /// Benzene's atoms.
    pub group: Vec<bool>,
}

impl Solvated {
    /// `per_side³` lattice waters at 0.997 g/cm³ turned by `seed`, benzene at the centre, every
    /// water with an atom within `clearance` of a benzene atom taken out; Ewald at `cutoff` and
    /// `accuracy`.
    pub fn new(
        benzene: &Benzene,
        per_side: usize,
        seed: u64,
        clearance: f64,
        cutoff: f64,
        accuracy: f64,
    ) -> Solvated {
        Solvated::at_density(
            benzene,
            per_side,
            water::LIQUID_DENSITY,
            seed,
            clearance,
            cutoff,
            accuracy,
        )
    }

    /// The same with the lattice at `density` kg m⁻³.
    pub fn at_density(
        benzene: &Benzene,
        per_side: usize,
        density: f64,
        seed: u64,
        clearance: f64,
        cutoff: f64,
        accuracy: f64,
    ) -> Solvated {
        let lattice = WaterBox::lattice_at(per_side, density, seed);
        let cell = lattice.cell();
        let side = cell.lengths()[0];
        let solute = benzene.centred_at([0.5 * side; 3]);
        let wb = lattice.without_overlaps(&solute, clearance);
        let p = EwaldParameters::for_accuracy(&cell, cutoff, accuracy);
        let field =
            PeriodicForceField::solvated(&benzene.component, &benzene.types, wb.count(), cell, p)
                .expect("benzene in water")
                .with_solute_charges(benzene.charges.clone());
        let mut at = solute;
        at.extend_from_slice(wb.positions());
        let mut masses = benzene.masses.clone();
        masses.extend(wb.masses());
        let mut group = vec![false; at.len()];
        group[..benzene.at.len()].fill(true);
        Solvated {
            field,
            at,
            masses,
            waters: wb.count(),
            group,
        }
    }

    pub fn cell(&self) -> PeriodicBox {
        self.field.cell()
    }

    /// Every water rigid and every C–H held, at rest.
    pub fn dynamics(&self, benzene: &Benzene) -> MolecularDynamics {
        MolecularDynamics::new(self.masses.clone())
            .with_constraints(Settle::tip3p(self.field.rigid_waters().to_vec()))
            .with_bond_constraints(benzene.shake())
    }

    /// Benzene decoupled from the water.
    pub fn decoupling(&self) -> PeriodicDecoupling {
        PeriodicDecoupling::new(&self.field, &self.group).expect("benzene decouples")
    }

    /// The water alone in the same box, with the same Ewald parameters.
    pub fn water_alone(&self) -> PeriodicForceField {
        PeriodicForceField::tip3p(self.cell(), self.field.ewald().parameters(), self.waters)
    }
}

/// The TIP3P masses, oxygen first: the crate's.
pub fn water_masses() -> [f64; 3] {
    water::masses()
}

/// W3's lattice seed: the measurement's box is turned by it.
pub const W3_SEED: u64 = 0x3_3A7E;

/// The density W3's box is built at: 33.00 waters per nm³, 0.9872 g/cm³, where Yeh and Hummer
/// (*J. Phys. Chem. B* **108**, 15873 (2004), read in W2) measure TIP3P's pressure under Ewald at
/// −2.9 ± 2.5 bar. W2's liquid at 0.997 g/cm³ read +139 to +180 bar; a hydration free energy at
/// that pressure would carry the solute's partial molar volume times it, a few tenths of a kcal/mol,
/// so the box is at the model's own density at one atmosphere instead of experiment's.
pub fn w3_density() -> f64 {
    33.0e27 * water::molecular_mass()
}

/// W3's box: 512 lattice waters at [`w3_density`] (24.94 Å, the smallest lattice box past `2 r_c`
/// plus benzene's 5.0 Å span at `r_c` = 9 Å), less every one with an atom within 1.4 Å of a benzene
/// atom — five, about benzene's volume in waters (liquid benzene's 89 cm³/mol is 4.9 of water's
/// 18.07) — and Ewald at δ = 10⁻⁵. Shared by W3's measurement and W4's decomposition, so that the
/// two are the same box.
pub fn w3_solvated(benzene: &Benzene) -> Solvated {
    Solvated::at_density(
        benzene,
        8,
        w3_density(),
        W3_SEED,
        1.4 * ANGSTROM,
        9.0 * ANGSTROM,
        1e-5,
    )
}

/// W3's melt of its box, in place, at `temperature`: 0.5 ps at 0.5 fs in a 50 ps⁻¹ bath at each of
/// `λ_v` = 0.5 with no charges — the soft core finite where a lattice water overlaps benzene —
/// `λ_v` = 1 with no charges, and fully coupled; then 20 ps at 2 fs in a 5 ps⁻¹ bath, fully
/// coupled. W2's melt with the soft core in front of it; every constraint of `template` on
/// throughout.
pub fn w3_melt(
    decoupling: &PeriodicDecoupling,
    template: &MolecularDynamics,
    temperature: f64,
    at: &mut [[f64; 3]],
) {
    let potential = AtLambda {
        hamiltonian: decoupling,
        lambda: Lambda::COUPLED,
    };
    let bath = |friction: f64, seed: u64| Bath::Langevin {
        temperature,
        friction,
        seed,
    };
    let mut md = template
        .clone()
        .with_bath(bath(50e12, 0xE01))
        .thermalised(at, temperature, 0xE01);
    for (k, lambda) in [Lambda::new(0.0, 0.0, 0.5), Lambda::new(0.0, 0.0, 1.0)]
        .into_iter()
        .enumerate()
    {
        let growing = AtLambda {
            hamiltonian: decoupling,
            lambda,
        };
        md = md.with_bath(bath(50e12, 0xE10 + k as u64));
        md.forget_forces();
        md.run(&growing, at, 0.5 * FS, 1000);
    }
    md.forget_forces();
    md.run(&potential, at, 0.5 * FS, 1000);
    let mut md = md.with_bath(bath(5e12, 0xE02));
    md.run(&potential, at, 2.0 * FS, 10_000);
}
