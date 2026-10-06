//! What the constrained-bond and hydration tests share: benzene from the dictionary, typed, its
//! QEq charges, relaxed in vacuum and put on its constraints; and benzene in a box of TIP3P water
//! with every water rigid and every C–H bond held.

#![allow(dead_code)]

use pantometry_forcefield::uff::{self, UffType};
use pantometry_forcefield::water::{self, Settle, WaterBox};
use pantometry_forcefield::{
    qeq, Component, EwaldParameters, ForceField, MolecularDynamics, PeriodicBox,
    PeriodicDecoupling, PeriodicForceField, Shake,
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
