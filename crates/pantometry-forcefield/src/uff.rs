//! Universal Force Field atom types, and the parameters the paper gives for them.
//!
//! The Universal Force Field (UFF) is one parameter set for the whole periodic table, built from
//! per-element rules rather than fitted per molecule:
//!
//! > A. K. Rappé, C. J. Casewit, K. S. Colwell, W. A. Goddard III and W. M. Skiff,
//! > "UFF, a full periodic table force field for molecular mechanics and molecular dynamics
//! > simulations", *J. Am. Chem. Soc.* **114**, 10024–10035 (1992).
//! > [doi:10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040)
//!
//! An atom's UFF type is its element and its hybridisation — `C_3` is sp³ carbon, `C_R` resonant
//! (aromatic) carbon, `C_2` sp², `C_1` sp — and every energy term is built from the six numbers
//! Table I gives per type. This module holds those numbers for the twenty-two types the ten
//! elements in [`Element`] need, and the rules ([`assign`]) that read a type off a
//! [`Component`]'s bonds — every one from the atom's own bonds, except the resonant oxygen or
//! sulfur, which also reads its neighbours' (an ether, ester or thioether O or S on an sp² atom is
//! `O_R` or `S_R`, as the paper types methyl vinyl ether's oxygen).
//!
//! # Table I, as printed
//!
//! The paper's units. [`UffType::parameters`] converts to SI, which is what everything past this
//! table uses; [`UffType::table_i`] returns the row as printed, so a reader can compare the two
//! columns of numbers against the paper without doing the conversion in their head.
//!
//! | type | r₁ (Å) | θ₀ (°) | x₁ (Å) | D₁ (kcal/mol) | ζ | Z₁ (e) |
//! | --- | --- | --- | --- | --- | --- | --- |
//! | `H_` | 0.354 | 180.0 | 2.886 | 0.044 | 12.0 | 0.712 |
//! | `C_3` | 0.757 | 109.47 | 3.851 | 0.105 | 12.73 | 1.912 |
//! | `C_R` | 0.729 | 120.0 | 3.851 | 0.105 | 12.73 | 1.912 |
//! | `C_2` | 0.732 | 120.0 | 3.851 | 0.105 | 12.73 | 1.912 |
//! | `C_1` | 0.706 | 180.0 | 3.851 | 0.105 | 12.73 | 1.912 |
//! | `N_3` | 0.700 | 106.7 | 3.660 | 0.069 | 13.407 | 2.544 |
//! | `N_R` | 0.699 | 120.0 | 3.660 | 0.069 | 13.407 | 2.544 |
//! | `N_2` | 0.685 | 111.2 | 3.660 | 0.069 | 13.407 | 2.544 |
//! | `N_1` | 0.656 | 180.0 | 3.660 | 0.069 | 13.407 | 2.544 |
//! | `O_3` | 0.658 | 104.51 | 3.500 | 0.060 | 14.085 | 2.300 |
//! | `O_R` | 0.680 | 110.0 | 3.500 | 0.060 | 14.085 | 2.300 |
//! | `O_2` | 0.634 | 120.0 | 3.500 | 0.060 | 14.085 | 2.300 |
//! | `O_1` | 0.639 | 180.0 | 3.500 | 0.060 | 14.085 | 2.300 |
//! | `F_` | 0.668 | 180.0 | 3.364 | 0.050 | 14.762 | 1.735 |
//! | `P_3+3` | 1.101 | 93.8 | 4.147 | 0.305 | 13.072 | 2.863 |
//! | `P_3+5` | 1.056 | 109.47 | 4.147 | 0.305 | 13.072 | 2.863 |
//! | `S_3+2` | 1.064 | 92.1 | 4.035 | 0.274 | 13.969 | 2.703 |
//! | `S_R` | 1.077 | 92.2 | 4.035 | 0.274 | 13.969 | 2.703 |
//! | `S_2` | 0.854 | 120.0 | 4.035 | 0.274 | 13.969 | 2.703 |
//! | `Cl` | 1.044 | 180.0 | 3.947 | 0.227 | 14.866 | 2.348 |
//! | `Br` | 1.192 | 180.0 | 4.189 | 0.251 | 15.0 | 2.519 |
//! | `I_` | 1.382 | 180.0 | 4.50 | 0.339 | 15.0 | 2.65 |
//!
//! r₁ is the bond radius, θ₀ the natural angle, x₁ the van der Waals distance, D₁ the van der
//! Waals well depth, ζ the van der Waals shape (scale) parameter, and Z₁ the effective charge used
//! by the bond-stretch and angle-bend force constants. Transcribed from the scanned original and
//! checked against Open Babel's `data/UFF.prm`, which carries the same values.
//!
//! # The GMP electronegativity χ, and where it comes from
//!
//! UFF's natural bond length carries an electronegativity correction `r_EN` (eq 4, p. 10027)
//! built from the generalized Mulliken–Pauling electronegativity χ of each element.
//! [`gmp_electronegativity`] gives it, in eV, for the ten elements here:
//!
//! | H | C | N | O | F | P | S | Cl | Br | I |
//! | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
//! | 4.528 | 5.343 | 6.899 | 8.741 | 10.874 | 5.463 | 6.928 | 8.564 | 7.790 | 6.822 |
//!
//! **Provenance.** These numbers were first transcribed from Open Babel's `data/UFF.prm`, which
//! says it copies RDKit's, while the paper they come from had not been read. It has now: the UFF
//! paper takes χ from the charge-equilibration paper — A. K. Rappé and W. A. Goddard III,
//! *J. Phys. Chem.* **95**, 3358 (1991) — whose **Table I (p. 3359) prints all ten, equal to these
//! exactly**, citing the GMP paper (its ref 9). Hydrogen's 4.528 is Table I's 4.5280, fitted in
//! that paper to experimental charges (eq 22, p. 3361). The test
//! `uffs_electronegativities_are_table_i` holds the two tables to each other: [`crate::qeq`] reads
//! the same χ from its own transcription of Table I. The partial check made before the paper was
//! read still stands: the UFF paper's worked example (p. 10027) gives the Si–`O_3_z` correction as
//! 0.0533 Å, and eq 4 with χ_Si 4.168 (Table I's) and χ_O 8.741 and Table I's r_Si 1.117 and
//! r_O_3_z 0.528 gives 0.05325 Å. The test `the_papers_silicon_oxygen_example` carries it.
//!
//! # Torsion barriers
//!
//! The sp³ barrier `V` is the paper's Table III (p. 10028) — `C_3` 2.119, `N_3` 0.450, `O_3`
//! 0.018, `P_3` 2.400, `S_3` 0.484 kcal/mol — and the sp² constant `U` is not in a table but in
//! the text under its eq 17 on the same page, assigned by period: 2 for the second row, 1.25 for
//! the third, 0.7 for the fourth, 0.2 for the fifth. Both read against the scanned paper on
//! 2026-10-01, and both are what Open Babel carries as `Vi` and `Uj`. See
//! [`UffType::sp3_torsion_barrier`] and [`UffType::sp2_torsion_barrier`].
//!
//! # Units
//!
//! Lengths in metres, angles in radians, energies in **joules per molecule** (kcal/mol divided
//! by Avogadro's number, as `pantometry-molecular` stores a Lennard-Jones well), and charges in
//! coulombs.

use crate::ccd::{BondOrder, Component, Element};
use std::fmt;

/// Avogadro's number, mol⁻¹ (exact by definition since 2019).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// One thermochemical kilocalorie per mole, in joules per molecule: `4184 J / N_A`.
pub const KCAL_PER_MOL: f64 = 4184.0 / AVOGADRO;

/// The elementary charge, in coulombs (exact by definition since 2019).
pub const ELEMENTARY_CHARGE: f64 = 1.602_176_634e-19;

/// One degree, in radians.
const DEGREE: f64 = std::f64::consts::PI / 180.0;

/// A Universal Force Field atom type: an element and a hybridisation.
///
/// Spelled in Rust case; [`UffType::label`] gives the paper's spelling (`C_3`, `P_3+5`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum UffType {
    /// `H_`, hydrogen.
    H,
    /// `C_3`, sp³ carbon.
    C3,
    /// `C_R`, resonant (aromatic) carbon.
    CR,
    /// `C_2`, sp² carbon.
    C2,
    /// `C_1`, sp carbon.
    C1,
    /// `N_3`, sp³ nitrogen.
    N3,
    /// `N_R`, resonant nitrogen: aromatic, or an amide's.
    NR,
    /// `N_2`, sp² nitrogen.
    N2,
    /// `N_1`, sp nitrogen.
    N1,
    /// `O_3`, sp³ oxygen.
    O3,
    /// `O_R`, aromatic oxygen.
    OR,
    /// `O_2`, sp² oxygen.
    O2,
    /// `O_1`, sp oxygen.
    O1,
    /// `F_`, fluorine.
    F,
    /// `P_3+3`, trivalent sp³ phosphorus.
    P3Trivalent,
    /// `P_3+5`, pentavalent sp³ phosphorus — a phosphate's.
    P3Pentavalent,
    /// `S_3+2`, divalent sp³ sulfur.
    S3Divalent,
    /// `S_R`, aromatic sulfur.
    SR,
    /// `S_2`, sp² sulfur.
    S2,
    /// `Cl`, chlorine.
    Cl,
    /// `Br`, bromine.
    Br,
    /// `I_`, iodine.
    I,
}

/// A hybridisation as the torsion rules use it. See [`UffType::hybridisation`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Hybridisation {
    /// Linear, `X_1`. A torsion about a bond to one is zero (p. 10029).
    Sp,
    /// Trigonal: `X_2`, and the resonant `X_R`.
    Sp2,
    /// Tetrahedral, `X_3`.
    Sp3,
}

/// One row of Table I in the paper's own units. See [`UffType::table_i`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TableI {
    /// r₁, bond radius, Å.
    pub r1: f64,
    /// θ₀, natural angle, degrees.
    pub theta0: f64,
    /// x₁, van der Waals distance, Å.
    pub x1: f64,
    /// D₁, van der Waals energy, kcal/mol.
    pub d1: f64,
    /// ζ, van der Waals scale, dimensionless.
    pub zeta: f64,
    /// Z₁, effective charge, elementary charges.
    pub z1: f64,
}

/// One type's parameters in SI. See [`UffType::parameters`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Parameters {
    /// Bond radius r₁, in metres.
    pub bond_radius: f64,
    /// Natural bond angle θ₀, in radians.
    pub angle: f64,
    /// Van der Waals distance x₁, in metres.
    pub vdw_distance: f64,
    /// Van der Waals well depth D₁, in joules per molecule.
    pub vdw_energy: f64,
    /// Van der Waals shape parameter ζ, dimensionless.
    pub vdw_scale: f64,
    /// Effective charge Z₁, in coulombs.
    pub effective_charge: f64,
}

impl UffType {
    /// Every type this crate has parameters for, in Table I's order.
    pub const ALL: [UffType; 22] = [
        UffType::H,
        UffType::C3,
        UffType::CR,
        UffType::C2,
        UffType::C1,
        UffType::N3,
        UffType::NR,
        UffType::N2,
        UffType::N1,
        UffType::O3,
        UffType::OR,
        UffType::O2,
        UffType::O1,
        UffType::F,
        UffType::P3Trivalent,
        UffType::P3Pentavalent,
        UffType::S3Divalent,
        UffType::SR,
        UffType::S2,
        UffType::Cl,
        UffType::Br,
        UffType::I,
    ];

    /// The paper's spelling: `H_`, `C_3`, `P_3+5`, `Cl`.
    pub fn label(self) -> &'static str {
        match self {
            UffType::H => "H_",
            UffType::C3 => "C_3",
            UffType::CR => "C_R",
            UffType::C2 => "C_2",
            UffType::C1 => "C_1",
            UffType::N3 => "N_3",
            UffType::NR => "N_R",
            UffType::N2 => "N_2",
            UffType::N1 => "N_1",
            UffType::O3 => "O_3",
            UffType::OR => "O_R",
            UffType::O2 => "O_2",
            UffType::O1 => "O_1",
            UffType::F => "F_",
            UffType::P3Trivalent => "P_3+3",
            UffType::P3Pentavalent => "P_3+5",
            UffType::S3Divalent => "S_3+2",
            UffType::SR => "S_R",
            UffType::S2 => "S_2",
            UffType::Cl => "Cl",
            UffType::Br => "Br",
            UffType::I => "I_",
        }
    }

    /// The element this type is a hybridisation of.
    pub fn element(self) -> Element {
        match self {
            UffType::H => Element::H,
            UffType::C3 | UffType::CR | UffType::C2 | UffType::C1 => Element::C,
            UffType::N3 | UffType::NR | UffType::N2 | UffType::N1 => Element::N,
            UffType::O3 | UffType::OR | UffType::O2 | UffType::O1 => Element::O,
            UffType::F => Element::F,
            UffType::P3Trivalent | UffType::P3Pentavalent => Element::P,
            UffType::S3Divalent | UffType::SR | UffType::S2 => Element::S,
            UffType::Cl => Element::Cl,
            UffType::Br => Element::Br,
            UffType::I => Element::I,
        }
    }

    /// This type's row of Table I, as printed — see the module documentation for the table.
    pub fn table_i(self) -> TableI {
        let (r1, theta0, x1, d1, zeta, z1) = match self {
            UffType::H => (0.354, 180.0, 2.886, 0.044, 12.0, 0.712),
            UffType::C3 => (0.757, 109.47, 3.851, 0.105, 12.73, 1.912),
            UffType::CR => (0.729, 120.0, 3.851, 0.105, 12.73, 1.912),
            UffType::C2 => (0.732, 120.0, 3.851, 0.105, 12.73, 1.912),
            UffType::C1 => (0.706, 180.0, 3.851, 0.105, 12.73, 1.912),
            UffType::N3 => (0.700, 106.7, 3.660, 0.069, 13.407, 2.544),
            UffType::NR => (0.699, 120.0, 3.660, 0.069, 13.407, 2.544),
            UffType::N2 => (0.685, 111.2, 3.660, 0.069, 13.407, 2.544),
            UffType::N1 => (0.656, 180.0, 3.660, 0.069, 13.407, 2.544),
            UffType::O3 => (0.658, 104.51, 3.500, 0.060, 14.085, 2.300),
            UffType::OR => (0.680, 110.0, 3.500, 0.060, 14.085, 2.300),
            UffType::O2 => (0.634, 120.0, 3.500, 0.060, 14.085, 2.300),
            UffType::O1 => (0.639, 180.0, 3.500, 0.060, 14.085, 2.300),
            UffType::F => (0.668, 180.0, 3.364, 0.050, 14.762, 1.735),
            UffType::P3Trivalent => (1.101, 93.8, 4.147, 0.305, 13.072, 2.863),
            UffType::P3Pentavalent => (1.056, 109.47, 4.147, 0.305, 13.072, 2.863),
            UffType::S3Divalent => (1.064, 92.1, 4.035, 0.274, 13.969, 2.703),
            UffType::SR => (1.077, 92.2, 4.035, 0.274, 13.969, 2.703),
            UffType::S2 => (0.854, 120.0, 4.035, 0.274, 13.969, 2.703),
            UffType::Cl => (1.044, 180.0, 3.947, 0.227, 14.866, 2.348),
            UffType::Br => (1.192, 180.0, 4.189, 0.251, 15.0, 2.519),
            UffType::I => (1.382, 180.0, 4.50, 0.339, 15.0, 2.65),
        };
        TableI {
            r1,
            theta0,
            x1,
            d1,
            zeta,
            z1,
        }
    }

    /// This type's Table I row in SI. See the module documentation for the units.
    pub fn parameters(self) -> Parameters {
        let t = self.table_i();
        Parameters {
            bond_radius: t.r1 * crate::ccd::ANGSTROM,
            angle: t.theta0 * DEGREE,
            vdw_distance: t.x1 * crate::ccd::ANGSTROM,
            vdw_energy: t.d1 * KCAL_PER_MOL,
            vdw_scale: t.zeta,
            effective_charge: t.z1 * ELEMENTARY_CHARGE,
        }
    }

    /// The sp³ torsion barrier `V`, in joules per molecule, for the sp³ types the paper gives one
    /// for: `C_3` 2.119, `N_3` 0.450, `O_3` 0.018, `P_3` 2.400, `S_3` 0.484 kcal/mol.
    ///
    /// The paper gives `P_3` without the oxidation-state suffix, so both `P_3+3` and `P_3+5`
    /// carry it — the reading Open Babel's `UFF.prm` makes too. `None` for every type that is not
    /// sp³ or has no tabulated value.
    pub fn sp3_torsion_barrier(self) -> Option<f64> {
        let kcal = match self {
            UffType::C3 => 2.119,
            UffType::N3 => 0.450,
            UffType::O3 => 0.018,
            UffType::P3Trivalent | UffType::P3Pentavalent => 2.400,
            UffType::S3Divalent => 0.484,
            _ => return None,
        };
        Some(kcal * KCAL_PER_MOL)
    }

    /// The hybridisation the torsion rules (p. 10028–29) read this type as, or `None` for a type
    /// that has none — hydrogen and the halogens, which are never the centre of a torsion.
    ///
    /// The paper's own suffix decides it: `_3` is sp³, `_1` sp, and both `_2` and `_R` are sp².
    /// **Resonant counts as sp²**, a reading and not a quotation: eq 17, the sp²–sp² barrier, was
    /// fitted to benzene's vibrational modes among others (p. 10028), which only makes sense if an
    /// aromatic `C_R` is an sp² centre; and the group-6 exception says "an sp² or resonant atom",
    /// which is the one place the paper names the two separately. So `C_R`, `N_R`, `O_R` and `S_R`
    /// are [`Hybridisation::Sp2`], as are the amide's `C_R` and `N_R`. `P_3+3`, `P_3+5` and
    /// `S_3+2` are sp³.
    pub fn hybridisation(self) -> Option<Hybridisation> {
        match self {
            UffType::C3
            | UffType::N3
            | UffType::O3
            | UffType::P3Trivalent
            | UffType::P3Pentavalent
            | UffType::S3Divalent => Some(Hybridisation::Sp3),
            UffType::CR
            | UffType::C2
            | UffType::NR
            | UffType::N2
            | UffType::OR
            | UffType::O2
            | UffType::SR
            | UffType::S2 => Some(Hybridisation::Sp2),
            UffType::C1 | UffType::N1 | UffType::O1 => Some(Hybridisation::Sp),
            UffType::H | UffType::F | UffType::Cl | UffType::Br | UffType::I => None,
        }
    }

    /// The sp² torsion constant `U`, in joules per molecule, by the element's period: 2, 1.25,
    /// 0.7 and 0.2 kcal/mol for periods two to five.
    ///
    /// Given for every type of a period-two-to-five element, because the paper tabulates it by
    /// row and not by hybridisation; which torsions it enters is the energy expression's business
    /// and not this table's. `None` for hydrogen, which has no row.
    pub fn sp2_torsion_barrier(self) -> Option<f64> {
        let kcal = match self.element().period() {
            2 => 2.0,
            3 => 1.25,
            4 => 0.7,
            5 => 0.2,
            _ => return None,
        };
        Some(kcal * KCAL_PER_MOL)
    }
}

/// The GMP electronegativity χ of `element`, in **electronvolts**.
///
/// Left in eV rather than converted, because the one place it enters — eq 4's `r_EN` — is
/// homogeneous of degree zero in χ: scaling every χ by one factor leaves `r_EN` unchanged, so the
/// unit cancels and a conversion would only add a rounding. See the module documentation for the
/// table and its source, the QEq paper's Table I.
pub fn gmp_electronegativity(element: Element) -> f64 {
    match element {
        Element::H => 4.528,
        Element::C => 5.343,
        Element::N => 6.899,
        Element::O => 8.741,
        Element::F => 10.874,
        Element::P => 5.463,
        Element::S => 6.928,
        Element::Cl => 8.564,
        Element::Br => 7.790,
        Element::I => 6.822,
    }
}

impl fmt::Display for UffType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.label())
    }
}

/// What one atom's bonds look like, which is all the typing rules read.
#[derive(Clone, Copy, Debug, Default)]
struct Bonding {
    aromatic: bool,
    doubles: usize,
    triples: usize,
    double_to_oxygen: bool,
}

fn bonding(component: &Component, i: usize) -> Bonding {
    let atoms = component.atoms();
    let mut b = Bonding {
        aromatic: atoms[i].aromatic,
        ..Bonding::default()
    };
    for (other, bond) in component.neighbours(i) {
        b.aromatic |= bond.aromatic;
        match bond.order {
            BondOrder::Single => {}
            BondOrder::Double => {
                b.doubles += 1;
                b.double_to_oxygen |= atoms[other].element == Element::O;
            }
            BondOrder::Triple => b.triples += 1,
        }
    }
    b
}

/// The UFF type of every atom of `component`, in its order.
///
/// Read from the element, the bond orders and the aromatic flags — **no ring perception**,
/// because the dictionary already states aromaticity per atom and per bond, and a second opinion
/// computed here could only disagree with it. An atom counts as aromatic when its own flag is set
/// *or* any of its bonds' is; the two agree in every well-formed entry, and taking either is a
/// choice made so that one missing flag does not silently demote a ring atom to sp².
///
/// The rules, each with its reason:
///
/// - **H → `H_`; F → `F_`; Cl → `Cl`; Br → `Br`; I → `I_`.** One type per element in UFF.
/// - **Carbon.** Aromatic → `C_R`. Otherwise any triple bond, or two double bonds (an allene's
///   centre, CO₂), → `C_1`: both are linear, two-coordinate sp carbon. Otherwise any double bond
///   → `C_2`. Otherwise `C_3`. The aromatic test comes first because the CCD writes a ring in
///   Kekulé form, so a ring carbon *has* a double bond and would read as `C_2` without it.
///   **An amide carbon is `C_R`**: a carbon with a double bond to oxygen and a single,
///   non-aromatic bond to a nitrogen the amide rule below types `N_R` (non-aromatic, single bonds
///   only). The paper's amide bond order on p. 10026 is derived "from the C_R and N_R single bond
///   radii", so UFF treats both ends of the amide bond as resonant; typing the carbon `C_2` would
///   give the amide C–N a radius the paper did not use. An N-acyl *aromatic* nitrogen (an
///   N-acetylimidazole) does not make its carbon `C_R` — see [`crate::energy`] for why that bond
///   is not treated as an amide.
/// - **Nitrogen.** Aromatic → `N_R`. Triple → `N_1`. Double → `N_2`. Otherwise, **an amide
///   nitrogen** — one with only single bonds, next to a carbon that has a double bond to oxygen —
///   → `N_R`; else `N_3`. The amide case is the one rule here that bond orders alone cannot
///   give: the C–N bond of an amide is drawn single and is partly double by resonance, which is
///   why an amide is planar, and the paper's resonant type is how UFF represents that: p. 10026
///   gives the amide C–N bond order as 1.41, chosen to reproduce the 1.366 Å bond of
///   N-methylformamide from "the C_R and N_R single bond radii" — read against the scanned paper
///   on 2026-10-01. Only C=O counts as the carbonyl — a thioamide (C=S) is typed
///   `N_3`, a choice the paper does not dictate.
/// - **Oxygen.** Aromatic → `O_R` (a furan's). Double → `O_2`. **Triple → `O_1`** (carbon
///   monoxide) — a choice: the specification this was written to gives no rule for `O_1`, and a
///   triple-bonded oxygen is what the type's 180° angle describes. Otherwise `O_3` — **unless it
///   is resonant**, below.
/// - **Sulfur.** Aromatic → `S_R` (a thiophene's). Double → `S_2`. Otherwise `S_3+2` — **unless
///   it is resonant**, below.
/// - **A resonant oxygen or sulfur: a divalent O or S with only single bonds, bonded to at least
///   one sp² or resonant atom, is `O_R` or `S_R`** — an ether, ester, enol or thioether oxygen or
///   sulfur whose lone pair is conjugated with a π system. "Sp² or resonant" is the neighbour's
///   type before this rule ([`UffType::hybridisation`] `Sp2`: `C_2`, `C_R`, `N_2`, `N_R`, `O_2`,
///   `S_2`, and the aromatic `O_R` and `S_R`), so the rule does not propagate along a chain.
///   **The paper's, for an oxygen on a C=C and on an aromatic ring**: `O_R`'s radius was "fit to
///   … methyl vinyl ether … C–O single … 1.428 Å" (p. 10025) and its angle to methyl vinyl
///   ether's 118.3° (p. 10028), which only enter if that ether oxygen is `O_R`; and anisole and
///   thioanisole are the paper's "tests of eq 17" (p. 10029), the sp²–sp² torsion, which only
///   applies if their ring–O and ring–S bonds join two sp² atoms. Measured before this rule was
///   adopted, typed so they relax to the paper's numbers: anisole's barrier 3.628 kcal/mol against
///   3.6, thioanisole's 1.666 against 1.7, methyl vinyl ether's C–O–C 118.05° against 118.3°;
///   typed `O_3`/`S_3+2`, 19.9, 14.5 and 107.8°. The paper also types a salicylidene ligand's
///   ring oxygen "resonating" in a nickel complex (p. 10034). **Choices, the paper having no
///   example**, each the same rule applied without exception:
///   - **an O bearing H** — a phenol's, a carboxylic acid's OH — is resonant too: its lone pair is
///     conjugated whatever its other neighbour is, and the paper's resonant phenolate oxygen is
///     the nearest it comes;
///   - **an ester's O** between the C=O and an alkyl carbon is resonant (methyl formate's
///     C–O–C relaxes to 115.6° against the paper's 113.6° so, and to 107.7° as `O_3`);
///   - **an O between two sp² atoms** — a diaryl ether, an anhydride — is resonant;
///   - **an O on an sp² nitrogen** — an oxime's — is resonant: the rule names the neighbour's
///     hybridisation, not its element;
///   - **on a ring, on a C=C or on a C=O alike**: the rule does not distinguish them;
///   - **an O bonded only to sp atoms and sp³ atoms stays `O_3`** — on a nitrile carbon, say:
///     an sp centre has no torsion, and the paper names sp², not sp;
///   - **sulfur follows oxygen** — a thioester's or thioanisole's S is `S_R`; the paper's `S_R`
///     radius is thiophene's, so a thioether's resonant S is the rule's, not the paper's.
///
///   The bond across an `O_R`–C or `S_R`–C bond keeps the dictionary's order — 1 for a single bond
///   — and only a bond the dictionary flags aromatic is 1.5, as before. A three-coordinate
///   oxygen (an oxonium) is never resonant; it is the one way left to reach the group-6 sp³–sp²
///   torsion rule (see [`crate::angular`]). Nitrogen is unchanged.
///   **Known gap:** a sulfoxide or sulfone sulfur — common in drugs — is hypervalent and UFF's
///   own types for it are `S_3+4` and `S_3+6`, which are not in this table. It is typed `S_2`
///   here because it has a double bond, and that is wrong for it. [`ForceField::new`] refuses
///   one — any sulfur whose bond orders sum past two — with
///   [`Unsupported::HypervalentSulfur`], rather than compute an energy from the wrong radius.
///
/// [`ForceField::new`]: crate::energy::ForceField::new
/// [`Unsupported::HypervalentSulfur`]: crate::energy::Unsupported::HypervalentSulfur
/// - **Phosphorus.** A double bond to oxygen (a phosphate, a phosphine oxide) → `P_3+5`;
///   otherwise `P_3+3`. The suffix is the oxidation state, and P=O is what puts phosphorus at +5
///   in every common drug motif.
pub fn assign(component: &Component) -> Vec<UffType> {
    let first = by_own_bonds(component);
    (0..first.len())
        .map(|i| resonant_heteroatom(component, &first, i).unwrap_or(first[i]))
        .collect()
}

/// Every atom's type from its own element and bonds — every rule of [`assign`] but the resonant
/// heteroatom's, which reads its neighbours' types from this.
fn by_own_bonds(component: &Component) -> Vec<UffType> {
    let atoms = component.atoms();
    (0..atoms.len())
        .map(|i| {
            let b = bonding(component, i);
            match atoms[i].element {
                Element::H => UffType::H,
                Element::F => UffType::F,
                Element::Cl => UffType::Cl,
                Element::Br => UffType::Br,
                Element::I => UffType::I,
                Element::C => {
                    if b.aromatic {
                        UffType::CR
                    } else if b.triples > 0 || b.doubles >= 2 {
                        UffType::C1
                    } else if is_amide_carbon(component, i) {
                        UffType::CR
                    } else if b.doubles > 0 {
                        UffType::C2
                    } else {
                        UffType::C3
                    }
                }
                Element::N => {
                    if b.aromatic {
                        UffType::NR
                    } else if b.triples > 0 {
                        UffType::N1
                    } else if b.doubles > 0 {
                        UffType::N2
                    } else if is_amide_nitrogen(component, i) {
                        UffType::NR
                    } else {
                        UffType::N3
                    }
                }
                Element::O => {
                    if b.aromatic {
                        UffType::OR
                    } else if b.triples > 0 {
                        UffType::O1
                    } else if b.doubles > 0 {
                        UffType::O2
                    } else {
                        UffType::O3
                    }
                }
                Element::S => {
                    if b.aromatic {
                        UffType::SR
                    } else if b.doubles > 0 {
                        UffType::S2
                    } else {
                        UffType::S3Divalent
                    }
                }
                Element::P => {
                    if b.double_to_oxygen {
                        UffType::P3Pentavalent
                    } else {
                        UffType::P3Trivalent
                    }
                }
            }
        })
        .collect()
}

/// `O_R` or `S_R` for a divalent oxygen or sulfur with only single bonds next to an sp² or
/// resonant atom — `None` for every other atom. The neighbours' types are `first`, the types
/// before this rule, so the rule does not feed on itself: in a peroxy ester `C(=O)–O–O–C` the
/// oxygen on the carbonyl is resonant and the one beyond it is not. See [`assign`] for why, and
/// for which of its branches are the paper's and which are choices.
fn resonant_heteroatom(component: &Component, first: &[UffType], i: usize) -> Option<UffType> {
    let resonant = match first[i] {
        UffType::O3 => UffType::OR,
        UffType::S3Divalent => UffType::SR,
        _ => return None,
    };
    let neighbours: Vec<usize> = component.neighbours(i).map(|(j, _)| j).collect();
    if neighbours.len() != 2 {
        return None;
    }
    let conjugated = neighbours
        .iter()
        .any(|&j| first[j].hybridisation() == Some(Hybridisation::Sp2));
    if !conjugated {
        return None;
    }
    Some(resonant)
}

/// A carbon with a C=O and a single, non-aromatic bond to a non-aromatic nitrogen whose bonds are
/// all single — the nitrogen [`assign`] types `N_R` by its amide rule.
fn is_amide_carbon(component: &Component, c: usize) -> bool {
    let atoms = component.atoms();
    let carbonyl = component
        .neighbours(c)
        .any(|(o, bond)| atoms[o].element == Element::O && bond.order == BondOrder::Double);
    carbonyl
        && component.neighbours(c).any(|(n, bond)| {
            if atoms[n].element != Element::N || bond.order != BondOrder::Single || bond.aromatic {
                return false;
            }
            let nb = bonding(component, n);
            !nb.aromatic && nb.doubles == 0 && nb.triples == 0
        })
}

/// A nitrogen with only single bonds, next to a carbon that carries a C=O.
fn is_amide_nitrogen(component: &Component, i: usize) -> bool {
    let atoms = component.atoms();
    component.neighbours(i).any(|(c, _)| {
        atoms[c].element == Element::C
            && component
                .neighbours(c)
                .any(|(o, bond)| atoms[o].element == Element::O && bond.order == BondOrder::Double)
    })
}
