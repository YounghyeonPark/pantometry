//! **Table I of Rappé et al. (1992), and what geometry and the force field's own construction say
//! about it.**
//!
//! Two kinds of check. A spot-check of rows against the paper, typed here a second time and
//! separately from `uff.rs`, so a slip in either copy shows up as a disagreement — which catches a
//! typing error and not a misreading of the scan, since both copies came from one transcription.
//! And properties the paper's numbers must have
//! whatever their digits: the natural angles of `C_3`, `C_R` and `C_1` are the tetrahedral,
//! trigonal and linear angles; UFF's non-bonded parameters belong to the element, not the
//! hybridisation; and the SI conversion inverts.

use pantometry_forcefield::uff::{AVOGADRO, ELEMENTARY_CHARGE, KCAL_PER_MOL};
use pantometry_forcefield::{Element, UffType};

/// Rows of Table I as printed: r1 Å, θ0 °, x1 Å, D1 kcal/mol, ζ, Z1 e.
const PAPER: [(UffType, [f64; 6]); 6] = [
    (UffType::H, [0.354, 180.0, 2.886, 0.044, 12.0, 0.712]),
    (UffType::C3, [0.757, 109.47, 3.851, 0.105, 12.73, 1.912]),
    (UffType::NR, [0.699, 120.0, 3.660, 0.069, 13.407, 2.544]),
    (UffType::O3, [0.658, 104.51, 3.500, 0.060, 14.085, 2.300]),
    (
        UffType::P3Pentavalent,
        [1.056, 109.47, 4.147, 0.305, 13.072, 2.863],
    ),
    (UffType::I, [1.382, 180.0, 4.50, 0.339, 15.0, 2.65]),
];

#[test]
fn rows_match_the_paper_as_printed() {
    for (t, row) in PAPER {
        let p = t.table_i();
        assert_eq!([p.r1, p.theta0, p.x1, p.d1, p.zeta, p.z1], row, "{t}");
    }
}

/// **SI is the paper's row times exact factors**: 1 Å = 1e-10 m, 1 kcal/mol = 4184 J / N_A, and
/// 1 e = 1.602 176 634e-19 C. Relative 1e-15 is a few ulps, which is what two multiplications
/// cost.
#[test]
fn si_is_the_printed_row_converted_exactly() {
    assert!((KCAL_PER_MOL * AVOGADRO / 4184.0 - 1.0).abs() < 1e-15);
    let close = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-15;
    for t in UffType::ALL {
        let (row, si) = (t.table_i(), t.parameters());
        assert!(close(si.bond_radius, row.r1 * 1e-10), "{t} r1");
        assert!(close(si.angle, row.theta0.to_radians()), "{t} theta0");
        assert!(close(si.vdw_distance, row.x1 * 1e-10), "{t} x1");
        assert!(
            close(si.vdw_energy, row.d1 * 4184.0 / 6.022_140_76e23),
            "{t} D1"
        );
        assert_eq!(si.vdw_scale, row.zeta, "{t} zeta");
        assert!(
            close(si.effective_charge, row.z1 * ELEMENTARY_CHARGE),
            "{t} Z1"
        );
    }
    // One absolute anchor: carbon's well is 0.105 kcal/mol = 7.295e-22 J per molecule.
    assert!((UffType::C3.parameters().vdw_energy - 7.2951e-22).abs() < 1e-26);
}

/// **The natural angles are geometry.** sp³ is the tetrahedral angle `acos(−1/3)` = 109.4712°,
/// printed to two decimals, so it agrees to the printing's half-unit 0.005°. sp² is 120° and sp
/// is 180°, exactly. `P_3+5` is tetrahedral too: a phosphate is.
#[test]
fn the_natural_angles_are_the_angles_of_the_hybridisation() {
    let tetrahedral = (-1.0f64 / 3.0).acos().to_degrees();
    for t in [UffType::C3, UffType::P3Pentavalent] {
        assert!((t.table_i().theta0 - tetrahedral).abs() < 0.005, "{t}");
    }
    for t in [
        UffType::CR,
        UffType::C2,
        UffType::NR,
        UffType::O2,
        UffType::S2,
    ] {
        assert_eq!(t.table_i().theta0, 120.0, "{t}");
    }
    for t in [UffType::C1, UffType::N1, UffType::O1] {
        assert_eq!(t.table_i().theta0, 180.0, "{t}");
    }
    // Water's and ammonia's measured angles are what `O_3` and `N_3` take, which is the paper's
    // stated choice: 104.51 deg and 106.7 deg.
    assert_eq!(UffType::O3.table_i().theta0, 104.51);
    assert_eq!(UffType::N3.table_i().theta0, 106.7);
}

/// **Non-bonded parameters belong to the element.** UFF's van der Waals distance, well depth and
/// shape, and its effective charge, are per element: every hybridisation of carbon has the same
/// four. A transcription slip in one row breaks this.
#[test]
fn every_hybridisation_of_an_element_shares_its_nonbonded_parameters() {
    for element in Element::ALL {
        let rows: Vec<_> = UffType::ALL
            .into_iter()
            .filter(|t| t.element() == element)
            .map(|t| {
                let r = t.table_i();
                [r.x1, r.d1, r.zeta, r.z1]
            })
            .collect();
        assert!(!rows.is_empty(), "{element} has no type");
        assert!(rows.iter().all(|r| *r == rows[0]), "{element}: {rows:?}");
    }
}

/// **Twenty-two types, each with its own label, and the torsion barriers the paper gives.**
#[test]
fn the_types_and_the_torsion_barriers() {
    let mut labels: Vec<&str> = UffType::ALL.iter().map(|t| t.label()).collect();
    labels.sort_unstable();
    labels.dedup();
    assert_eq!(labels.len(), 22);

    let kcal = |t: UffType| t.sp3_torsion_barrier().map(|v| v / KCAL_PER_MOL);
    let near = |a: Option<f64>, b: f64| a.is_some_and(|a| (a - b).abs() < 1e-12);
    assert!(near(kcal(UffType::C3), 2.119));
    assert!(near(kcal(UffType::N3), 0.450));
    assert!(near(kcal(UffType::O3), 0.018));
    assert!(near(kcal(UffType::P3Trivalent), 2.400));
    assert!(near(kcal(UffType::P3Pentavalent), 2.400));
    assert!(near(kcal(UffType::S3Divalent), 0.484));
    assert_eq!(UffType::CR.sp3_torsion_barrier(), None);

    let u = |t: UffType| t.sp2_torsion_barrier().map(|v| v / KCAL_PER_MOL);
    assert!(near(u(UffType::C2), 2.0));
    assert!(near(u(UffType::S2), 1.25));
    assert!(near(u(UffType::Br), 0.7));
    assert!(near(u(UffType::I), 0.2));
    assert_eq!(UffType::H.sp2_torsion_barrier(), None);
}
