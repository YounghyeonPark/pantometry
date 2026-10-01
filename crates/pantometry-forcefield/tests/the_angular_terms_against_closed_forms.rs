//! **Angle bend, torsion and inversion, each against a number that does not come from this
//! crate's code.**
//!
//! The angle forms' minimum and curvature at θ₀, which eq 12 fixes; eq 13 by hand for methane and
//! for an unsymmetric angle; ethane's nine torsions summing to Table III's barrier; eq 17 by hand
//! for ethylene and benzene; the collinear switch against its documented formula; the inversion
//! zero when planar, by hand when not, and phosphorus's 22 kcal/mol at planarity. Hand values were
//! evaluated separately from the crate (the expressions are written beside them) and typed here
//! as literals.

use pantometry_forcefield::angular::{
    bend_form, pyramid_inversion_angle, switch, torsion_parameters, Bend, BendForm, TorsionCase,
    TORSION_SWITCH_ANGLE,
};
use pantometry_forcefield::energy::natural_length;
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{Component, ForceField, UffType};
use std::f64::consts::PI;

const ANGSTROM: f64 = 1e-10;
const EPS: f64 = f64::EPSILON;
const DEG: f64 = PI / 180.0;
const AIN: &str = include_str!("../components/AIN.cif");

/// A one-component CCD file: atoms as `(name, element, aromatic)` at made-up places (only the
/// topology matters; every test passes its own positions), bonds as `(first, second, order,
/// aromatic)`.
fn entry(atoms: &[(&str, &str, bool)], bonds: &[(&str, &str, &str, bool)]) -> Component {
    let yn = |b: bool| if b { "Y" } else { "N" };
    let mut s = String::from(
        "data_TST\n_chem_comp.id TST\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\n",
    );
    for (k, (name, element, aromatic)) in atoms.iter().enumerate() {
        s += &format!("TST {name} {element} 0 {} {k} 0 0\n", yn(*aromatic));
    }
    s += "loop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n\
          _chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\n";
    for (a, b, order, aromatic) in bonds {
        s += &format!("{a} {b} {order} {}\n", yn(*aromatic));
    }
    Component::from_ccd(&s).expect("the test entry parses")
}

fn force_field(c: &Component) -> ForceField {
    ForceField::new(c, &uff::assign(c)).expect("supported")
}

fn metres(p: [f64; 3]) -> [f64; 3] {
    [p[0] * ANGSTROM, p[1] * ANGSTROM, p[2] * ANGSTROM]
}

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

// ---------------------------------------------------------------- angle bend

/// **Which form each centre takes**, by its θ₀ — the decision `angular` documents. The ones that
/// would differ under the other reading, "by hybridisation": `N_2`, `O_R` and `S_R` are sp² types
/// whose θ₀ is not 120°, and p. 10028 says `N_2`'s and `O_R`'s were fitted, so they take the
/// general form.
#[test]
fn each_centre_takes_the_form_its_natural_angle_names() {
    for t in UffType::ALL {
        let form = bend_form(t);
        let expected = match t {
            UffType::C1
            | UffType::N1
            | UffType::O1
            | UffType::H
            | UffType::F
            | UffType::Cl
            | UffType::Br
            | UffType::I => "linear",
            UffType::CR | UffType::C2 | UffType::NR | UffType::O2 | UffType::S2 => "trigonal",
            _ => "general",
        };
        let got = match form {
            BendForm::Linear => "linear",
            BendForm::Trigonal => "trigonal",
            BendForm::General { .. } => "general",
        };
        assert_eq!(got, expected, "{t}");
    }
}

/// The largest `|E⁽ᵏ⁾(θ)|` of a form, k = 3 or 4, as a multiple of K: for the general form
/// `|c₁| + 2ᵏ |c₂|`; for the linear `1`; for the trigonal `3ᵏ / 9`.
fn derivative_bound(form: BendForm, k: i32) -> f64 {
    match form {
        BendForm::General { c1, c2, .. } => c1.abs() + 2f64.powi(k) * c2.abs(),
        BendForm::Linear => 1.0,
        BendForm::Trigonal => 3f64.powi(k) / 9.0,
    }
}

/// The sum of the magnitudes of a form's parts, as a multiple of K — what its rounding scales
/// with.
fn magnitude(form: BendForm) -> f64 {
    match form {
        BendForm::General { c0, c1, c2 } => c0.abs() + c1.abs() + 2.0 * c2.abs(),
        BendForm::Linear => 2.0,
        BendForm::Trigonal => 6.0 / 9.0,
    }
}

/// **Every form is zero, flat and of curvature K at θ₀, for every centre type.** Eq 12's
/// coefficients make the general form satisfy all three — `C₀ + C₁ c + C₂ (2c² − 1) = 0`,
/// `C₁ s + 4 C₂ s c = 0`, and `−C₁ c − 4 C₂ (2c² − 1) = 4 C₂ s² = K/K` — and the special forms do
/// by their own algebra (`d²/dθ² [1 + cos θ] = 1` at π, `d²/dθ² [(1 − cos 3θ)/9] = 1` at 2π/3).
/// Checked numerically, not assumed:
///
/// - `E(θ₀)`: a cancellation among parts of total size `M K` (`magnitude`), each good to a few ε;
///   `8 ε M K`.
/// - The slope, `(E(θ₀ + h) − E(θ₀ − h)) / 2h`: truncation `h²/6 |E‴|` plus rounding
///   `8 ε M K / h`, h = 1e-4.
/// - The curvature, `(E(θ₀ + h) − 2 E(θ₀) + E(θ₀ − h)) / h²` against K: truncation
///   `h²/12 |E⁗|` plus rounding `4 · 8 ε M K / h²`, h = 1e-3.
///
/// And it is a minimum: 0.1 rad either side is higher.
#[test]
fn every_form_is_zero_flat_and_of_curvature_k_at_its_natural_angle() {
    for t in UffType::ALL {
        let r = natural_length(t, UffType::H, 1.0);
        let b = Bend::new([0, 1, 2], [UffType::H, t, UffType::H], r, r);
        let (k, t0) = (b.force_constant, b.natural_angle);
        let m = magnitude(b.form) * k;
        let e0 = b.energy(t0);
        assert!(e0.abs() <= 8.0 * EPS * m, "{t}: E(θ0) = {e0:e}");

        let h = 1e-4;
        let slope = (b.energy(t0 + h) - b.energy(t0 - h)) / (2.0 * h);
        let tol = h * h / 6.0 * derivative_bound(b.form, 3) * k + 8.0 * EPS * m / h;
        assert!(slope.abs() <= tol, "{t}: slope {slope:e} > {tol:e}");

        let h = 1e-3;
        let curvature = (b.energy(t0 + h) - 2.0 * e0 + b.energy(t0 - h)) / (h * h);
        let tol = h * h / 12.0 * derivative_bound(b.form, 4) * k + 32.0 * EPS * m / (h * h);
        assert!(
            (curvature - k).abs() <= tol,
            "{t}: curvature {curvature:e} vs K {k:e}, tol {tol:e}"
        );

        assert!(b.energy(t0 - 0.1) > e0, "{t}");
        if t0 < PI {
            assert!(b.energy(t0 + 0.1) > e0, "{t}");
        }
    }
}

/// **The linear form has its minimum at 180°, as Table I's θ₀ says it must — eq 10 as printed
/// does not.** With n = 1, eq 10 reads `K (1 − cos θ)`: zero at 0° and `2K` at 180°. The crate's
/// `K (1 + cos θ)` is exactly zero at π (`cos π` is −1 in floating point) and exactly `2K` at 0.
/// And the general form has the maximum at 180° the paper asks of it for water (p. 10027):
/// `E″(π) = K (C₁ − 4 C₂) = −4 K C₂ (1 + cos θ₀) < 0`.
#[test]
fn a_linear_centre_is_straight_and_a_bent_one_is_not() {
    let r = natural_length(UffType::C1, UffType::H, 1.0);
    let b = Bend::new([0, 1, 2], [UffType::H, UffType::C1, UffType::H], r, r);
    assert_eq!(b.form, BendForm::Linear);
    assert_eq!(b.energy(PI), 0.0);
    assert_eq!(b.energy(0.0), 2.0 * b.force_constant);
    assert!(b.energy(PI - 0.1) > 0.0);

    let r = natural_length(UffType::O3, UffType::H, 1.0);
    let water = Bend::new([0, 1, 2], [UffType::H, UffType::O3, UffType::H], r, r);
    let BendForm::General { c1, c2, .. } = water.form else {
        panic!("O_3 is general");
    };
    assert!(c1 - 4.0 * c2 < 0.0, "a maximum at 180°");
    assert!(water.energy(PI) > water.energy(PI - 0.1));
}

/// **Eq 13 by hand**, worked outside the crate from Table I, the χ of the bond step, and β =
/// 664.12 / (r_IJ r_JK):
///
/// - H–`C_3`–H with `r_CH` = 1.109400794877744 Å (eq 2, the bond step's value) and θ₀ = 109.47°:
///   `r_IK² = 2 r² (1 − cos θ₀)`, `K = 664.12 Z_H² / r_IK⁵ [3 r² (1 − cos² θ₀) − r_IK² cos θ₀]` =
///   75.4987656584479 kcal mol⁻¹ rad⁻².
/// - H–`C_3`–Cl, `r_CCl` = 1.7779854853733603 Å: 102.20457149190283. The Z* are the two *ends'*;
///   using the centre's for one of them (`Z_H Z_C`) would give 83.226, far outside.
///
/// And methane through [`ForceField`] has six bends, each with the first K. 1e-13 relative is a
/// few dozen roundings.
#[test]
fn the_angle_force_constant_by_hand() {
    let close = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-13;
    let k = |b: &Bend| b.force_constant / KCAL_PER_MOL;
    let ch = natural_length(UffType::C3, UffType::H, 1.0);
    let ccl = natural_length(UffType::C3, UffType::Cl, 1.0);
    let hch = Bend::new([0, 1, 2], [UffType::H, UffType::C3, UffType::H], ch, ch);
    assert!(close(k(&hch), 75.498_765_658_447_9), "{}", k(&hch));
    let hccl = Bend::new([0, 1, 2], [UffType::H, UffType::C3, UffType::Cl], ch, ccl);
    assert!(close(k(&hccl), 102.204_571_491_902_83), "{}", k(&hccl));
    let clch = Bend::new([0, 1, 2], [UffType::Cl, UffType::C3, UffType::H], ccl, ch);
    assert!(close(k(&clch), k(&hccl)), "symmetric in the two ends");

    let methane = entry(
        &[
            ("C", "C", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
            ("H4", "H", false),
        ],
        &[
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
            ("C", "H3", "SING", false),
            ("C", "H4", "SING", false),
        ],
    );
    let ff = force_field(&methane);
    assert_eq!(ff.bends().len(), 6);
    for b in ff.bends() {
        assert_eq!(b.atoms[1], 0, "every angle is centred on the carbon");
        assert!(close(k(b), 75.498_765_658_447_9));
    }
    assert!(ff.torsions().is_empty() && ff.inversions().is_empty());
}

// ---------------------------------------------------------------- torsion

/// **Every row of the torsion table, by hand.** V in kcal/mol, worked outside the crate:
/// sp³–sp³ `C_3` `√(2.119²)` = 2.119; `O_3`–`O_3` group 6, `√(2 · 2)` = 2; `O_3`–`S_3+2`,
/// `√(2 · 6.8)` = 3.687817782917155; eq 17 with U = 2 for the second row and 1.25 for the third:
/// `C_2`=`C_2` at order 2, `10 (1 + 4.18 ln 2)` = 38.97355214740571; `C_R`–`C_R` at 1.5,
/// 26.94844151892127; the amide `C_R`–`N_R` at 1.41, 24.362049643505213; `O_3`–`C_R`, group 6
/// sp³ on sp², order 1, 10; `S_3+2`–`C_2`, `5 √(1.25 · 2)` = 7.905694150420948; propene 2;
/// otherwise sp³–sp² 1. And none about a bond to an sp atom or to a hydrogen.
#[test]
fn the_torsion_rules_by_hand() {
    let c = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-14;
    let rows = [
        (
            UffType::C3,
            UffType::C3,
            1.0,
            false,
            false,
            TorsionCase::Sp3Sp3,
            2.119,
            3,
            180.0,
        ),
        (
            UffType::N3,
            UffType::C3,
            1.0,
            false,
            false,
            TorsionCase::Sp3Sp3,
            (0.45f64 * 2.119).sqrt(),
            3,
            180.0,
        ),
        (
            UffType::O3,
            UffType::O3,
            1.0,
            false,
            false,
            TorsionCase::Group6Sp3Sp3,
            2.0,
            2,
            90.0,
        ),
        (
            UffType::O3,
            UffType::S3Divalent,
            1.0,
            false,
            false,
            TorsionCase::Group6Sp3Sp3,
            3.687_817_782_917_155,
            2,
            90.0,
        ),
        (
            UffType::C2,
            UffType::C2,
            2.0,
            true,
            true,
            TorsionCase::Sp2Sp2,
            38.973_552_147_405_71,
            2,
            180.0,
        ),
        (
            UffType::CR,
            UffType::CR,
            1.5,
            true,
            true,
            TorsionCase::Sp2Sp2,
            26.948_441_518_921_27,
            2,
            180.0,
        ),
        (
            UffType::CR,
            UffType::NR,
            1.41,
            true,
            true,
            TorsionCase::Sp2Sp2,
            24.362_049_643_505_213,
            2,
            180.0,
        ),
        (
            UffType::O3,
            UffType::CR,
            1.0,
            false,
            true,
            TorsionCase::Group6Sp3Sp2,
            10.0,
            2,
            90.0,
        ),
        (
            UffType::C2,
            UffType::S3Divalent,
            1.0,
            true,
            false,
            TorsionCase::Group6Sp3Sp2,
            7.905_694_150_420_948,
            2,
            90.0,
        ),
        (
            UffType::C3,
            UffType::C2,
            1.0,
            false,
            true,
            TorsionCase::Propene,
            2.0,
            3,
            180.0,
        ),
        (
            UffType::C2,
            UffType::C3,
            1.0,
            true,
            false,
            TorsionCase::Propene,
            2.0,
            3,
            180.0,
        ),
        (
            UffType::C3,
            UffType::C2,
            1.0,
            false,
            false,
            TorsionCase::Sp2Sp3,
            1.0,
            6,
            0.0,
        ),
        // The group-6 row is tried before the propene row.
        (
            UffType::O3,
            UffType::C2,
            1.0,
            false,
            true,
            TorsionCase::Group6Sp3Sp2,
            10.0,
            2,
            90.0,
        ),
        // The group-6 row is for an sp2 atom of *another* column (p. 10028); on a group-6 sp2
        // atom the sp3 oxygen or sulfur falls through to the propene or general row.
        (
            UffType::O3,
            UffType::S2,
            1.0,
            false,
            true,
            TorsionCase::Propene,
            2.0,
            3,
            180.0,
        ),
        (
            UffType::S3Divalent,
            UffType::OR,
            1.0,
            false,
            false,
            TorsionCase::Sp2Sp3,
            1.0,
            6,
            0.0,
        ),
    ];
    for (j, k, order, js, ks, case, v, n, phi0) in rows {
        let p = torsion_parameters(j, k, order, js, ks).unwrap_or_else(|| panic!("{j}-{k}"));
        assert_eq!(p.case, case, "{j}-{k}");
        assert!(c(kcal(p.barrier), v), "{j}-{k}: {} vs {v}", kcal(p.barrier));
        assert_eq!(p.periodicity, n, "{j}-{k}");
        assert_eq!(p.equilibrium, phi0 * DEG, "{j}-{k}");
    }
    for (j, k) in [
        (UffType::C1, UffType::C3),
        (UffType::C3, UffType::N1),
        (UffType::O1, UffType::C2),
        (UffType::H, UffType::C3),
        (UffType::C3, UffType::Cl),
    ] {
        assert!(
            torsion_parameters(j, k, 1.0, true, true).is_none(),
            "{j}-{k}"
        );
    }
}

/// **The propene test is read off the molecule: propene's methyl–CH bond is the exception, and
/// methyl isocyanate's methyl–N bond is the general sp²–sp³ case**, because the isocyanate
/// nitrogen (`N_2`)'s other neighbour is the sp carbon of N=C=O, not an sp² atom. Propene's bond
/// is shared six ways (three hydrogens times the CH's hydrogen and CH₂ carbon), the isocyanate's
/// three; and nothing turns about the N=C or C=O bonds of the sp carbon.
#[test]
fn the_propene_exception_reads_the_neighbours() {
    let propene = entry(
        &[
            ("C1", "C", false),
            ("C2", "C", false),
            ("C3", "C", false),
            ("H11", "H", false),
            ("H12", "H", false),
            ("H13", "H", false),
            ("H2", "H", false),
            ("H31", "H", false),
            ("H32", "H", false),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C2", "C3", "DOUB", false),
            ("C1", "H11", "SING", false),
            ("C1", "H12", "SING", false),
            ("C1", "H13", "SING", false),
            ("C2", "H2", "SING", false),
            ("C3", "H31", "SING", false),
            ("C3", "H32", "SING", false),
        ],
    );
    let ff = force_field(&propene);
    let about = |ff: &ForceField, j: usize, k: usize| -> Vec<_> {
        ff.torsions()
            .iter()
            .filter(|t| [t.atoms[1], t.atoms[2]] == [j, k])
            .map(|t| (t.parameters.case, t.torsions_about_bond))
            .collect()
    };
    assert_eq!(about(&ff, 0, 1), vec![(TorsionCase::Propene, 6); 6]);

    let isocyanate = entry(
        &[
            ("C1", "C", false),
            ("N", "N", false),
            ("C2", "C", false),
            ("O", "O", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
        ],
        &[
            ("C1", "N", "SING", false),
            ("N", "C2", "DOUB", false),
            ("C2", "O", "DOUB", false),
            ("C1", "H1", "SING", false),
            ("C1", "H2", "SING", false),
            ("C1", "H3", "SING", false),
        ],
    );
    assert_eq!(
        uff::assign(&isocyanate)[..4],
        [UffType::C3, UffType::N2, UffType::C1, UffType::O2]
    );
    let ff = force_field(&isocyanate);
    assert_eq!(about(&ff, 0, 1), vec![(TorsionCase::Sp2Sp3, 3); 3]);
    assert_eq!(ff.torsions().len(), 3, "none about the sp carbon's bonds");
}

/// **In a three-membered ring, I = L is not a torsion.** About each ring bond of cyclopropane,
/// I can be the third ring carbon or either of J's hydrogens and L the third ring carbon or
/// either of K's: nine pairs, one of which puts the third carbon at both ends. So eight torsions
/// per ring bond, twenty-four in all, each sharing its bond's barrier eight ways.
#[test]
fn a_three_ring_has_eight_torsions_per_bond() {
    let mut atoms = vec![("C1", "C", false), ("C2", "C", false), ("C3", "C", false)];
    let hs = ["H11", "H12", "H21", "H22", "H31", "H32"];
    atoms.extend(hs.iter().map(|h| (*h, "H", false)));
    let mut bonds = vec![
        ("C1", "C2", "SING", false),
        ("C2", "C3", "SING", false),
        ("C3", "C1", "SING", false),
    ];
    for (k, h) in hs.iter().enumerate() {
        bonds.push((["C1", "C2", "C3"][k / 2], *h, "SING", false));
    }
    let ff = force_field(&entry(&atoms, &bonds));
    assert_eq!(ff.torsions().len(), 24);
    for t in ff.torsions() {
        assert_eq!(t.torsions_about_bond, 8);
        assert_ne!(t.atoms[0], t.atoms[3]);
    }
}

/// Ethane's topology: C1, C2, then H11, H12, H13 on C1 and H21, H22, H23 on C2.
fn ethane() -> Component {
    entry(
        &[
            ("C1", "C", false),
            ("C2", "C", false),
            ("H11", "H", false),
            ("H12", "H", false),
            ("H13", "H", false),
            ("H21", "H", false),
            ("H22", "H", false),
            ("H23", "H", false),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C1", "H11", "SING", false),
            ("C1", "H12", "SING", false),
            ("C1", "H13", "SING", false),
            ("C2", "H21", "SING", false),
            ("C2", "H22", "SING", false),
            ("C2", "H23", "SING", false),
        ],
    )
}

/// Ethane at its natural lengths and angles — C–C 1.514 Å along z, C–H `r_CH`, every H–C–C at
/// `C_3`'s θ₀ — with the second methyl turned by `psi` degrees from eclipsed.
fn ethane_at(psi: f64) -> Vec<[f64; 3]> {
    let cc = natural_length(UffType::C3, UffType::C3, 1.0) / ANGSTROM;
    let ch = natural_length(UffType::C3, UffType::H, 1.0) / ANGSTROM;
    let beta = UffType::C3.parameters().angle;
    let (sb, cb) = beta.sin_cos();
    let mut at = vec![metres([0.0, 0.0, 0.0]), metres([0.0, 0.0, cc])];
    for a in [0.0, 120.0, 240.0] {
        let (s, c) = (a * DEG).sin_cos();
        at.push(metres([ch * sb * c, ch * sb * s, ch * cb]));
    }
    for a in [0.0, 120.0, 240.0] {
        let (s, c) = ((a + psi) * DEG).sin_cos();
        at.push(metres([ch * sb * c, ch * sb * s, cc - ch * cb]));
    }
    at
}

/// **Ethane's nine H–C–C–H torsions share one barrier, and sum to exactly Table III's 2.119
/// kcal/mol eclipsed and to zero staggered.** Each is `½ (V/9) [1 − cos 540° cos 3φ]`, and with
/// the second methyl turned by ψ every one of the nine has `cos 3φ = cos 3ψ`, so the sum is
/// `½ V (1 + cos 3ψ)` — V at ψ = 0, 0 at 60°, V/2 at 30°, and 1.785 kcal/mol at 17°.
///
/// Tolerance: each term's `cos φ` comes from two cross products of bond vectors that meet at
/// 109.47° (sin 0.943), good to about `20 ε / 0.943`; `|dE/d cos φ| ≤ ½ (V/9) n²` with n = 3; nine
/// terms give `V · 9/2 · 21.2 ε` ≈ `95 ε V`, and the positions, built from sines and cosines of ψ,
/// are good to a few ε more. `200 ε V`.
#[test]
fn ethanes_nine_torsions_sum_to_the_tables_barrier() {
    let ff = force_field(&ethane());
    let v = 2.119;
    assert_eq!(ff.torsions().len(), 9);
    for t in ff.torsions() {
        assert_eq!(t.torsions_about_bond, 9);
        assert_eq!(t.parameters.case, TorsionCase::Sp3Sp3);
        assert!((kcal(t.parameters.barrier) / v - 1.0).abs() < 1e-15);
    }
    let tol = 200.0 * EPS * v;
    for psi in [0.0, 60.0, 30.0, 17.0, 180.0, 300.0] {
        let e = kcal(ff.energy(&ethane_at(psi)).torsion);
        let expected = 0.5 * v * (1.0 + (3.0 * psi * DEG).cos());
        assert!(
            (e - expected).abs() <= tol,
            "ψ = {psi}: {e} vs {expected}, tol {tol:e}"
        );
    }
}

/// **Rigid rotation of ethane's methyl: the torsion term's barrier is 2.119 kcal/mol, and the
/// total is reported — not compared with the paper's 2.90.** Table II's "calculated" barriers are
/// relaxed: the saddle found by hill climbing, every other coordinate minimised (p. 10031,
/// III.A). A rigid rotation holds the bonds and angles at their natural values, so only the
/// torsion and the six 1-4 H···H van der Waals pairs change. The relaxed comparison needs the
/// minimiser and is step 1d's.
///
/// Asserted: the torsion difference is V (the tolerance of the test above, twice); the bond,
/// angle and inversion terms do not change beyond rounding (a rigid rotation moves no bond length
/// or angle; a few hundred ε of their own size); and the total difference is the torsion plus the
/// van der Waals difference, to the summation's rounding.
#[test]
fn ethanes_rigid_barrier() {
    let ff = force_field(&ethane());
    let (ecl, stag) = (ff.energy(&ethane_at(0.0)), ff.energy(&ethane_at(60.0)));
    let d = |f: fn(&pantometry_forcefield::Energy) -> f64| kcal(f(&ecl) - f(&stag));
    let torsion = d(|e| e.torsion);
    let vdw = d(|e| e.van_der_waals);
    let total = d(|e| e.total);
    let fixed = d(|e| e.bond + e.angle + e.inversion);
    println!(
        "ethane, rigid rotation at natural lengths and angles: eclipsed − staggered = \
         torsion {torsion:.6} + van der Waals {vdw:.6} + bond/angle/inversion {fixed:.1e} = \
         total {total:.6} kcal/mol (Table II: experiment 2.93, UFF relaxed 2.90)"
    );
    assert!((torsion - 2.119).abs() <= 400.0 * EPS * 2.119, "{torsion}");
    let size = kcal(stag.bond.abs() + stag.angle.abs()) + 1.0;
    assert!(fixed.abs() <= 400.0 * EPS * size, "{fixed:e}");
    let s = kcal(ecl.total.abs() + stag.total.abs()) + 10.0;
    assert!((total - torsion - vdw - fixed).abs() <= 100.0 * EPS * s);
}

/// **Ethylene's four torsions sum to eq 17's V at a 90° twist, and are exactly zero planar.**
/// `5 √(2 · 2) (1 + 4.18 ln 2)` = 38.97355214740571 kcal/mol, each torsion carrying a quarter as
/// `½ (V/4) (1 − cos 2φ)`; twisted by τ every φ is τ or τ + 180°, so the sum is `½ V (1 − cos 2τ)`.
/// Planar, with every z exactly 0, both cross products lie exactly along z and `cos φ` is exactly
/// ±1, so the sum is exactly 0. Twisted, the tolerance is the ethane test's with n = 2 and four
/// terms: `100 ε V`.
#[test]
fn ethylenes_barrier_is_eq_17() {
    let c = entry(
        &[
            ("C1", "C", false),
            ("C2", "C", false),
            ("H11", "H", false),
            ("H12", "H", false),
            ("H21", "H", false),
            ("H22", "H", false),
        ],
        &[
            ("C1", "C2", "DOUB", false),
            ("C1", "H11", "SING", false),
            ("C1", "H12", "SING", false),
            ("C2", "H21", "SING", false),
            ("C2", "H22", "SING", false),
        ],
    );
    let ff = force_field(&c);
    let v = 38.973_552_147_405_71;
    assert_eq!(ff.torsions().len(), 4);
    let cc = natural_length(UffType::C2, UffType::C2, 2.0) / ANGSTROM;
    let ch = natural_length(UffType::C2, UffType::H, 1.0) / ANGSTROM;
    let at = |tau: f64| {
        let (s, co) = (tau * DEG).sin_cos();
        let (h, w) = (ch * 0.5, ch * 3f64.sqrt() / 2.0);
        vec![
            metres([0.0, 0.0, 0.0]),
            metres([cc, 0.0, 0.0]),
            metres([-h, w, 0.0]),
            metres([-h, -w, 0.0]),
            metres([cc + h, w * co, w * s]),
            metres([cc + h, -w * co, -w * s]),
        ]
    };
    assert_eq!(ff.energy(&at(0.0)).torsion, 0.0);
    for tau in [90.0, 45.0, 30.0] {
        let e = kcal(ff.energy(&at(tau)).torsion);
        let expected = 0.5 * v * (1.0 - (2.0 * tau * DEG).cos());
        assert!(
            (e - expected).abs() <= 100.0 * EPS * v,
            "τ = {tau}: {e} vs {expected}"
        );
    }
}

/// **Benzene's ring torsions carry eq 17 at order 1.5, and the flat ring has exactly no torsion or
/// inversion energy.** `10 (1 + 4.18 ln 1.5)` = 26.94844151892127 kcal/mol, shared four ways per
/// ring bond. The ring and its hydrogens at z = 0 make every `cos φ` exactly ±1 and every
/// inversion's out-of-plane sine exactly 0.
#[test]
fn benzenes_ring_torsion_is_eq_17_at_order_one_and_a_half() {
    let names: Vec<String> = (1..=6).map(|k| format!("C{k}")).collect();
    let hs: Vec<String> = (1..=6).map(|k| format!("H{k}")).collect();
    let mut atoms: Vec<(&str, &str, bool)> =
        names.iter().map(|n| (n.as_str(), "C", true)).collect();
    atoms.extend(hs.iter().map(|n| (n.as_str(), "H", false)));
    let mut bonds = Vec::new();
    for k in 0..6 {
        let order = if k % 2 == 0 { "DOUB" } else { "SING" };
        bonds.push((names[k].as_str(), names[(k + 1) % 6].as_str(), order, true));
        bonds.push((names[k].as_str(), hs[k].as_str(), "SING", false));
    }
    let c = entry(&atoms, &bonds);
    let ff = force_field(&c);
    let ring = ff
        .torsions()
        .iter()
        .filter(|t| t.atoms[1] < 6 && t.atoms[2] < 6)
        .count();
    assert_eq!(ring, 24, "four per ring bond");
    for t in ff.torsions() {
        assert_eq!(t.parameters.case, TorsionCase::Sp2Sp2);
        assert_eq!(t.torsions_about_bond, 4);
        assert!((kcal(t.parameters.barrier) / 26.948_441_518_921_27 - 1.0).abs() < 1e-14);
    }
    assert_eq!(ff.inversions().len(), 18, "three per ring carbon");
    let at: Vec<[f64; 3]> = (0..12)
        .map(|i| {
            let r = if i < 6 { 1.39 } else { 2.47 };
            let a = (60 * (i % 6)) as f64 * DEG;
            metres([r * a.cos(), r * a.sin(), 0.0])
        })
        .collect();
    let e = ff.energy(&at);
    assert_eq!((e.torsion, e.inversion), (0.0, 0.0));
}

/// **The torsion is switched off continuously as a central angle straightens** — the convention
/// p. 10029 states, done without a jump. A four-carbon chain C1–C2–C3–C4, one torsion about
/// C2–C3 with V = 2.119 kcal/mol and n = 3, at dihedral 40° so that unswitched it is
/// `½ V (1 + cos 120°)` = V/4, with the angle C1–C2–C3 at θ:
///
/// - up to 170° it is V/4, unswitched;
/// - past it, V/4 times `s = t² (3 − 2t)`, `t = (1 − |cos θ|) / (1 − cos 10°)`, typed here from
///   the documentation and not taken from the crate;
/// - straight, exactly 0;
/// - across 170° ± 1e-7 rad it moves by no more than `s`'s own change there, `3 Δt²` with
///   `Δt = 1e-7 sin 10° / (1 − cos 10°)` — about 4e-12 of V/4, so the switch adds no step;
/// - at 0.001 rad short of straight it is below `3 t²` of V/4, about 3e-9 of it.
///
/// Rounding throughout: `cos φ` from cross products of bonds meeting at up to 179.9°, good to
/// `20 ε / sin θ`, times `|dE/d cos φ| ≤ ½ V n²`.
#[test]
fn the_torsion_fades_continuously_as_an_angle_straightens() {
    let c = entry(
        &[
            ("C1", "C", false),
            ("C2", "C", false),
            ("C3", "C", false),
            ("C4", "C", false),
        ],
        &[
            ("C1", "C2", "SING", false),
            ("C2", "C3", "SING", false),
            ("C3", "C4", "SING", false),
        ],
    );
    let ff = force_field(&c);
    assert_eq!(ff.torsions().len(), 1);
    let t = ff.torsions()[0];
    let v = kcal(t.parameters.barrier);
    let r = 1.514;
    let phi = 40.0 * DEG;
    let at = |c1: [f64; 3]| {
        let (s, co) = (110.0 * DEG).sin_cos();
        vec![
            metres(c1),
            metres([0.0, 0.0, 0.0]),
            metres([r, 0.0, 0.0]),
            metres([r - r * co, r * s * phi.cos(), r * s * phi.sin()]),
        ]
    };
    let bent = |theta: f64| at([r * theta.cos(), r * theta.sin(), 0.0]);
    let e = |p: &[[f64; 3]]| kcal(t.energy_at(p));
    let quarter = v / 4.0;
    let on = (10.0 * DEG).cos();
    let s_hand = |theta: f64| {
        let tt = (1.0 - theta.cos().abs()) / (1.0 - on);
        if tt >= 1.0 {
            1.0
        } else {
            tt * tt * (3.0 - 2.0 * tt)
        }
    };
    let round = |theta: f64| 20.0 * EPS / theta.sin().abs() * 0.5 * v * 9.0 + 4.0 * EPS * v;
    for deg in [100.0, 150.0, 169.0, 169.999] {
        let th = deg * DEG;
        assert!((e(&bent(th)) - quarter).abs() <= round(th), "{deg}°");
    }
    assert_eq!(TORSION_SWITCH_ANGLE, 170.0 * DEG, "the documented edge");
    // The same chain read backwards puts the bent angle at K, J–K–L: the other half of the
    // switch. The term is symmetric under reversal, so the expectations are the same.
    let reversed = |p: Vec<[f64; 3]>| p.into_iter().rev().collect::<Vec<_>>();
    for deg in [170.5, 172.0, 175.0, 178.0, 179.0, 179.9] {
        let th = deg * DEG;
        let expected = quarter * s_hand(th);
        for (half, got) in [("I-J-K", e(&bent(th))), ("J-K-L", e(&reversed(bent(th))))] {
            assert!(got > 0.0 && got < quarter, "{half} {deg}°: {got}");
            assert!(
                (got - expected).abs() <= round(th),
                "{half} {deg}°: {got} vs {expected}"
            );
        }
    }
    assert_eq!(e(&at([-r, 0.0, 0.0])), 0.0, "I-J-K straight");
    assert_eq!(e(&reversed(at([-r, 0.0, 0.0]))), 0.0, "J-K-L straight");

    let delta = 1e-7;
    let dt = delta * (10.0 * DEG).sin() / (1.0 - on);
    let th = 170.0 * DEG;
    let step = (e(&bent(th - delta)) - e(&bent(th + delta))).abs();
    assert!(
        step <= quarter * 3.0 * dt * dt + 2.0 * round(th),
        "a step of {step:e} at 170°"
    );

    let th = PI - 1e-3;
    let tt = (1.0 - th.cos().abs()) / (1.0 - on);
    let near = e(&bent(th));
    assert!(
        near >= 0.0 && near <= quarter * 3.0 * tt * tt + round(th),
        "{near:e} at 0.001 rad from straight"
    );
}

// ---------------------------------------------------------------- inversion

/// Three bonds of lengths `r` from the origin at azimuths 0°, 120°, 240°, each raised `alpha`
/// out of the xy plane: a three-fold pyramid, flat at `alpha` = 0.
fn pyramid(alpha: f64, r: [f64; 3]) -> Vec<[f64; 3]> {
    let mut at = vec![[0.0; 3]];
    let (sa, ca) = alpha.sin_cos();
    for (k, a) in [0.0, 120.0, 240.0].iter().enumerate() {
        let (s, c) = (a * DEG).sin_cos();
        at.push(metres([r[k] * ca * c, r[k] * ca * s, r[k] * sa]));
    }
    at
}

/// **Formaldehyde's carbon has K = 50 kcal/mol and methanimine's K = 6, both exactly zero planar,
/// and both by hand when pyramidal** (p. 10029). At z = 0 the plane normal is exactly along z and
/// each bond exactly in the plane, so each sin ω is exactly 0. Raised 10° out of the plane, each
/// bond makes ω = 29.425400140682804° with the plane of the other two (`asin` of the triple
/// product over the cross product's length, computed outside the crate), and the three terms
/// together are `K (1 − cos ω)`: 6.450195000624487 and 0.7740234000749384 kcal/mol. The bond
/// lengths are unequal on purpose: ω does not depend on them.
///
/// Tolerance: sin ω from a cross product of bonds 117° apart (sin 0.89), good to `20 ε / 0.89`;
/// `|dE/d sin ω| = K tan ω`, so `K · 0.56 · 22.5 ε` ≈ `13 ε K`, about `100 ε` of the energy; the
/// positions from sines and cosines add a few ε. `500 ε` relative.
#[test]
fn a_trigonal_carbon_resists_pyramidalisation_by_hand() {
    let formaldehyde = entry(
        &[
            ("C", "C", false),
            ("O", "O", false),
            ("H1", "H", false),
            ("H2", "H", false),
        ],
        &[
            ("C", "O", "DOUB", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
        ],
    );
    let methanimine = entry(
        &[
            ("C", "C", false),
            ("N", "N", false),
            ("H1", "H", false),
            ("H2", "H", false),
        ],
        &[
            ("C", "N", "DOUB", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
        ],
    );
    for (c, k, by_hand) in [
        (formaldehyde, 50.0, 6.450_195_000_624_487),
        (methanimine, 6.0, 0.774_023_400_074_938_4),
    ] {
        let ff = force_field(&c);
        assert_eq!(ff.inversions().len(), 3);
        for v in ff.inversions() {
            assert_eq!(v.atoms[0], 0, "centred on the carbon");
            assert!((kcal(v.force_constant) / k - 1.0).abs() < 1e-15);
            assert_eq!(v.coefficients, [1.0, -1.0, 0.0]);
        }
        let r = [1.21, 1.10, 1.09];
        assert_eq!(ff.energy(&pyramid(0.0, r)).inversion, 0.0);
        let e = kcal(ff.energy(&pyramid(10.0 * DEG, r)).inversion);
        assert!(
            (e / by_hand - 1.0).abs() <= 500.0 * EPS,
            "K {k}: {e} vs {by_hand}"
        );
    }
}

/// **The three inversions are taken about three different axes.** In a three-fold pyramid all
/// three ω are equal, so a term counted twice about one bond and never about another would pass
/// the test above. Here formaldehyde's bonds are at azimuth 0°, 115° and 245° and elevation 25°,
/// −5° and 0° (lengths 1.21, 1.10, 1.09 Å), which gives three different ω — 19.0677°, 15.7618°,
/// 16.0375°, outside the crate — and `50/3 Σ (1 − cos ωᵢ)` = 2.1897720417560977 kcal/mol. Counting
/// the O axis twice in place of either H's would give 2.478 or 2.456; the tolerance is the one
/// above.
#[test]
fn the_three_inversions_use_three_axes() {
    let c = entry(
        &[
            ("C", "C", false),
            ("O", "O", false),
            ("H1", "H", false),
            ("H2", "H", false),
        ],
        &[
            ("C", "O", "DOUB", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
        ],
    );
    let bond = |r: f64, az: f64, el: f64| {
        let (a, e) = (az * DEG, el * DEG);
        metres([r * e.cos() * a.cos(), r * e.cos() * a.sin(), r * e.sin()])
    };
    let at = vec![
        [0.0; 3],
        bond(1.21, 0.0, 25.0),
        bond(1.10, 115.0, -5.0),
        bond(1.09, 245.0, 0.0),
    ];
    let e = kcal(force_field(&c).energy(&at).inversion);
    let by_hand = 2.189_772_041_756_097_7;
    assert!((e / by_hand - 1.0).abs() <= 500.0 * EPS, "{e} vs {by_hand}");
}

/// **Phosphine's inversion: zero at PH₃'s own angle and 22 kcal/mol flat** — the paper's two
/// conditions (p. 10029), reached through the crate's ω₀. ω₀ for a pyramid with bond angle 93.8°
/// is 84.43386488851856° (`asin(√(1 − 3c² + 2c³) / sin θ)`, computed outside the crate). A pyramid
/// whose bonds meet at θ has elevation `α` with `sin² α = (cos θ + ½) / (3/2)`. At that geometry
/// the three terms are a cancellation among parts of size `K (|C₀| + |C₁| + |C₂|)`, good to
/// `30 ε` of it; flat, they are `K (C₀ + C₁ + C₂)` = 22 to the same.
#[test]
fn a_phosphine_inverts_over_twenty_two_kcal_per_mol() {
    let w0 = pyramid_inversion_angle(93.8 * DEG);
    assert!(
        (w0 / DEG - 84.433_864_888_518_56).abs() < 1e-12,
        "{}",
        w0 / DEG
    );
    let c = entry(
        &[
            ("P", "P", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
        ],
        &[
            ("P", "H1", "SING", false),
            ("P", "H2", "SING", false),
            ("P", "H3", "SING", false),
        ],
    );
    assert_eq!(uff::assign(&c)[0], UffType::P3Trivalent);
    let ff = force_field(&c);
    assert_eq!(ff.inversions().len(), 3);
    let v = ff.inversions()[0];
    let size = kcal(v.force_constant) * v.coefficients.iter().map(|x| x.abs()).sum::<f64>();
    let alpha = ((((93.8 * DEG).cos() + 0.5) / 1.5).sqrt()).asin();
    let r = [1.42, 1.42, 1.42];
    let rest = kcal(ff.energy(&pyramid(alpha, r)).inversion);
    assert!(rest.abs() <= 30.0 * EPS * size, "{rest:e}");
    let flat = kcal(ff.energy(&pyramid(0.0, r)).inversion);
    assert!((flat - 22.0).abs() <= 30.0 * EPS * size, "{flat}");
}

/// **No nitrogen has an inversion term** — "E_barrier = 0 for nitrogen" (p. 10029) — neither
/// ammonia's `N_3` nor an amide's `N_R`; formamide's carbon, `C_R` bonded to `O_2`, has K = 50.
/// And aspirin has eight trigonal carbons, the two carbonyl ones at 50 and the six ring ones at 6.
#[test]
fn nitrogen_has_no_inversion_and_aspirin_has_eight_centres() {
    let ammonia = entry(
        &[
            ("N", "N", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
        ],
        &[
            ("N", "H1", "SING", false),
            ("N", "H2", "SING", false),
            ("N", "H3", "SING", false),
        ],
    );
    assert!(force_field(&ammonia).inversions().is_empty());
    let formamide = entry(
        &[
            ("C", "C", false),
            ("O", "O", false),
            ("N", "N", false),
            ("H", "H", false),
            ("HN1", "H", false),
            ("HN2", "H", false),
        ],
        &[
            ("C", "O", "DOUB", false),
            ("C", "N", "SING", false),
            ("C", "H", "SING", false),
            ("N", "HN1", "SING", false),
            ("N", "HN2", "SING", false),
        ],
    );
    assert_eq!(uff::assign(&formamide)[2], UffType::NR);
    let ff = force_field(&formamide);
    assert_eq!(ff.inversions().len(), 3);
    assert!(ff.inversions().iter().all(|v| v.atoms[0] == 0));
    assert!(ff
        .inversions()
        .iter()
        .all(|v| (kcal(v.force_constant) - 50.0).abs() < 1e-12));

    let c = Component::from_ccd(AIN).expect("AIN parses");
    let ff = force_field(&c);
    assert_eq!(ff.inversions().len(), 24);
    let at50: Vec<&str> = ff
        .inversions()
        .iter()
        .filter(|v| (kcal(v.force_constant) - 50.0).abs() < 1e-12)
        .map(|v| c.atoms()[v.atoms[0]].name.as_str())
        .collect();
    assert_eq!(at50, ["C7", "C7", "C7", "C8", "C8", "C8"]);
}

/// **The trigonal form away from its minimum: `E(180°) = E(60°) = 2K/9` and `E(0) = 0`** —
/// `K/9 (1 − cos 3θ)` at 3θ = 540°, 180° and 0. The tests at θ₀ pin only the neighbourhood of
/// 120°, where eq 11 with θ₀ = 120° (C₀ = ½, C₁ = ⅔, C₂ = ⅓) agrees to second order; it gives
/// K/6 at 180° and 3K/2 at 0. At π and 0 the cosine is exactly −1 and 1 in floating point; at
/// π/3 it is good to an ulp. `8 ε K`.
#[test]
fn the_trigonal_form_away_from_its_minimum() {
    let r = natural_length(UffType::CR, UffType::CR, 1.5);
    let b = Bend::new([0, 1, 2], [UffType::CR; 3], r, r);
    assert_eq!(b.form, BendForm::Trigonal);
    let k = b.force_constant;
    let tol = 8.0 * EPS * k;
    for theta in [PI, PI / 3.0] {
        let e = b.energy(theta);
        assert!((e - 2.0 * k / 9.0).abs() <= tol, "{theta}: {e:e}");
    }
    assert!(b.energy(0.0).abs() <= tol);
}

/// **The bend's two lengths are the bonds' natural lengths at their own orders**, read off the
/// force field rather than built by hand. Eq 2 and eq 13 worked outside the crate:
///
/// - aspirin's ring angles, `C_R`–`C_R`–`C_R` at order 1.5: r = 1.379256405400789 Å each and
///   K = 222.59501669477606 kcal mol⁻¹ rad⁻² (188.44 if the bonds were taken single);
/// - formaldehyde's H–C=O, `C_2` with r_CH = 1.0844161478393872 Å and r_C=O = 1.2194547241544453
///   Å at order 2: K = 170.395684267617 (143.69 with a single C–O).
///
/// 1e-13 relative, as in the other eq 13 test.
#[test]
fn the_bend_reads_its_bonds_orders() {
    let close = |a: f64, b: f64| (a / b - 1.0).abs() < 1e-13;
    let c = Component::from_ccd(AIN).expect("AIN parses");
    let types = uff::assign(&c);
    let ff = force_field(&c);
    let ring: Vec<&Bend> = ff
        .bends()
        .iter()
        .filter(|b| b.atoms.iter().all(|&a| types[a] == UffType::CR))
        .collect();
    assert_eq!(ring.len(), 6, "one per ring carbon");
    for b in ring {
        assert!(close(kcal(b.force_constant), 222.595_016_694_776_06));
    }

    let formaldehyde = entry(
        &[
            ("C", "C", false),
            ("O", "O", false),
            ("H1", "H", false),
            ("H2", "H", false),
        ],
        &[
            ("C", "O", "DOUB", false),
            ("C", "H1", "SING", false),
            ("C", "H2", "SING", false),
        ],
    );
    let ff = force_field(&formaldehyde);
    let hco: Vec<&Bend> = ff.bends().iter().filter(|b| b.atoms.contains(&1)).collect();
    assert_eq!(hco.len(), 2);
    for b in hco {
        assert!(close(kcal(b.force_constant), 170.395_684_267_617), "{b:?}");
    }
}

/// **The n = 6 row evaluated: methyl isocyanate's three H–C–N=C torsions sum to
/// `½ V (1 − cos 6ψ)`**, V = 1 kcal/mol, with the methyl turned by ψ from eclipsing the N=C bond:
/// each is `½ (V/3) [1 − cos 0° cos 6φ]` and every φ is ψ plus a multiple of 120°, so every
/// `cos 6φ` is `cos 6ψ`. 0 at ψ = 0, 0.25 kcal/mol at 10°, 0.9330127 at 25°. Tolerance as for
/// ethane with n = 6: three terms, `|dE/d cos φ| ≤ ½ (V/3) 36`, `cos φ` good to `20 ε / 0.91`,
/// so about `400 ε V`, and the positions' sines and cosines a few ε more: `1000 ε V`.
#[test]
fn methyl_isocyanates_six_fold_torsion() {
    let c = entry(
        &[
            ("C1", "C", false),
            ("N", "N", false),
            ("C2", "C", false),
            ("O", "O", false),
            ("H1", "H", false),
            ("H2", "H", false),
            ("H3", "H", false),
        ],
        &[
            ("C1", "N", "SING", false),
            ("N", "C2", "DOUB", false),
            ("C2", "O", "DOUB", false),
            ("C1", "H1", "SING", false),
            ("C1", "H2", "SING", false),
            ("C1", "H3", "SING", false),
        ],
    );
    let ff = force_field(&c);
    assert!(ff
        .torsions()
        .iter()
        .all(|t| t.parameters.case == TorsionCase::Sp2Sp3));
    let (sb, cb) = (109.47 * DEG).sin_cos();
    let at = |psi: f64| {
        let mut p = vec![
            metres([-1.47, 0.0, 0.0]),
            metres([0.0, 0.0, 0.0]),
            metres([0.55, 1.05, 0.0]),
            metres([1.093, 2.087, 0.0]),
        ];
        for a in [0.0, 120.0, 240.0] {
            let (s, co) = ((psi + a) * DEG).sin_cos();
            p.push(metres([-1.47 + 1.09 * cb, 1.09 * sb * co, 1.09 * sb * s]));
        }
        p
    };
    for psi in [0.0, 10.0, 25.0] {
        let e = kcal(ff.energy(&at(psi)).torsion);
        let expected = 0.5 * (1.0 - (6.0 * psi * DEG).cos());
        assert!(
            (e - expected).abs() <= 1000.0 * EPS,
            "ψ = {psi}: {e} vs {expected}"
        );
    }
}

/// **The switch's slope against a central difference of the switch, on both ends** — near
/// `x = +1` (a straight angle) and near `x = −1` (an angle folded to 0°), where the slope has the
/// opposite sign. `s` is `t² (3 − 2t)` in `t = (1 − |x|) / (1 − cos 10°)`, so `|s‴| ≤ 12 / (1 −
/// cos 10°)³` ≈ 3.4e6. The difference at h = 1e-6, divided by the step actually taken (`x ± h`
/// rounds by ε|x|, which at a slope of 90 is 1e-7 of error at h = 1e-7), is good to
/// `h²/6 |s‴|` for truncation plus `2 · 8 ε / h` for two evaluations of `s ≤ 1` at a few roundings
/// each — about 6e-7, against slopes of order 50 to 90.
#[test]
fn the_switch_slope_on_both_ends() {
    let h = 1e-6;
    let s3 = 12.0 / (1.0 - (10.0 * DEG).cos()).powi(3);
    let tol = h * h / 6.0 * s3 + 16.0 * EPS / h;
    for x in [0.99, 0.995, -0.99, -0.995] {
        let (_, ds) = switch(x);
        let (up, down) = (x + h, x - h);
        let numeric = (switch(up).0 - switch(down).0) / (up - down);
        assert!(
            (ds - numeric).abs() <= tol,
            "x = {x}: {ds} vs {numeric}, tol {tol:e}"
        );
        assert!(
            ds.abs() > 10.0 && ds.signum() == -x.signum(),
            "x = {x}: {ds}"
        );
    }
}

/// **The planar reading of the group-6 rule moves its minimum and nothing else**: φ₀ = 180°
/// instead of 90°, with the same n = 2 and eq 17's V; and on every other row the two readings
/// agree. It exists to measure the open question (see `the_papers_barriers`), so what it changes
/// is pinned here.
#[test]
fn the_planar_reading_moves_only_the_minimum() {
    use pantometry_forcefield::angular::{torsion_parameters_with, Group6OnSp2};
    for (j, k, order, js, ks) in [
        (UffType::O3, UffType::CR, 1.0, false, true),
        (UffType::C2, UffType::S3Divalent, 1.0, true, false),
        (UffType::C3, UffType::C3, 1.0, false, false),
        (UffType::CR, UffType::CR, 1.5, true, true),
        (UffType::C3, UffType::C2, 1.0, false, true),
    ] {
        let printed = torsion_parameters_with(j, k, order, js, ks, Group6OnSp2::AsPrinted)
            .expect("a torsion");
        let planar =
            torsion_parameters_with(j, k, order, js, ks, Group6OnSp2::Planar).expect("a torsion");
        assert_eq!(printed.case, planar.case);
        assert_eq!(printed.barrier, planar.barrier);
        assert_eq!(printed.periodicity, planar.periodicity);
        if printed.case == TorsionCase::Group6Sp3Sp2 {
            assert_eq!(printed.equilibrium, 90.0 * DEG);
            assert_eq!(planar.equilibrium, 180.0 * DEG);
        } else {
            assert_eq!(printed.equilibrium, planar.equilibrium);
        }
    }
}
