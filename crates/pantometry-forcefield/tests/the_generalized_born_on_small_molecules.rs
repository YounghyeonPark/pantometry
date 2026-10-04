//! **What OBC II with QEq charges says about real molecules' hydration, beside experiment.**
//!
//! Each molecule is relaxed in vacuum by UFF from the dictionary's ideal coordinates (zero charges,
//! the crate's default), given QEq charges at that minimum, and then scored without moving: the
//! GB electrostatic solvation energy of eq 2, the nonpolar `0.005 kcal mol⁻¹ Å⁻² × SASA` of p. 392
//! (Bondi radii + 1.4 Å probe, 1600 points), and their sum — against the experimental hydration
//! free energy in **FreeSolv** (D. L. Mobley and J. P. Guthrie, *J. Comput.-Aided Mol. Des.* **28**,
//! 711 (2014); database v0.52, `database.txt` as fetched from
//! `github.com/MobleyLab/FreeSolv` at commit `6c7d19b4`, 2026-10-04; each value's own experimental
//! reference is beside it below).
//!
//! **Reported, not asserted against experiment.** Nothing in this model was fitted to these
//! molecules: QEq's charges are not GB's (OBC was fitted with AMBER's), and the radii were fitted to
//! Poisson–Boltzmann energies of proteins. No agreement is earned in advance, so none is asserted.
//! What is asserted is what holds by construction: every number finite, the polar part negative,
//! and the nonpolar part positive.

use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry_forcefield::solvation::{self, GeneralizedBorn, PROBE_RADIUS};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::{qeq, Component, Element, ForceField, Status};

/// `(file, name, FreeSolv id, experiment kcal/mol, its uncertainty, its reference)`.
const MOLECULES: [(&str, &str, &str, f64, f64, &str); 18] = [
    (
        "MOH",
        "methanol",
        "mobley_1636752",
        -5.10,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "NME",
        "methylamine",
        "mobley_6714389",
        -4.55,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "DMN",
        "dimethylamine",
        "mobley_5692472",
        -4.29,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "KEN",
        "trimethylamine",
        "mobley_9209581",
        -3.20,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "MEE",
        "methanethiol",
        "mobley_525934",
        -1.20,
        0.20,
        "10.1039/P29900000291",
    ),
    (
        "2F2",
        "dimethyl ether",
        "mobley_7015518",
        -1.91,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "2ME",
        "methyl ethyl ether",
        "mobley_6911232",
        -2.10,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "ACE",
        "acetaldehyde",
        "mobley_1967551",
        -3.50,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "ACN",
        "acetone",
        "mobley_3867265",
        -3.80,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "CCN",
        "acetonitrile",
        "mobley_7532833",
        -3.88,
        0.60,
        "10.1021/jp0264477",
    ),
    (
        "ACM",
        "acetamide",
        "mobley_8048190",
        -9.71,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "DMF",
        "N,N-dimethylformamide",
        "mobley_8011706",
        -7.81,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "TME",
        "propane",
        "mobley_2068538",
        2.00,
        0.20,
        "10.1039/P29900000291",
    ),
    (
        "61G",
        "isoprene",
        "mobley_8765203",
        0.68,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "PYJ",
        "ethylbenzene",
        "mobley_8127829",
        -0.79,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "A1JFW",
        "anisole",
        "mobley_7295828",
        -2.45,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "16R",
        "thioanisole",
        "mobley_1189457",
        -2.73,
        0.60,
        "10.1021/ct050097l",
    ),
    (
        "AIN",
        "aspirin",
        "mobley_2913224",
        -9.94,
        0.18,
        "10.1007/s10822-010-9350-8",
    ),
];

fn text(file: &str) -> &'static str {
    match file {
        "MOH" => include_str!("../components/MOH.cif"),
        "NME" => include_str!("../components/NME.cif"),
        "DMN" => include_str!("../components/DMN.cif"),
        "KEN" => include_str!("../components/KEN.cif"),
        "MEE" => include_str!("../components/MEE.cif"),
        "2F2" => include_str!("../components/2F2.cif"),
        "2ME" => include_str!("../components/2ME.cif"),
        "ACE" => include_str!("../components/ACE.cif"),
        "ACN" => include_str!("../components/ACN.cif"),
        "CCN" => include_str!("../components/CCN.cif"),
        "ACM" => include_str!("../components/ACM.cif"),
        "DMF" => include_str!("../components/DMF.cif"),
        "TME" => include_str!("../components/TME.cif"),
        "61G" => include_str!("../components/61G.cif"),
        "PYJ" => include_str!("../components/PYJ.cif"),
        "A1JFW" => include_str!("../components/A1JFW.cif"),
        "16R" => include_str!("../components/16R.cif"),
        "AIN" => include_str!("../components/AIN.cif"),
        other => panic!("no component {other}"),
    }
}

/// Points per atom for the surface: from the two-sphere measurement in
/// `the_generalized_born_against_closed_forms.rs`, within 0.3% of the exact area of eight
/// overlapping pairs there: 0.005 kcal/mol on these molecules' nonpolar terms.
const POINTS: usize = 1600;

struct Row {
    polar: f64,
    nonpolar: f64,
    area: f64,
}

fn score(file: &str) -> Row {
    let c = Component::from_ccd(text(file)).expect("parses");
    let ff = ForceField::new(&c, &uff::assign(&c)).expect("supported");
    let mut at: Vec<[f64; 3]> = c.atoms().iter().map(|a| a.at).collect();
    let p = ff.minimise(&mut at, 50_000, 1e-4 * KCAL_PER_MOL_ANGSTROM);
    assert_eq!(p.status, Status::Converged, "{file} relaxes");
    let q = qeq::charges(&c, &at).expect("QEq").charges;
    let elements: Vec<Element> = c.atoms().iter().map(|a| a.element).collect();
    let gb = GeneralizedBorn::new(&elements);
    let polar = gb.energy(&q, &at) / KCAL_PER_MOL;
    let radii: Vec<f64> = elements
        .iter()
        .map(|&e| solvation::intrinsic_radius(e))
        .collect();
    let area: f64 = solvation::surface_area(&radii, &at, PROBE_RADIUS, POINTS)
        .iter()
        .sum();
    let nonpolar = solvation::nonpolar_energy(area) / KCAL_PER_MOL;
    // Through the force field too, the way a caller would: the same number.
    let wet = ff
        .with_qeq_charges(&at)
        .expect("QEq")
        .with_generalized_born();
    let through = wet.evaluate(&at).energy.solvation / KCAL_PER_MOL;
    assert!((through - polar).abs() <= 1e-12 * polar.abs());
    Row {
        polar,
        nonpolar,
        area: area / 1e-20,
    }
}

#[test]
fn hydration_free_energies_beside_freesolv() {
    println!(
        "{:22} {:>9} {:>8} {:>9} {:>9} {:>11} {:>8}",
        "molecule", "ΔG_GB", "SASA Å²", "ΔG_surf", "ΔG_solv", "experiment", "error"
    );
    let mut errors = Vec::new();
    let mut polar_errors = Vec::new();
    let mut pairs = Vec::new();
    for (file, name, id, experiment, uncertainty, reference) in MOLECULES {
        let r = score(file);
        assert!(r.polar.is_finite() && r.polar < 0.0, "{name}: {}", r.polar);
        assert!(r.nonpolar > 0.0);
        let total = r.polar + r.nonpolar;
        println!(
            "{name:22} {:9.3} {:8.1} {:9.3} {:9.3} {experiment:6.2} ±{uncertainty:.2} {:8.3}   {id} {reference}",
            r.polar,
            r.area,
            r.nonpolar,
            total,
            total - experiment
        );
        errors.push(total - experiment);
        polar_errors.push(r.polar - experiment);
        pairs.push((total, experiment));
    }
    let stats = |e: &[f64]| {
        let n = e.len() as f64;
        let mean = e.iter().sum::<f64>() / n;
        let rms = (e.iter().map(|x| x * x).sum::<f64>() / n).sqrt();
        (mean, rms)
    };
    let (mse, rmse) = stats(&errors);
    let (pmse, prmse) = stats(&polar_errors);
    let n = pairs.len() as f64;
    let (mx, my) = (
        pairs.iter().map(|p| p.0).sum::<f64>() / n,
        pairs.iter().map(|p| p.1).sum::<f64>() / n,
    );
    let cov: f64 = pairs.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum();
    let vx: f64 = pairs.iter().map(|p| (p.0 - mx).powi(2)).sum();
    let vy: f64 = pairs.iter().map(|p| (p.1 - my).powi(2)).sum();
    println!(
        "{} molecules: ΔG_solv − experiment mean {mse:.3}, RMS {rmse:.3} kcal/mol, Pearson r {:.3}; \
         ΔG_GB alone mean {pmse:.3}, RMS {prmse:.3}",
        pairs.len(),
        cov / (vx * vy).sqrt()
    );
}
