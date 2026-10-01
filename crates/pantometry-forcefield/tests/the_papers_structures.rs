//! **The paper's minimised structures, Figures 4–7 (pp. 10032–33), against this crate's, under
//! both signs of `r_EN`.**
//!
//! The figures give UFF's calculated bond lengths and angles for small organic molecules, beside
//! experiment in parentheses; read off the scanned pages, enlarged. Each molecule here is
//! minimised from the dictionary's ideal coordinates (or a hand-built start, for the ones the
//! dictionary lacks — `tests/hand_built/`, each file says so) to 1e-5 kcal mol⁻¹ Å⁻¹, which leaves
//! every length within about 1e-7 Å and every angle within 1e-5° of the minimum — far below the
//! paper's three decimals and one.
//!
//! **The question this answers** is the one `energy` left open: whether eq 2's `r_EN` is
//! subtracted (this crate, Open Babel, RDKit) or added (as printed). Its answer here is not close.
//! For every bond the sign changes the natural length of — every bond between different elements —
//! the minus is nearer the paper's own number than the plus: 24 of 24, and the plus is off by up
//! to 0.04 Å (C=O, where `r_EN` is largest). **And it is the sign, not only the size, that is
//! measured**: a comparison of minus with plus alone would also come out 24 of 24 if the paper had
//! used no correction at all, so the correction dropped entirely is the third column. The minus is
//! nearer than no correction on all 24, and no correction nearer than the plus on all 24: the
//! paper's numbers sit where the full correction subtracted puts them. All three counts are
//! asserted. Two rows are
//! left out of that count and say why: methyl vinyl ether and methyl formate typed as this crate
//! types them, whose ether oxygen the paper types `O_R` (see [`Row::compare_sign`]); with it typed
//! so they are in, and the minus wins those too. With them typed `O_3`, the plus is nearer on two
//! of their six bonds — what the 0.022 Å between the two oxygen radii does, not the sign.
//!
//! Agreement within half the last printed figure (0.0005 Å, 0.05°) is asserted for the
//! measures that have it — a list decided after measuring, with the tolerance fixed before — and
//! every measure is printed. Means of a carbon's C–H bonds are printed against the figure's
//! value for one of them and are never asserted.

mod common;

use common::{angle, distance, index, positions, relaxed};
use pantometry_forcefield::energy::{ElectronegativitySign, Variant};
use pantometry_forcefield::{uff, Component, Element, ForceField, UffType};

#[derive(Clone, Copy, Debug)]
enum Measure {
    Bond(&'static str, &'static str),
    Angle(&'static str, &'static str, &'static str),
    /// The mean of a carbon's C–H bonds.
    MeanCh(&'static str),
    /// The shortest of a carbon's C–H bonds, for a methyl the figure labels with two values.
    MinCh(&'static str),
    /// The longest of a carbon's C–H bonds, likewise.
    MaxCh(&'static str),
}

struct Row {
    molecule: &'static str,
    text: &'static str,
    retype: Option<(&'static str, UffType)>,
    /// Whether this row enters the comparison of the two signs. Two rows do not: methyl vinyl
    /// ether and methyl formate with their ether oxygen typed `O_3`, as this crate types it. The
    /// paper types such an oxygen `O_R` — its `O_R` angle was fitted to methyl vinyl ether's C–O–C
    /// (p. 10028), and typed so, that angle comes out 118.05° against the paper's 118.3°, typed
    /// `O_3` 107.8° — and `O_R`'s radius is 0.022 Å longer, so their C–O bonds measure the typing
    /// before the sign. The same molecules typed `O_R` are in the comparison instead.
    compare_sign: bool,
    /// The measure, the paper's value as printed, and whether this crate's value is asserted to
    /// agree with it to half the last figure.
    measures: Vec<(Measure, &'static str, bool)>,
}

fn rows() -> Vec<Row> {
    use Measure::*;
    let r = |molecule, text, measures| Row {
        molecule,
        text,
        retype: None,
        compare_sign: true,
        measures,
    };
    vec![
        r(
            "propane",
            include_str!("../components/TME.cif"),
            vec![
                (Bond("C1", "C2"), "1.526", true),
                (Angle("C1", "C2", "C3"), "111.3", true),
                (MeanCh("C1"), "1.111", false),
            ],
        ),
        r(
            "propene",
            include_str!("hand_built/ppe.cif"),
            vec![
                (Bond("C1", "C2"), "1.336", true),
                (Bond("C2", "C3"), "1.500", true),
                (Angle("C1", "C2", "C3"), "122.0", true),
                (MeanCh("C3"), "1.111", false),
            ],
        ),
        r(
            "butadiene",
            include_str!("hand_built/bde.cif"),
            vec![
                (Bond("C1", "C2"), "1.330", false),
                (Bond("C2", "C3"), "1.471", false),
                (Angle("C1", "C2", "C3"), "121.1", true),
            ],
        ),
        r(
            "propyne",
            include_str!("hand_built/pyn.cif"),
            vec![
                (Bond("C2", "C3"), "1.463", true),
                (Bond("C1", "C2"), "1.205", true),
                (Bond("C1", "H1"), "1.058", true),
                (MeanCh("C3"), "1.109", false),
            ],
        ),
        r(
            "dimethylamine",
            include_str!("../components/DMN.cif"),
            vec![
                (Bond("N1", "HN1"), "1.047", true),
                (Bond("N1", "C2"), "1.464", true),
                (Angle("C2", "N1", "C3"), "110.2", true),
                (MeanCh("C2"), "1.110", false),
            ],
        ),
        r(
            "trimethylamine",
            include_str!("../components/KEN.cif"),
            vec![
                (Bond("CA1", "NE1"), "1.471", true),
                (Angle("CA1", "NE1", "CB1"), "110.0", true),
                (MeanCh("CA1"), "1.112", false),
            ],
        ),
        r(
            "dimethyldiazene",
            include_str!("hand_built/dmd.cif"),
            vec![
                (Bond("N1", "N2"), "1.247", true),
                (Bond("N1", "C1"), "1.442", true),
                (Angle("C1", "N1", "N2"), "112.3", false),
                (MeanCh("C1"), "1.112", false),
            ],
        ),
        r(
            "acetonitrile",
            include_str!("../components/CCN.cif"),
            vec![
                (Bond("C1", "C2"), "1.463", true),
                (Bond("N", "C1"), "1.157", true),
                (MeanCh("C2"), "1.109", false),
            ],
        ),
        r(
            "acrylonitrile",
            include_str!("../components/6AC.cif"),
            vec![
                (Bond("C2", "C3"), "1.334", true),
                (Bond("C1", "C2"), "1.441", true),
                (Bond("N1", "C1"), "1.157", true),
                (Angle("C3", "C2", "C1"), "120.8", true),
            ],
        ),
        r(
            "dimethyl ether",
            include_str!("../components/2F2.cif"),
            vec![
                (Bond("C", "O"), "1.410", true),
                (Angle("C", "O", "C1"), "109.2", true),
                // Fig 6 labels this methyl's C–H 1.109 and 1.113, not one value.
                (MinCh("C"), "1.109", false),
                (MaxCh("C"), "1.113", false),
            ],
        ),
        r(
            "methyl ethyl ether",
            include_str!("../components/2ME.cif"),
            vec![
                (Bond("CA'", "CB'"), "1.521", true),
                (Bond("CB'", "OC'"), "1.411", true),
                (Bond("OC'", "CD'"), "1.410", true),
                (Angle("CB'", "OC'", "CD'"), "109.3", true),
            ],
        ),
        Row {
            compare_sign: false,
            ..r(
                "methyl vinyl ether, O typed O_3",
                include_str!("hand_built/mve.cif"),
                vec![
                    (Bond("O", "CM"), "1.428", false),
                    (Angle("C2", "O", "CM"), "118.3", false),
                    (Bond("C2", "O"), "1.413", false),
                    (Bond("C1", "C2"), "1.343", false),
                ],
            )
        },
        Row {
            retype: Some(("O", UffType::OR)),
            ..r(
                "methyl vinyl ether, O typed O_R",
                include_str!("hand_built/mve.cif"),
                vec![
                    (Bond("O", "CM"), "1.428", true),
                    (Angle("C2", "O", "CM"), "118.3", false),
                    (Bond("C2", "O"), "1.413", false),
                    (Bond("C1", "C2"), "1.343", true),
                ],
            )
        },
        r(
            "acetaldehyde",
            include_str!("../components/ACE.cif"),
            vec![
                (Bond("C", "O"), "1.221", false),
                (Bond("C", "CH3"), "1.494", false),
                (Angle("CH3", "C", "O"), "119.8", false),
                (Bond("C", "H"), "1.085", true),
            ],
        ),
        r(
            "acetone",
            include_str!("../components/ACN.cif"),
            vec![
                (Bond("C", "O"), "1.222", true),
                (Bond("C", "C1"), "1.498", true),
                (Angle("C1", "C", "C2"), "120.1", false),
                (MeanCh("C1"), "1.111", false),
            ],
        ),
        Row {
            compare_sign: false,
            ..r(
                "methyl formate, O typed O_3",
                include_str!("hand_built/mfo.cif"),
                vec![
                    (Bond("C", "O1"), "1.219", false),
                    (Bond("C", "O2"), "1.401", false),
                    (Bond("O2", "CM"), "1.425", false),
                    (Angle("C", "O2", "CM"), "113.6", false),
                ],
            )
        },
        Row {
            retype: Some(("O2", UffType::OR)),
            ..r(
                "methyl formate, O typed O_R",
                include_str!("hand_built/mfo.cif"),
                vec![
                    (Bond("C", "O1"), "1.219", false),
                    (Bond("C", "O2"), "1.401", false),
                    (Bond("O2", "CM"), "1.425", false),
                    (Angle("C", "O2", "CM"), "113.6", false),
                ],
            )
        },
        r(
            "acetamide",
            include_str!("../components/ACM.cif"),
            vec![
                (Bond("C1", "O"), "1.222", false),
                (Bond("C1", "C2"), "1.498", false),
                (Bond("C1", "N"), "1.365", false),
                (Angle("C2", "C1", "N"), "120.3", true),
                (Bond("N", "HN1"), "1.045", true),
            ],
        ),
        r(
            "N-methylformamide",
            include_str!("hand_built/nmf.cif"),
            vec![
                (Bond("C", "O"), "1.218", false),
                (Bond("C", "N"), "1.365", false),
                (Bond("N", "CM"), "1.459", true),
                (Angle("O", "C", "N"), "120.1", false),
                (Angle("C", "N", "CM"), "122.1", false),
                (MeanCh("CM"), "1.111", false),
            ],
        ),
    ]
}

/// Every measure of `row` at its minimum under `variant`.
fn measure(row: &Row, variant: Variant) -> Vec<f64> {
    let c = Component::from_ccd(row.text).expect("parses");
    let mut types = uff::assign(&c);
    if let Some((name, t)) = row.retype {
        types[index(&c, name)] = t;
    }
    let ff = ForceField::with_variant(&c, &types, variant).expect("supported");
    let (_, at) = relaxed(&ff, &positions(&c));
    row.measures
        .iter()
        .map(|(m, _, _)| match *m {
            Measure::Bond(a, b) => distance(&at, index(&c, a), index(&c, b)),
            Measure::Angle(a, b, d) => angle(&at, index(&c, a), index(&c, b), index(&c, d)),
            Measure::MinCh(a) | Measure::MaxCh(a) => {
                let i = index(&c, a);
                let hs = c
                    .neighbours(i)
                    .filter(|(j, _)| c.atoms()[*j].element == Element::H)
                    .map(|(j, _)| distance(&at, i, j));
                if matches!(m, Measure::MinCh(_)) {
                    hs.fold(f64::INFINITY, f64::min)
                } else {
                    hs.fold(0.0, f64::max)
                }
            }
            Measure::MeanCh(a) => {
                let i = index(&c, a);
                let hs: Vec<f64> = c
                    .neighbours(i)
                    .filter(|(j, _)| c.atoms()[*j].element == Element::H)
                    .map(|(j, _)| distance(&at, i, j))
                    .collect();
                hs.iter().sum::<f64>() / hs.len() as f64
            }
        })
        .collect()
}

/// Whether a bond's natural length depends on the sign of `r_EN`: its two atoms are different
/// elements.
fn sign_matters(c: &Component, m: Measure) -> bool {
    match m {
        Measure::Bond(a, b) => c.atoms()[index(c, a)].element != c.atoms()[index(c, b)].element,
        Measure::MeanCh(_) | Measure::MinCh(_) | Measure::MaxCh(_) => true,
        Measure::Angle(..) => false,
    }
}

fn half_digit(printed: &str) -> f64 {
    let decimals = printed.split('.').nth(1).map_or(0, str::len) as i32;
    0.5 * 10f64.powi(-decimals)
}

/// The minimisation's own error on a length (Å) or angle (°): 1e-6 covers both — see the file's
/// documentation.
const NUMERICAL: f64 = 1e-6;

#[test]
fn figures_4_to_7_under_both_signs() {
    let plus = Variant {
        electronegativity: ElectronegativitySign::AddedAsPrinted,
        ..Variant::default()
    };
    let zero = Variant {
        electronegativity: ElectronegativitySign::Dropped,
        ..Variant::default()
    };
    println!("| molecule | measure | paper | −r_EN (this crate) | +r_EN (as printed) | no r_EN | −r_EN within ±½ digit |");
    println!("| --- | --- | --- | --- | --- | --- | --- |");
    let mut nearer = 0;
    let mut nearer_than_zero = 0;
    let mut zero_nearer_than_plus = 0;
    let mut sign_bonds = 0;
    let mut asserted = 0;
    for row in rows() {
        let c = Component::from_ccd(row.text).expect("parses");
        let minus = measure(&row, Variant::default());
        let added = measure(&row, plus);
        let none = measure(&row, zero);
        for (k, (m, printed, agrees)) in row.measures.iter().enumerate() {
            let paper: f64 = printed.parse().expect("a number");
            let within = (minus[k] - paper).abs() <= half_digit(printed) + NUMERICAL;
            println!(
                "| {} | {m:?} | {printed} | {:.4} | {:.4} | {:.4} | {} |",
                row.molecule,
                minus[k],
                added[k],
                none[k],
                if within { "yes" } else { "no" }
            );
            let a_c_h = matches!(
                m,
                Measure::MeanCh(_) | Measure::MinCh(_) | Measure::MaxCh(_)
            );
            if *agrees {
                assert!(!a_c_h, "a C–H mean or range is never asserted");
                assert!(
                    within,
                    "{} {m:?}: {} against {printed}",
                    row.molecule, minus[k]
                );
                asserted += 1;
            }
            if row.compare_sign && sign_matters(&c, *m) && !a_c_h {
                sign_bonds += 1;
                let (dm, dp, dz) = (
                    (minus[k] - paper).abs(),
                    (added[k] - paper).abs(),
                    (none[k] - paper).abs(),
                );
                if dm < dp {
                    nearer += 1;
                }
                if dm < dz {
                    nearer_than_zero += 1;
                }
                if dz < dp {
                    zero_nearer_than_plus += 1;
                }
            }
        }
    }
    println!("{asserted} measures asserted to agree");
    println!(
        "of {sign_bonds} heteronuclear bonds: −r_EN nearer than +r_EN for {nearer}, −r_EN nearer \
         than no r_EN for {nearer_than_zero}, no r_EN nearer than +r_EN for {zero_nearer_than_plus}"
    );
    assert_eq!(nearer, sign_bonds);
    assert_eq!(nearer_than_zero, sign_bonds);
    assert_eq!(zero_nearer_than_plus, sign_bonds);
}
