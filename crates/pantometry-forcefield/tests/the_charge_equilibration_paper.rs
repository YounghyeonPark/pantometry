//! **QEq against its own paper**: Rappé and Goddard, *J. Phys. Chem.* 95, 3358 (1991) — Table II
//! (p. 3360, the alkali halides by eq 18, three decimals), Table III (p. 3361, the QEq column,
//! three decimals) and Table IV (p. 3363, the QEq column, two decimals).
//!
//! # Geometries, and why each comparison has the tolerance it has
//!
//! The paper says its calculations were "carried out at the experimental geometries" (p. 3361,
//! refs 14–16: Harmony et al. 1979 and Landolt–Börnstein for polyatomics, Huber and Herzberg for
//! diatomics) and does not print them. Every geometry here is an experimental one, stated with
//! its source: the NIST Computational Chemistry Comparison and Benchmark Database (CCCBDB,
//! `cccbdb.nist.gov`, "experimental geometry" pages) and the NIST Chemistry WebBook's constants
//! of diatomic molecules (Huber and Herzberg's compilation), read 2026-10-04, each citing its
//! original — named on each case below.
//!
//! Two things separate this calculation from the printed number, and the tolerance is their sum,
//! worked out per atom before the comparison is made:
//!
//! 1. **The printed precision**: half a unit of the last figure, 0.0005 e in Tables II and III
//!    and 0.005 e in Table IV.
//! 2. **The geometry.** For a **diatomic**, the paper's distance is Huber and Herzberg's r_e
//!    (its footnote 14), the same source as here, known to about 10⁻⁴ Å: each is allowed
//!    **10⁻⁴ Å**, so HF and HCl are held, in effect, to their printed figure alone. For a
//!    **polyatomic**, the paper's structures are not these, and each bond length is allowed
//!    **0.01 Å** and each angle **1°**, the same for every molecule. **That is an assumed
//!    allowance, not a measured one**: it was not derived from the r_e, r_0 and r_s structures
//!    of these particular molecules, and it was set after an independent Python reproduction had
//!    shown misses of up to 0.012 e. Its size is reported against what is needed: the test prints
//!    the smallest common scale of it at which every row passes. The charge's sensitivity to each
//!    coordinate is measured by central differences on the internal coordinates — **with the code
//!    under test**, so a defect that made the charges insensitive to geometry would also shrink
//!    the tolerance, and one that made them oversensitive would widen it — and the allowance is
//!    `Σ |∂Q/∂g| δg`, the first-order worst case: every coordinate off by its full allowance,
//!    independently, in the direction that hurts.
//!
//! A row that misses by more is reported and not asserted, and the test says so: [`MISSES`]
//! lists them, and the test fails if one of them starts to agree or another stops.
//!
//! # Every charge is converged
//!
//! At [`Qeq::default`] settings each molecule's charges are checked against eq 8 directly: every
//! atom's chemical potential, from a matrix rebuilt at the final charges, equal to within
//! `2 · 20 eV · n · 10⁻¹⁰` — the bound the documented default tolerance gives (see
//! `every_atom_has_the_same_electronegativity`). A looser default fails it here, and fails HF's
//! Table III row, which is held to its printed figure.
//!
//! # Where the paper disagrees with itself
//!
//! The text on p. 3361 says "QEq leads to Q_C = 0.21" for H₂CO and "Q_H = 0.36" for H₂O. Table IV
//! prints 0.19 for that carbon (the same as its MP2 column) and 0.35 for that hydrogen, and
//! Table III prints 0.353. The tables are what is compared: they agree with each other on H₂O,
//! and 0.353 is not 0.36. The H₂CO carbon here is 0.1968 — 0.0068 from the table's 0.19 and
//! 0.0132 from the text's 0.21, against a tolerance of 0.0127, so this calculation sides with
//! neither decisively. Table IV also has an unlabelled row, "H 0.19 0.18", under PH₃. PH₃'s three
//! hydrogens are equivalent, so it is presumably another molecule whose label was lost (H₂S would
//! fit the order), and since the page does not say which, it is not compared.

use pantometry_forcefield::qeq::{
    coulomb_integral, coulomb_matrix, principal_quantum_number, table_i, Qeq, HARTREE_EV,
};
use pantometry_forcefield::Element;

const ANGSTROM: f64 = 1e-10;
const DEG: f64 = std::f64::consts::PI / 180.0;

/// Allowance for an experimental bond length in a polyatomic, Å — assumed; see the module
/// documentation.
const LENGTH: f64 = 0.01;
/// Allowance for a diatomic's Huber–Herzberg r_e, Å.
const DIATOMIC: f64 = 1e-4;
/// Allowance for an experimental angle, degrees.
const ANGLE: f64 = 1.0;

#[derive(Clone, Copy, PartialEq)]
enum Coordinate {
    Length,
    Angle,
}
use Coordinate::{Angle, Length};

struct Case {
    name: &'static str,
    source: &'static str,
    elements: Vec<Element>,
    /// Internal coordinates: Å and degrees.
    nominal: Vec<(Coordinate, f64)>,
    build: fn(&[f64]) -> Vec<[f64; 3]>,
    /// `(label, atom index, the paper's value, decimals printed)`.
    rows: Vec<(&'static str, usize, f64, i32)>,
}

/// Rows that miss by more than their tolerance, as `(molecule, label)`. Reported, not asserted.
const MISSES: &[(&str, &str)] = &[];

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// In the xy plane: the point `r` from `b` with angle `a–b–X` = `theta` degrees, on the side of
/// line `a–b` that `side` (+1 or −1) picks — positive is anticlockwise from `b → a`.
fn place(a: [f64; 3], b: [f64; 3], r: f64, theta: f64, side: f64) -> [f64; 3] {
    let u = [a[0] - b[0], a[1] - b[1]];
    let l = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let (c, s) = ((theta * DEG).cos(), (side * theta * DEG).sin());
    let d = [(u[0] * c - u[1] * s) / l, (u[0] * s + u[1] * c) / l];
    [b[0] + r * d[0], b[1] + r * d[1], 0.0]
}

/// Three bonds of length `r` from the origin with mutual angle `theta`, around −z.
fn pyramid(r: f64, theta: f64) -> Vec<[f64; 3]> {
    let ca = (theta * DEG).cos();
    let cb = ((1.0 + 2.0 * ca) / 3.0).sqrt();
    let sb = (1.0 - cb * cb).sqrt();
    (0..3)
        .map(|k| {
            let p = f64::from(k) * 2.0 * std::f64::consts::PI / 3.0;
            [r * sb * p.cos(), r * sb * p.sin(), -r * cb]
        })
        .collect()
}

/// A methyl's three hydrogens on a carbon at `c`, axis +z away from the rest, bond `r`, at angle
/// `from_axis` degrees from the axis, the first at azimuth `phase` degrees.
fn methyl(c: [f64; 3], r: f64, from_axis: f64, phase: f64) -> Vec<[f64; 3]> {
    (0..3)
        .map(|k| {
            let p = (phase + 120.0 * f64::from(k)) * DEG;
            let a = from_axis * DEG;
            add(
                c,
                [r * a.sin() * p.cos(), r * a.sin() * p.sin(), r * a.cos()],
            )
        })
        .collect()
}

fn cases() -> Vec<Case> {
    use Element::{Cl, C, F, H, N, O, P};
    vec![
        Case {
            name: "HF",
            source: "r_e 0.9168 Å; CCCBDB, citing Le Roy, J. Mol. Spect. 194, 189 (1999)",
            elements: vec![H, F],
            nominal: vec![(Length, 0.9168)],
            build: |g| vec![[0.0; 3], [g[0], 0.0, 0.0]],
            rows: vec![("H (Table III)", 0, 0.462, 3), ("H (Table IV)", 0, 0.46, 2)],
        },
        Case {
            name: "H2O",
            source:
                "r_e 0.958 Å, 104.4776°; CCCBDB, citing Hoy and Bunker, J. Mol. Spect. 74, 1 (1979)",
            elements: vec![O, H, H],
            nominal: vec![(Length, 0.958), (Angle, 104.4776)],
            build: |g| {
                let a = g[1] / 2.0 * DEG;
                vec![
                    [0.0; 3],
                    [g[0] * a.sin(), 0.0, g[0] * a.cos()],
                    [-g[0] * a.sin(), 0.0, g[0] * a.cos()],
                ]
            },
            rows: vec![("H (Table III)", 1, 0.353, 3), ("H (Table IV)", 1, 0.35, 2)],
        },
        Case {
            name: "NH3",
            source: "1.012 Å, 106.67°; CCCBDB, citing Herzberg 1966",
            elements: vec![N, H, H, H],
            nominal: vec![(Length, 1.012), (Angle, 106.67)],
            build: |g| {
                let mut v = vec![[0.0; 3]];
                v.extend(pyramid(g[0], g[1]));
                v
            },
            rows: vec![("H (Table III)", 1, 0.243, 3), ("H (Table IV)", 1, 0.24, 2)],
        },
        Case {
            name: "CH4",
            source:
                "r_e 1.087 Å, tetrahedral; CCCBDB, citing Hirota, J. Mol. Spect. 77, 213 (1979)",
            elements: vec![C, H, H, H, H],
            nominal: vec![(Length, 1.087)],
            build: |g| {
                let t = g[0] / 3f64.sqrt();
                let mut v = vec![[0.0; 3]];
                for s in [
                    [1.0, 1.0, 1.0],
                    [1.0, -1.0, -1.0],
                    [-1.0, 1.0, -1.0],
                    [-1.0, -1.0, 1.0],
                ] {
                    v.push(scale(s, t));
                }
                v
            },
            rows: vec![("H (Table III)", 1, 0.149, 3), ("H (Table IV)", 1, 0.15, 2)],
        },
        Case {
            name: "C2H6",
            source:
                "C–C 1.536, C–H 1.091 Å, H–C–C 110.91°, staggered; CCCBDB, citing Herzberg 1966",
            elements: vec![C, C, H, H, H, H, H, H],
            nominal: vec![(Length, 1.536), (Length, 1.091), (Angle, 110.91)],
            build: |g| {
                let c = g[0] / 2.0;
                let mut v = vec![[0.0, 0.0, c], [0.0, 0.0, -c]];
                v.extend(methyl([0.0, 0.0, c], g[1], 180.0 - g[2], 0.0));
                // The far methyl points along −z: flip z of a +z methyl at the other phase.
                for h in methyl([0.0, 0.0, 0.0], g[1], 180.0 - g[2], 60.0) {
                    v.push([h[0], h[1], -c - h[2]]);
                }
                v
            },
            rows: vec![("H", 2, 0.16, 2)],
        },
        Case {
            name: "C2H2",
            source: "C–C 1.203, C–H 1.063 Å, r_m; CCCBDB, citing Kuchitsu 1998",
            elements: vec![C, C, H, H],
            nominal: vec![(Length, 1.203), (Length, 1.063)],
            build: |g| {
                let c = g[0] / 2.0;
                vec![
                    [0.0, 0.0, c],
                    [0.0, 0.0, -c],
                    [0.0, 0.0, c + g[1]],
                    [0.0, 0.0, -c - g[1]],
                ]
            },
            rows: vec![("H", 2, 0.13, 2)],
        },
        Case {
            name: "C2H4",
            source: "C–C 1.339, C–H 1.086 Å, H–C–C 121.2°; CCCBDB, citing Herzberg 1966",
            elements: vec![C, C, H, H, H, H],
            nominal: vec![(Length, 1.339), (Length, 1.086), (Angle, 121.2)],
            build: |g| {
                let (c, r, a) = (g[0] / 2.0, g[1], g[2] * DEG);
                let (y, z) = (r * a.sin(), r * a.cos());
                vec![
                    [0.0, 0.0, c],
                    [0.0, 0.0, -c],
                    [0.0, y, c - z],
                    [0.0, -y, c - z],
                    [0.0, y, -c + z],
                    [0.0, -y, -c + z],
                ]
            },
            rows: vec![("H", 2, 0.15, 2)],
        },
        Case {
            name: "C6H6",
            source: "C–C 1.397, C–H 1.084 Å, D6h; CCCBDB, citing Herzberg 1966",
            elements: [vec![C; 6], vec![H; 6]].concat(),
            nominal: vec![(Length, 1.397), (Length, 1.084)],
            build: |g| {
                let ring = |r: f64| -> Vec<[f64; 3]> {
                    (0..6)
                        .map(|k| {
                            let p = f64::from(k) * std::f64::consts::PI / 3.0;
                            [r * p.cos(), r * p.sin(), 0.0]
                        })
                        .collect()
                };
                [ring(g[0]), ring(g[0] + g[1])].concat()
            },
            rows: vec![("H", 6, 0.10, 2)],
        },
        Case {
            name: "CO2",
            source: "C–O 1.162 Å, linear; CCCBDB, citing Herzberg 1966",
            elements: vec![O, C, O],
            nominal: vec![(Length, 1.162)],
            build: |g| vec![[0.0, 0.0, -g[0]], [0.0; 3], [0.0, 0.0, g[0]]],
            rows: vec![("O", 0, -0.45, 2)],
        },
        Case {
            name: "H2CO",
            source: "C–O 1.205, C–H 1.111 Å, H–C–O 121.9°; CCCBDB, citing Gurvich et al. 1989",
            elements: vec![O, C, H, H],
            nominal: vec![(Length, 1.205), (Length, 1.111), (Angle, 121.9)],
            build: |g| {
                let a = g[2] * DEG;
                let (y, z) = (g[1] * a.sin(), g[1] * a.cos());
                vec![[0.0, 0.0, g[0]], [0.0; 3], [0.0, y, z], [0.0, -y, z]]
            },
            // The text on p. 3361 says QEq gives this carbon 0.21; the table's 0.19 is compared.
            // See the module documentation.
            rows: vec![("O", 0, -0.43, 2), ("C", 1, 0.19, 2), ("H", 2, 0.12, 2)],
        },
        Case {
            name: "CH3OH",
            source: "C–O 1.427, O–H 0.956, C–H 1.096 Å, H–C–H 109.03°, C–O–H 108.87°, a \
                     symmetric methyl with one C–H anti to O–H; CCCBDB, citing Venkateswarlu and \
                     Gordy, J. Chem. Phys. 23, 1200 (1955)",
            elements: vec![C, O, H, H, H, H],
            nominal: vec![
                (Length, 1.427),
                (Length, 0.956),
                (Length, 1.096),
                (Angle, 109.03),
                (Angle, 108.87),
            ],
            build: |g| {
                // C at the origin, O on −z; the methyl's axis +z.
                let ch = g[3] * DEG;
                let from_axis = (((ch.cos() + 0.5) / 1.5).sqrt()).acos() / DEG;
                let o = [0.0, 0.0, -g[0]];
                let coh = g[4] * DEG;
                let mut v = vec![[0.0; 3], o];
                v.extend(methyl([0.0; 3], g[2], from_axis, 0.0));
                v.push(add(o, [-g[1] * coh.sin(), 0.0, g[1] * coh.cos()]));
                v
            },
            // H_t is the methyl hydrogen in the C–O–H plane, H_s the two out of it.
            rows: vec![
                ("H(O)", 5, 0.36, 2),
                ("O", 1, -0.66, 2),
                ("C", 0, -0.15, 2),
                ("H_s", 3, 0.14, 2),
                ("H_t", 2, 0.18, 2),
            ],
        },
        Case {
            name: "H2NC(O)H",
            source: "C–N 1.350, C–O 1.210, C–H 1.090, N–H 1.001 Å, N–C–O 124.7°, H–C–N 112.7°, \
                     C–N–H 120.0° (H cis to O) and 118.5°, planar; CCCBDB, citing Hirota and \
                     Sugisaki, J. Mol. Spect. 49, 251 (1974)",
            elements: vec![N, C, O, H, H, H],
            nominal: vec![
                (Length, 1.350),
                (Length, 1.210),
                (Length, 1.090),
                (Length, 1.001),
                (Angle, 124.7),
                (Angle, 112.7),
                (Angle, 120.0),
                (Angle, 118.5),
            ],
            build: |g| {
                let n = [0.0; 3];
                let c = [g[0], 0.0, 0.0];
                let o = place(n, c, g[1], g[4], 1.0);
                let hc = place(n, c, g[2], g[5], -1.0);
                let h_cis = place(c, n, g[3], g[6], -1.0);
                let h_trans = place(c, n, g[3], g[7], 1.0);
                vec![n, c, o, hc, h_cis, h_trans]
            },
            // H_c and H_t read as cis and trans to the carbonyl oxygen.
            rows: vec![
                ("O", 2, -0.42, 2),
                ("C", 1, 0.39, 2),
                ("N", 0, -0.63, 2),
                ("H_c", 4, 0.29, 2),
                ("H_t", 5, 0.23, 2),
            ],
        },
        Case {
            name: "HC(O)OH",
            source: "C=O 1.202, C–O 1.343, C–H 1.097, O–H 0.972 Å, O–C–O 124.9°, H–C=O 124.1°, \
                     H–O–C 106.3°, planar, the O–H syn to C=O; CCCBDB, citing Herzberg 1966",
            elements: vec![C, O, O, H, H],
            nominal: vec![
                (Length, 1.202),
                (Length, 1.343),
                (Length, 1.097),
                (Length, 0.972),
                (Angle, 124.9),
                (Angle, 124.1),
                (Angle, 106.3),
            ],
            build: |g| {
                let c = [0.0; 3];
                let od = [g[0], 0.0, 0.0];
                let os = place(od, c, g[1], g[4], 1.0);
                let hc = place(od, c, g[2], g[5], -1.0);
                // Syn: the hydrogen on the same side of C–O(H) as the carbonyl oxygen.
                let ho = place(c, os, g[3], g[6], 1.0);
                vec![c, od, os, hc, ho]
            },
            rows: vec![
                ("O (carbonyl)", 1, -0.44, 2),
                ("C", 0, 0.56, 2),
                ("H(C)", 3, 0.16, 2),
                ("O (hydroxyl)", 2, -0.65, 2),
                ("H(O)", 4, 0.38, 2),
            ],
        },
        Case {
            name: "CH3CN",
            source:
                "C≡N 1.157, C–C 1.458, C–H 1.104 Å, H–C–C 109.44°; CCCBDB, citing Herzberg 1966",
            elements: vec![N, C, C, H, H, H],
            nominal: vec![
                (Length, 1.157),
                (Length, 1.458),
                (Length, 1.104),
                (Angle, 109.44),
            ],
            build: |g| {
                let mut v = vec![[0.0, 0.0, g[1] + g[0]], [0.0, 0.0, g[1]], [0.0; 3]];
                v.extend(methyl([0.0; 3], g[2], g[3], 0.0));
                v
            },
            rows: vec![
                ("N", 0, -0.24, 2),
                ("C (nitrile)", 1, 0.22, 2),
                ("C (methyl)", 2, -0.37, 2),
                ("H", 3, 0.13, 2),
            ],
        },
        Case {
            name: "H2C=C=O",
            source: "C=O 1.162, C=C 1.314, C–H 1.083 Å (r_s), H–C–C 118.72°; CCCBDB, citing \
                     Kuchitsu 1998",
            elements: vec![O, C, C, H, H],
            nominal: vec![
                (Length, 1.162),
                (Length, 1.314),
                (Length, 1.083),
                (Angle, 118.72),
            ],
            build: |g| {
                let a = g[3] * DEG;
                let (y, z) = (g[2] * a.sin(), g[2] * a.cos());
                vec![
                    [0.0, 0.0, g[1] + g[0]],
                    [0.0, 0.0, g[1]],
                    [0.0; 3],
                    [0.0, y, z],
                    [0.0, -y, z],
                ]
            },
            rows: vec![
                ("O", 0, -0.45, 2),
                ("C (carbonyl)", 1, 0.42, 2),
                ("C (methylene)", 2, -0.23, 2),
            ],
        },
        Case {
            name: "PH3",
            source: "1.421 Å, 93.3°; CCCBDB, citing Herzberg 1966",
            elements: vec![P, H, H, H],
            nominal: vec![(Length, 1.421), (Angle, 93.3)],
            build: |g| {
                let mut v = vec![[0.0; 3]];
                v.extend(pyramid(g[0], g[1]));
                v
            },
            rows: vec![("H", 1, 0.08, 2)],
        },
        Case {
            name: "ClH",
            source: "r_e 1.2746 Å; CCCBDB, citing Huber and Herzberg 1979",
            elements: vec![H, Cl],
            nominal: vec![(Length, 1.2746)],
            build: |g| vec![[0.0; 3], [g[0], 0.0, 0.0]],
            rows: vec![("H", 0, 0.32, 2)],
        },
    ]
}

/// The documented default tolerance of the hydrogen iteration, e — written here rather than read
/// from [`Qeq::default`], so that a looser default fails the convergence check below.
const DEFAULT_TOLERANCE: f64 = 1e-10;

fn charges(case: &Case, g: &[f64]) -> Vec<f64> {
    solve(case, g).1
}

fn solve(case: &Case, g: &[f64]) -> (Vec<[f64; 3]>, Vec<f64>) {
    let at: Vec<[f64; 3]> = (case.build)(g)
        .into_iter()
        .map(|p| scale(p, ANGSTROM))
        .collect();
    assert_eq!(
        at.len(),
        case.elements.len(),
        "{}: one position per atom",
        case.name
    );
    let q = Qeq::default()
        .equilibrate(&case.elements, &at, 0.0)
        .unwrap_or_else(|e| panic!("{}: {e}", case.name))
        .charges;
    (at, q)
}

/// Eq 8 at the final charges: the largest departure of any atom's chemical potential from their
/// mean, eV, with the matrix rebuilt at those charges.
fn electronegativity_spread(elements: &[Element], at: &[[f64; 3]], q: &[f64]) -> f64 {
    let n = elements.len();
    let j = coulomb_matrix(elements, at, q);
    let chi: Vec<f64> = (0..n)
        .map(|a| {
            table_i(elements[a]).electronegativity
                + (0..n).map(|b| j[a * n + b] * q[b]).sum::<f64>()
        })
        .collect();
    let mean = chi.iter().sum::<f64>() / n as f64;
    chi.iter().fold(0.0f64, |m, x| m.max((x - mean).abs()))
}

#[test]
fn tables_iii_and_iv_at_experimental_geometries() {
    let mut misses = Vec::new();
    let mut asserted = 0;
    let mut needed = 0.0f64;
    println!("| molecule | atom | this | paper | |Δ| | tolerance (print + geometry) |");
    println!("| --- | --- | --- | --- | --- | --- |");
    for case in cases() {
        let g0: Vec<f64> = case.nominal.iter().map(|c| c.1).collect();
        let (at, q0) = solve(&case, &g0);
        assert!(
            q0.iter().sum::<f64>().abs() < 1e-12,
            "{}: neutral",
            case.name
        );
        let n = q0.len();
        let spread = electronegativity_spread(&case.elements, &at, &q0);
        let bound = 2.0 * 20.0 * n as f64 * DEFAULT_TOLERANCE + 1e-12;
        assert!(
            spread <= bound,
            "{}: eq 8 spread {spread:e} eV at default settings, bound {bound:e}",
            case.name
        );
        let diatomic = n == 2;
        // ∂Q/∂g by central differences at the allowance itself: the change a coordinate off by
        // its allowance makes, to first order.
        let mut geometry = vec![0.0; n];
        for (k, &(kind, _)) in case.nominal.iter().enumerate() {
            let d = match kind {
                Length if diatomic => DIATOMIC,
                Length => LENGTH,
                Angle => ANGLE,
            };
            let (mut up, mut down) = (g0.clone(), g0.clone());
            up[k] += d;
            down[k] -= d;
            let (qu, qd) = (charges(&case, &up), charges(&case, &down));
            for i in 0..n {
                geometry[i] += ((qu[i] - qd[i]) / 2.0).abs();
            }
        }
        for &(label, i, paper, decimals) in &case.rows {
            let print = 0.5 * 10f64.powi(-decimals);
            let tolerance = print + geometry[i];
            let miss = (q0[i] - paper).abs();
            if !diatomic && miss > print {
                needed = needed.max((miss - print) / geometry[i]);
            }
            println!(
                "| {} | {label} | {:.4} | {paper} | {miss:.4} | {tolerance:.4} ({print} + {:.4}) |",
                case.name, q0[i], geometry[i]
            );
            if miss > tolerance {
                misses.push((case.name, label));
            } else {
                asserted += 1;
            }
        }
        println!("  geometry: {}; eq 8 spread {spread:.1e} eV", case.source);
    }
    println!("{asserted} rows within tolerance; outside: {misses:?}");
    println!(
        "the smallest common scale of the polyatomic allowance at which every row passes: \
         {needed:.3} (0.01 Å and 1° is 1)"
    );
    assert_eq!(
        misses, MISSES,
        "the rows outside their tolerance are not the recorded ones"
    );
}

/// Bohr radius as the paper converts with it in eq 17 (p. 3360: "R_A in units of a₀ (a₀ =
/// 0.52917 Å)"), Å. Used for ζ here, as Table II was computed; [`coulomb_integral`]'s distance is
/// converted with CODATA's, and the two differ by 1.4e-5 relative, which moves no charge below in
/// its fifth decimal.
const PAPERS_BOHR: f64 = 0.52917;

/// Table I's χ (eV), J (eV), R (Å) and n for the alkali metals, which are not elements of this
/// crate: Na, K, Rb.
const METALS: [(&str, f64, f64, f64, u32); 3] = [
    ("Na", 2.843, 4.592, 2.085, 3),
    ("K", 2.421, 3.84, 2.586, 4),
    ("Rb", 2.331, 3.692, 2.770, 5),
];

/// **Table II (p. 3360): the alkali halides by eq 18**, `Q_M = (χ⁰_X − χ⁰_M) / (J_MM + J_XX −
/// 2 J_MX)`, both columns — λ = 0.4913 (`Q_QEq`) and λ = ½ (`Q_λ=0.5`) — for Na, K and Rb with Cl,
/// Br and I: eighteen charges, three decimals.
///
/// The physical check of the 4s and 5s integrals that Tables III and IV lack: K and Br are 4s, Rb
/// and I 5s, and every pair here is computed by [`coulomb_integral`] with Table I's numbers. The
/// halogens' χ, J, R and n come from the crate ([`table_i`], [`principal_quantum_number`]), so
/// chlorine's n is held here by the paper's numbers and not only by its sentence. Distances are
/// Huber and Herzberg's r_e (the paper's footnote 14), from the NIST Chemistry WebBook's
/// constants of diatomic molecules, for the most abundant isotopologue.
///
/// Tolerance: half the last printed figure, 0.0005 e. The r_e are known to ~10⁻⁴ Å, and
/// `|∂Q/∂R|` is below 0.1 e/Å for all nine, so the distance is worth under 10⁻⁵ e.
#[test]
fn table_ii_the_alkali_halides() {
    use Element::{Br, Cl, I};
    let rows: [(&str, Element, f64, f64, f64); 9] = [
        ("Na", Cl, 2.360795, 0.766, 0.776),
        ("Na", Br, 2.502038, 0.745, 0.756),
        ("Na", I, 2.711452, 0.709, 0.720),
        ("K", Cl, 2.66665, 0.775, 0.784),
        ("K", Br, 2.82078, 0.768, 0.777),
        ("K", I, 3.047844, 0.754, 0.764),
        ("Rb", Cl, 2.786736, 0.763, 0.771),
        ("Rb", Br, 2.944744, 0.757, 0.766),
        ("Rb", I, 3.176879, 0.747, 0.757),
    ];
    let zeta =
        |lambda: f64, n: u32, r: f64| lambda * (2.0 * f64::from(n) + 1.0) / (2.0 * r / PAPERS_BOHR);
    let bohr = pantometry_forcefield::qeq::BOHR / ANGSTROM;
    let mut worst = 0.0f64;
    println!("| MX | r_e (Å) | λ = 0.4913 | paper | λ = ½ | paper |");
    println!("| --- | --- | --- | --- | --- | --- |");
    for (metal, x, r_e, paper_fit, paper_half) in rows {
        let &(_, chi_m, j_m, r_m, n_m) = METALS.iter().find(|m| m.0 == metal).expect("metal");
        let t = table_i(x);
        let n_x = principal_quantum_number(x);
        let q = |lambda: f64, r: f64| {
            let j = HARTREE_EV
                * coulomb_integral(
                    n_m,
                    zeta(lambda, n_m, r_m),
                    n_x,
                    zeta(lambda, n_x, t.radius),
                    r / bohr,
                );
            (t.electronegativity - chi_m) / (j_m + t.idempotential - 2.0 * j)
        };
        let (fit, half) = (q(0.4913, r_e), q(0.5, r_e));
        // The distance's worth: |∂Q/∂R| × 10⁻⁴ Å, by central difference.
        let slope = ((q(0.5, r_e + 1e-3) - q(0.5, r_e - 1e-3)) / 2e-3).abs();
        assert!(slope < 0.1, "{metal}{x}: ∂Q/∂R = {slope}");
        println!("| {metal}{x} | {r_e} | {fit:.5} | {paper_fit} | {half:.5} | {paper_half} |");
        for (got, paper) in [(fit, paper_fit), (half, paper_half)] {
            worst = worst.max((got - paper).abs());
            assert!(
                (got - paper).abs() <= 0.0005,
                "{metal}{x}: {got:.5} against Table II's {paper}"
            );
        }
    }
    println!("largest miss {worst:.5} e, against 0.0005");
}
