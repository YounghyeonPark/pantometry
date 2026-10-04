//! **The congener series: nine ligands in T4 lysozyme L99A's cavity, built with the pipeline that
//! built 181L, and their rigid-pocket binding energies set beside experiment.** Step 2c-3a.
//!
//! The entries are the nine of A. Morton and B. W. Matthews, *Biochemistry* **34**, 8576 (1995):
//! 181L benzene, 182L benzofuran, 183L indene, 184L isobutylbenzene, 185L indole, 186L
//! n-butylbenzene, 187L p-xylene, 188L o-xylene and 1NHB ethylbenzene. Each is held, by string
//! operations on the file and not by the crate, to being the same protein as 181L: the same
//! `SEQRES`, the same three `SEQADV` conflicts (C54T, C97A, L99A), the same unmodelled residues,
//! and nothing besides waters, two chlorides, one 2-hydroxyethyl disulfide and its ligand. Each
//! ligand's dictionary entry is held to the molecule by its bond graph, written here by hand from
//! the chemistry. Then every entry is built, and the rule makes His31 +1 in each.
//!
//! **What is asserted about the binding energies is only what must hold**: van der Waals
//! attraction, buried area, frozen atoms, convergence, and 181L's numbers as 2c-2 committed them.
//! **No correlation with experiment is asserted.** The series is nine ligands of one size whose
//! measured ΔG spans 2.1 kcal/mol; a correlation measured on nine points has a 95% interval that
//! the series test computes and prints, and it is wide.
//!
//! The series takes minutes, so it is ignored by default:
//! `cargo test -p pantometry-forcefield --release --test the_congener_series -- --ignored`.

mod protein;

use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    Binding, Component, Element, Histidine, Minimiser, Part, Selection, Status, System,
};
use protein::*;
use std::collections::BTreeMap;

const ANGSTROM: f64 = 1e-10;

/// The minimiser's tolerance, 2c-2's: the largest force on a ligand atom at most
/// 1e-4 kcal mol⁻¹ Å⁻¹.
const TOLERANCE: f64 = 1e-4 * KCAL_PER_MOL_ANGSTROM;

/// Points per atom for the surface, 2c-2's.
const POINTS: usize = 1600;

/// One entry of the series.
struct Entry {
    /// The PDB code.
    code: &'static str,
    /// The ligand's dictionary code, as the entry's `HETNAM` gives it.
    ligand: &'static str,
    /// The molecule, as Mobley et al. name it.
    name: &'static str,
    /// The entry as fetched.
    pdb: &'static str,
    /// The ligand's dictionary entry as fetched.
    ccd: &'static str,
    /// ΔG°, kcal/mol, and its uncertainty: Mobley, Graves, Chodera, McReynolds, Shoichet and
    /// Dill, *J. Mol. Biol.* **371**, 1118 (2007), Table 1, read at PMC2104542, whose caption
    /// gives them as from its ref. 30: Morton, Baase and Matthews, *Biochemistry* **34**, 8564
    /// (1995), by isothermal titration calorimetry at 302 K. That paper itself was not read.
    dg: (f64, f64),
}

const ENTRIES: [Entry; 9] = [
    Entry {
        code: "181L",
        ligand: "BNZ",
        name: "benzene",
        pdb: include_str!("../components/181L.pdb"),
        ccd: include_str!("../components/BNZ.cif"),
        dg: (-5.19, 0.16),
    },
    Entry {
        code: "182L",
        ligand: "BZF",
        name: "2,3-benzofuran",
        pdb: include_str!("../components/182L.pdb"),
        ccd: include_str!("../components/BZF.cif"),
        dg: (-5.46, 0.03),
    },
    Entry {
        code: "183L",
        ligand: "DEN",
        name: "indene",
        pdb: include_str!("../components/183L.pdb"),
        ccd: include_str!("../components/DEN.cif"),
        dg: (-5.13, 0.01),
    },
    Entry {
        code: "184L",
        ligand: "I4B",
        name: "isobutylbenzene",
        pdb: include_str!("../components/184L.pdb"),
        ccd: include_str!("../components/I4B.cif"),
        dg: (-6.51, 0.06),
    },
    Entry {
        code: "185L",
        ligand: "IND",
        name: "indole",
        pdb: include_str!("../components/185L.pdb"),
        ccd: include_str!("../components/IND.cif"),
        dg: (-4.89, 0.06),
    },
    Entry {
        code: "186L",
        ligand: "N4B",
        name: "n-butylbenzene",
        pdb: include_str!("../components/186L.pdb"),
        ccd: include_str!("../components/N4B.cif"),
        dg: (-6.70, 0.02),
    },
    Entry {
        code: "187L",
        ligand: "PXY",
        name: "p-xylene",
        pdb: include_str!("../components/187L.pdb"),
        ccd: include_str!("../components/PXY.cif"),
        dg: (-4.67, 0.06),
    },
    Entry {
        code: "188L",
        ligand: "OXE",
        name: "o-xylene",
        pdb: include_str!("../components/188L.pdb"),
        ccd: include_str!("../components/OXE.cif"),
        dg: (-4.60, 0.06),
    },
    Entry {
        code: "1NHB",
        ligand: "PYJ",
        name: "ethylbenzene",
        pdb: include_str!("../components/1NHB.pdb"),
        ccd: include_str!("../components/PYJ.cif"),
        dg: (-5.76, 0.07),
    },
];

/// What every entry drops: its waters, its two chloride ions and its 2-hydroxyethyl disulfide.
const DROPPED: [&str; 3] = ["HOH", "CL", "HED"];

fn kcal(joules: f64) -> f64 {
    joules / KCAL_PER_MOL
}

fn ligand_template(e: &Entry) -> Component {
    Component::from_ccd(e.ccd).unwrap_or_else(|err| panic!("{}: {err}", e.ligand))
}

/// The twenty amino acids and this entry's ligand.
fn templates_for(e: &Entry) -> Vec<Component> {
    let mut t = templates();
    if e.ligand != "BNZ" {
        t.push(ligand_template(e));
    }
    t
}

fn build(e: &Entry) -> System {
    System::from_pdb(
        e.pdb,
        &templates_for(e),
        &Selection::new(e.ligand).dropping(&DROPPED),
    )
    .unwrap_or_else(|err| panic!("{}: {err}", e.code))
}

/// The `SEQADV` records of `text` as `(residue, number, database residue, comment)`, read with
/// string operations.
fn seqadv(text: &str) -> Vec<(String, i32, String, String)> {
    text.lines()
        .filter(|l| l.starts_with("SEQADV"))
        .map(|l| {
            (
                columns(l, 13, 15).trim().to_string(),
                columns(l, 19, 22).trim().parse().expect("a residue number"),
                columns(l, 40, 42).trim().to_string(),
                columns(l, 50, 70).trim().to_string(),
            )
        })
        .collect()
}

/// A heavy atom's place in its molecule's graph: element, aromatic flag, heavy neighbours and
/// hydrogens.
type Site = (Element, bool, usize, usize);

/// Every heavy atom's [`Site`], sorted.
fn fingerprint(c: &Component) -> Vec<Site> {
    let atoms = c.atoms();
    let mut out: Vec<Site> = (0..atoms.len())
        .filter(|&i| atoms[i].element != Element::H)
        .map(|i| {
            let (h, heavy) = c
                .neighbours(i)
                .partition::<Vec<_>, _>(|(j, _)| atoms[*j].element == Element::H);
            (atoms[i].element, atoms[i].aromatic, heavy.len(), h.len())
        })
        .collect();
    out.sort();
    out
}

/// The molecule each code is, written by hand from its structure: `n` copies of each site.
fn expected_fingerprint(code: &str) -> Vec<Site> {
    use Element::{C, N, O};
    let ring_ch = (C, true, 2, 1);
    let ring_c = (C, true, 3, 0);
    let ch3 = (C, false, 1, 3);
    let ch2 = (C, false, 2, 2);
    let spec: Vec<(Site, usize)> = match code {
        "BNZ" => vec![(ring_ch, 6)],
        // Furan's O and its two CH, fused to benzene at two carbons.
        "BZF" => vec![(ring_ch, 6), (ring_c, 2), ((O, true, 2, 0), 1)],
        // Benzene fused to a five-ring of CH2–CH=CH, which the dictionary does not flag aromatic.
        "DEN" => vec![(ring_ch, 4), (ring_c, 2), (ch2, 1), ((C, false, 2, 1), 2)],
        // Ph–CH2–CH(CH3)2.
        "I4B" => vec![
            (ring_ch, 5),
            (ring_c, 1),
            (ch2, 1),
            ((C, false, 3, 1), 1),
            (ch3, 2),
        ],
        // Pyrrole's NH and two CH, fused to benzene.
        "IND" => vec![(ring_ch, 6), (ring_c, 2), ((N, true, 2, 1), 1)],
        // Ph–CH2–CH2–CH2–CH3.
        "N4B" => vec![(ring_ch, 5), (ring_c, 1), (ch2, 3), (ch3, 1)],
        // Two methyls on a ring; which two carbons is the graph distance, checked separately.
        "PXY" | "OXE" => vec![(ring_ch, 4), (ring_c, 2), (ch3, 2)],
        // Ph–CH2–CH3.
        "PYJ" => vec![(ring_ch, 5), (ring_c, 1), (ch2, 1), (ch3, 1)],
        other => panic!("no structure written for {other}"),
    };
    let mut out: Vec<Site> = spec
        .into_iter()
        .flat_map(|(s, n)| (0..n).map(move |_| s))
        .collect();
    out.sort();
    out
}

/// Bonds between the two methyl carbons, by breadth-first search over heavy atoms: 3 for ortho, 5
/// for para.
fn methyl_separation(c: &Component) -> usize {
    let atoms = c.atoms();
    let methyls: Vec<usize> = (0..atoms.len())
        .filter(|&i| {
            atoms[i].element == Element::C
                && c.neighbours(i)
                    .filter(|(j, _)| atoms[*j].element == Element::H)
                    .count()
                    == 3
        })
        .collect();
    assert_eq!(methyls.len(), 2, "{}", c.id());
    let mut depth = vec![usize::MAX; atoms.len()];
    depth[methyls[0]] = 0;
    let mut queue = std::collections::VecDeque::from([methyls[0]]);
    while let Some(i) = queue.pop_front() {
        for (j, _) in c.neighbours(i) {
            if atoms[j].element != Element::H && depth[j] == usize::MAX {
                depth[j] = depth[i] + 1;
                queue.push_back(j);
            }
        }
    }
    depth[methyls[1]]
}

/// **Each ligand's dictionary entry is the molecule the entry names**: its stated formula, its
/// name, and its bond graph — every heavy atom's element, aromatic flag, heavy neighbours and
/// hydrogens, against the structure written by hand — and for the two xylenes the methyls 3 bonds
/// apart (ortho) or 5 (para), which a formula and a site count cannot tell apart.
#[test]
fn each_ligand_is_the_molecule_its_entry_names() {
    for e in &ENTRIES {
        let c = ligand_template(e);
        assert_eq!(c.id(), e.ligand);
        assert_eq!(
            fingerprint(&c),
            expected_fingerprint(e.ligand),
            "{}",
            e.ligand
        );
        let name = e
            .ccd
            .lines()
            .find(|l| l.starts_with("_chem_comp.name"))
            .expect("a name")
            .trim_start_matches("_chem_comp.name")
            .trim()
            .to_string();
        match e.ligand {
            "PXY" => assert_eq!(methyl_separation(&c), 5),
            "OXE" => assert_eq!(methyl_separation(&c), 3),
            _ => {}
        }
        eprintln!(
            "{} {}: {name}, {}, {} heavy atoms",
            e.code,
            e.ligand,
            c.formula().unwrap_or("?"),
            fingerprint(&c).len()
        );
    }
}

/// The atom named `name` of residue `number` in `s`, metres.
fn atom_at(s: &System, number: i32, name: &str) -> [f64; 3] {
    let r = s
        .residues()
        .iter()
        .find(|r| r.part == Part::Protein && r.number == number)
        .unwrap_or_else(|| panic!("no residue {number}"));
    let i = r
        .atoms
        .clone()
        .find(|&i| s.atom_name(i) == name)
        .unwrap_or_else(|| panic!("{} has no {name}", r.label()));
    s.component().atoms()[i].at
}

/// **Every entry is 181L's protein, and builds as 181L does.** Read from the file with string
/// operations: the `HEADER`'s code, `SEQRES` identical to 181L's, the same three `SEQADV`
/// conflicts against UniProt P00720 — Thr54 for Cys, Ala97 for Cys, Ala99 for Leu — and `HETATM`
/// residues that are the waters, two chlorides, one HED and the ligand only, each dropped and
/// counted as the file has them. Built: the same protein residues as 181L, residues 163 and 164
/// unmodelled, no alternate location, the ligand's every atom placed, His31 doubly protonated and
/// +1 by the rule, and the total +9 with 2604 protein atoms. Reported: how far each entry's Cα
/// atoms are from 181L's, without superposition, since the crystals are isomorphous.
#[test]
fn every_entry_is_181ls_protein_and_builds() {
    let reference = build(&ENTRIES[0]);
    let ca = |s: &System| -> BTreeMap<i32, [f64; 3]> {
        s.residues()
            .iter()
            .filter(|r| r.part == Part::Protein)
            .map(|r| (r.number, atom_at(s, r.number, "CA")))
            .collect()
    };
    let ca0 = ca(&reference);
    let labels = |s: &System| -> Vec<String> {
        s.residues()
            .iter()
            .filter(|r| r.part == Part::Protein)
            .map(|r| r.label())
            .collect()
    };
    let seq0 = seqres(ENTRIES[0].pdb, 'A');
    let conflicts = vec![
        (
            "THR".to_string(),
            54,
            "CYS".to_string(),
            "CONFLICT".to_string(),
        ),
        (
            "ALA".to_string(),
            97,
            "CYS".to_string(),
            "CONFLICT".to_string(),
        ),
        (
            "ALA".to_string(),
            99,
            "LEU".to_string(),
            "CONFLICT".to_string(),
        ),
    ];
    for e in &ENTRIES {
        let header = e.pdb.lines().next().expect("a first line");
        assert!(header.starts_with("HEADER"), "{}", e.code);
        assert_eq!(columns(header, 63, 66), e.code);
        assert_eq!(seqres(e.pdb, 'A'), seq0, "{}: SEQRES", e.code);
        assert_eq!(seqadv(e.pdb), conflicts, "{}: SEQADV", e.code);
        // The file's HETATM residues, counted here.
        let mut het: BTreeMap<String, (std::collections::BTreeSet<String>, usize)> =
            BTreeMap::new();
        for l in atom_lines(e.pdb).filter(|l| l.starts_with("HETATM")) {
            let x = het
                .entry(columns(l, 18, 20).trim().to_string())
                .or_default();
            x.0.insert(columns(l, 22, 26).to_string());
            x.1 += 1;
        }
        let mut names: Vec<&str> = het.keys().map(String::as_str).collect();
        names.sort_unstable();
        let mut want: Vec<&str> = vec!["CL", "HED", "HOH", e.ligand];
        want.sort_unstable();
        assert_eq!(names, want, "{}: HETATM residues", e.code);
        let template = ligand_template(e);
        let heavy = template
            .atoms()
            .iter()
            .filter(|a| a.element != Element::H)
            .count();
        assert_eq!(het[e.ligand].0.len(), 1, "{}: one ligand", e.code);
        assert_eq!(het[e.ligand].1, heavy, "{}: every heavy atom", e.code);

        let s = build(e);
        let dropped: BTreeMap<String, (usize, usize)> = s
            .dropped()
            .iter()
            .map(|d| (d.name.clone(), (d.residues, d.atoms)))
            .collect();
        let counted: BTreeMap<String, (usize, usize)> = het
            .iter()
            .filter(|(k, _)| k.as_str() != e.ligand)
            .map(|(k, (r, a))| (k.clone(), (r.len(), *a)))
            .collect();
        assert_eq!(dropped, counted, "{}: dropped", e.code);
        assert_eq!(s.alternates_dropped(), 0, "{}", e.code);
        assert_eq!(labels(&s), labels(&reference), "{}: residues", e.code);
        let unmodelled: Vec<(i32, &str)> = s
            .unmodelled()
            .iter()
            .map(|u| (u.number, u.name.as_str()))
            .collect();
        assert_eq!(unmodelled, [(163, "ASN"), (164, "LEU")], "{}", e.code);
        assert_eq!(
            s.atoms_in(Part::Ligand).len(),
            template.atoms().len(),
            "{}: the ligand with its hydrogens",
            e.code
        );
        let his = s
            .residues()
            .iter()
            .find(|r| r.name == "HIS")
            .expect("His31");
        assert_eq!(
            (his.number, his.histidine, his.charge),
            (31, Some(Histidine::Both), 1),
            "{}: His31",
            e.code
        );
        assert_eq!(
            s.residues().iter().filter(|r| r.name == "HIS").count(),
            1,
            "{}",
            e.code
        );
        assert_eq!(s.formal_charge(), 9, "{}", e.code);
        assert_eq!(s.atoms_in(Part::Protein).len(), 2604, "{}", e.code);

        let bridge = ["ND1", "NE2"]
            .iter()
            .flat_map(|n| ["OD1", "OD2"].map(|o| (*n, o)))
            .map(|(n, o)| len(sub(atom_at(&s, 31, n), atom_at(&s, 70, o))) / ANGSTROM)
            .fold(f64::INFINITY, f64::min);
        let ca1 = ca(&s);
        let moves: Vec<(i32, f64)> = ca0
            .iter()
            .map(|(k, p)| (*k, len(sub(ca1[k], *p)) / ANGSTROM))
            .collect();
        let rms = (moves.iter().map(|m| m.1 * m.1).sum::<f64>() / moves.len() as f64).sqrt();
        let worst = moves
            .iter()
            .fold((0, 0.0f64), |w, m| if m.1 > w.1 { *m } else { w });
        let over: Vec<i32> = moves.iter().filter(|m| m.1 > 0.5).map(|m| m.0).collect();
        let lig = |s: &System| {
            let at: Vec<[f64; 3]> = s
                .atoms_in(Part::Ligand)
                .iter()
                .map(|&i| s.component().atoms()[i].at)
                .collect();
            let n = at.len() as f64;
            at.iter().fold([0.0; 3], |c, p| {
                [c[0] + p[0] / n, c[1] + p[1] / n, c[2] + p[2] / n]
            })
        };
        // The closest ligand–protein contacts at the crystal pose, hydrogens included: the
        // ligand's hydrogens are placed by its template and the protein's by 2c-1's rules, and
        // neither placement sees the other partner's hydrogens.
        let atoms = s.component().atoms();
        let mut contacts: Vec<(f64, usize, usize)> = s
            .atoms_in(Part::Ligand)
            .iter()
            .flat_map(|&l| {
                s.atoms_in(Part::Protein)
                    .into_iter()
                    .map(move |p| (len(sub(atoms[l].at, atoms[p].at)) / ANGSTROM, l, p))
            })
            .collect();
        contacts.sort_by(|a, b| a.0.total_cmp(&b.0));
        let closest: Vec<String> = contacts[..3]
            .iter()
            .map(|(d, l, p)| format!("{} to {} {d:.2} Å", atoms[*l].name, atoms[*p].name))
            .collect();
        eprintln!(
            "{} {} ({}): {} atoms; dropped {dropped:?}; His31–Asp70 N···O {bridge:.2} Å; Cα from \
             181L's: RMS {rms:.3} Å, worst {:.3} Å at residue {}, over 0.5 Å at {over:?}; ligand \
             centroid {:.2} Å from benzene's; closest contacts {closest:?}",
            e.code,
            e.ligand,
            e.name,
            atoms.len(),
            worst.1,
            worst.0,
            len(sub(lig(&s), lig(&reference))) / ANGSTROM
        );
    }
}

// --- Statistics ---------------------------------------------------------------------------------

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

/// `Σ (x − x̄)(y − ȳ)`.
fn co(x: &[f64], y: &[f64]) -> f64 {
    let (mx, my) = (mean(x), mean(y));
    x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum()
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    co(x, y) / (co(x, x) * co(y, y)).sqrt()
}

/// Ranks from 1, ties given their mean rank.
fn ranks(x: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..x.len()).collect();
    order.sort_by(|&a, &b| x[a].total_cmp(&x[b]));
    let mut r = vec![0.0; x.len()];
    let mut k = 0;
    while k < order.len() {
        let mut m = k;
        while m + 1 < order.len() && x[order[m + 1]] == x[order[k]] {
            m += 1;
        }
        let rank = (k + m) as f64 / 2.0 + 1.0;
        for &i in &order[k..=m] {
            r[i] = rank;
        }
        k = m + 1;
    }
    r
}

fn spearman(x: &[f64], y: &[f64]) -> f64 {
    pearson(&ranks(x), &ranks(y))
}

/// The least-squares slope of `y` on `x`.
fn slope(x: &[f64], y: &[f64]) -> f64 {
    co(x, y) / co(x, x)
}

/// The RMS of `y − x` after its mean is removed.
fn rms_after_offset(x: &[f64], y: &[f64]) -> f64 {
    let d: Vec<f64> = x.iter().zip(y).map(|(a, b)| b - a).collect();
    let m = mean(&d);
    (d.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / d.len() as f64).sqrt()
}

/// The 95% interval of a correlation `r` on `n` points by Fisher's z-transform: `tanh(atanh r ±
/// 1.96 s)` with `s = √(c / (n − 3))`, `c` = 1 for Pearson's r. For Spearman's ρ the factor
/// `c` = 1.06 is the one usually attributed to Fieller, Hartley and Pearson, *Biometrika* **44**,
/// 470 (1957) — quoted from memory, that paper not read.
fn fisher_interval(r: f64, n: usize, c: f64) -> (f64, f64) {
    let z = r.atanh();
    let s = (c / (n as f64 - 3.0)).sqrt();
    ((z - 1.96 * s).tanh(), (z + 1.96 * s).tanh())
}

/// **The statistics against closed forms**, since the series reports them and asserts none:
/// `y = x²` on 1…5 has `Σ(x − x̄)(y − ȳ)` = 60, `Σ(x − x̄)²` = 10 and `Σ(y − ȳ)²` = 374, so
/// Pearson's r is `60/√3740`, the slope 6, Spearman's ρ exactly 1, and `y − x` = 0, 2, 6, 12, 20
/// has RMS `√52.8` about its mean; a reversed order is −1 both ways; tied values share their mean
/// rank; and Fisher's interval at r = 1/2, n = 12 is `tanh(ln 3 / 2 ± 0.6533…)`, with r = 0's
/// symmetric about zero.
#[test]
fn the_statistics_against_closed_forms() {
    let eps = 8.0 * f64::EPSILON;
    let x = [1.0, 2.0, 3.0, 4.0, 5.0];
    let y = [1.0, 4.0, 9.0, 16.0, 25.0];
    assert!((pearson(&x, &y) - 60.0 / 3740f64.sqrt()).abs() <= eps);
    assert!((slope(&x, &y) - 6.0).abs() <= eps * 6.0);
    assert_eq!(spearman(&x, &y), 1.0);
    assert!((rms_after_offset(&x, &y) - 52.8f64.sqrt()).abs() <= eps * 8.0);
    let rev = [5.0, 4.0, 3.0, 2.0, 1.0];
    assert!((pearson(&x, &rev) + 1.0).abs() <= eps);
    assert!((spearman(&x, &rev) + 1.0).abs() <= eps);
    assert_eq!(ranks(&[3.0, 2.0, 2.0, 1.0]), [4.0, 2.5, 2.5, 1.0]);
    let (lo, hi) = fisher_interval(0.5, 12, 1.0);
    let z = 3f64.ln() / 2.0;
    assert!((lo - (z - 1.96 / 3.0).tanh()).abs() <= eps);
    assert!((hi - (z + 1.96 / 3.0).tanh()).abs() <= eps);
    let (lo, hi) = fisher_interval(0.0, 12, 1.0);
    assert!((lo + hi).abs() <= eps && hi > 0.0);
}

// --- The series ---------------------------------------------------------------------------------

/// The binding terms at one pose, kcal/mol and Å².
#[derive(Clone, Copy, Debug)]
struct Terms {
    vdw: f64,
    elec: f64,
    polar: f64,
    area: f64,
    nonpolar: f64,
}

impl Terms {
    fn of(b: &Binding) -> Terms {
        let i = b.interaction();
        let d = b.desolvation(POINTS);
        Terms {
            vdw: kcal(i.van_der_waals),
            elec: kcal(i.electrostatic),
            polar: kcal(d.polar),
            area: d.buried_area / (ANGSTROM * ANGSTROM),
            nonpolar: kcal(d.nonpolar),
        }
    }

    fn vacuum(&self) -> f64 {
        self.vdw + self.elec
    }

    fn total(&self) -> f64 {
        self.vacuum() + self.polar + self.nonpolar
    }
}

/// A contact: its distance in Å and the two atoms' names in the system.
type Contact = (f64, String, String);

/// The closest pair of a ligand hydrogen and a pocket hydrogen at `at` (complex order): the
/// distance in Å and the two atoms' names in the system.
fn closest_hydrogens(s: &System, b: &Binding, at: &[[f64; 3]]) -> Contact {
    let (el, ids) = (b.elements(), b.system_atoms());
    let mut best = (f64::INFINITY, 0, 0);
    for l in b.ligand_range().filter(|&l| el[l] == Element::H) {
        for p in (0..b.pocket_len()).filter(|&p| el[p] == Element::H) {
            let d = len(sub(at[l], at[p])) / ANGSTROM;
            if d < best.0 {
                best = (d, l, p);
            }
        }
    }
    let name = |k: usize| s.component().atoms()[ids[k]].name.clone();
    (best.0, name(best.1), name(best.2))
}

/// The crystal's heavy atoms with **every hydrogen relaxed**, protein's and ligand's, on the
/// complex's whole force field: `(van der Waals, Coulomb)` of ΔE_bind there, kcal/mol, from the
/// three force fields' term-by-term difference at the relaxed positions, and the closest H–H
/// contact before and after. Not a [`Binding`] (whose pocket cannot move), so not solvated.
///
/// **Why it is here.** At the crystal pose both partners' hydrogens are placed — the ligand's by
/// its template, the protein's by 2c-1's rules, whose rotor step clears heavy atoms only — and
/// neither placement sees the other's hydrogens. On three of the nine entries that leaves an H–H
/// contact UFF scores in the hundreds of kcal/mol (183L's indene H12 1.41 Å from Val111's HG13),
/// so the crystal pose's vacuum ΔE measures the placement, not the binding. Relaxing hydrogens
/// alone keeps every measured coordinate.
fn hydrogens_relaxed(s: &System, b: &Binding) -> ((f64, f64), Contact, Contact) {
    let frozen: Vec<bool> = b.elements().iter().map(|&e| e != Element::H).collect();
    let mut at = b.positions().to_vec();
    let before = closest_hydrogens(s, b, &at);
    let ff = b.force_field();
    let mut m = Minimiser::new(TOLERANCE).with_frozen(frozen);
    let mut p = m.step(ff, &[], &mut at);
    while p.status == Status::Running && p.steps < 50000 {
        p = m.step(ff, &[], &mut at);
    }
    assert_eq!(p.status, Status::Converged, "hydrogens relaxed: {p:?}");
    let n0 = b.pocket_len();
    let (c, q, l) = (
        ff.energy(&at),
        b.pocket_force_field().energy(&at[..n0]),
        b.ligand_force_field().energy(&at[n0..]),
    );
    let after = closest_hydrogens(s, b, &at);
    (
        (
            kcal(c.van_der_waals - q.van_der_waals - l.van_der_waals),
            kcal(c.electrostatic - q.electrostatic - l.electrostatic),
        ),
        before,
        after,
    )
}

/// One entry at one cutoff: the crystal pose, the crystal's heavy atoms with every hydrogen
/// relaxed, and the ligand minimised in the frozen pocket.
struct Measured {
    residues: usize,
    atoms: usize,
    pocket_charge: i32,
    crystal: Terms,
    relaxed: (f64, f64),
    contact: Contact,
    relaxed_contact: f64,
    minimised: Terms,
    rmsd: f64,
    steps: usize,
}

/// Builds the pocket, measures the crystal pose, relaxes every hydrogen and measures the vacuum
/// terms, minimises the ligand in the frozen pocket and measures again — asserting the facts that
/// must hold whatever the energies are: ΔSASA negative at both poses, every protein atom unmoved
/// to the bit, the minimiser converged, and van der Waals attractive at the minimised pose and
/// with the hydrogens relaxed. **Not at the crystal pose**, where it is not a fact: see
/// [`hydrogens_relaxed`].
fn measure(e: &Entry, s: &System, cutoff: f64) -> Measured {
    let b = Binding::new(s, cutoff * ANGSTROM).unwrap_or_else(|err| panic!("{}: {err}", e.code));
    let crystal = Terms::of(&b);
    let (relaxed, contact, after) = hydrogens_relaxed(s, &b);
    let mut m = b.clone();
    let p = m.minimise_ligand(20000, TOLERANCE);
    assert_eq!(p.status, Status::Converged, "{} {cutoff} Å: {p:?}", e.code);
    let n0 = b.pocket_len();
    for i in 0..n0 {
        assert_eq!(
            m.positions()[i].map(f64::to_bits),
            b.positions()[i].map(f64::to_bits),
            "{} {cutoff} Å: protein atom {i} moved",
            e.code
        );
    }
    let minimised = Terms::of(&m);
    for (pose, t) in [("crystal", &crystal), ("minimised", &minimised)] {
        assert!(
            t.area < 0.0,
            "{} {cutoff} Å {pose}: ΔSASA {}",
            e.code,
            t.area
        );
    }
    for (pose, vdw) in [
        ("minimised", minimised.vdw),
        ("hydrogens relaxed", relaxed.0),
    ] {
        assert!(vdw < 0.0, "{} {cutoff} Å {pose}: vdW {vdw}", e.code);
    }
    Measured {
        residues: b.residues().len(),
        atoms: b.positions().len(),
        pocket_charge: b.formal_charges().0,
        crystal,
        relaxed,
        contact,
        relaxed_contact: after.0,
        minimised,
        rmsd: m.ligand_rmsd() / ANGSTROM,
        steps: p.steps,
    }
}

/// 181L as 2c-2 committed it (`CHANGELOG.md`, the cutoff table): at 6 and 8 Å, van der Waals,
/// Coulomb, polar, ΔSASA, nonpolar, the minimised pose's RMSD and its vacuum ΔE — each held to
/// half a unit of the last digit printed there.
fn assert_181l_unchanged(cutoff: f64, m: &Measured) {
    let committed: [(f64, f64); 8] = match cutoff as i32 {
        6 => [
            (-11.295, 5e-4),
            (-1.133, 5e-4),
            (7.322, 5e-4),
            (-294.6, 5e-2),
            (-1.473, 5e-4),
            (0.809, 5e-4),
            (-22.346, 5e-4),
            (18.0, 0.0),
        ],
        8 => [
            (-13.048, 5e-4),
            (-1.104, 5e-4),
            (10.607, 5e-4),
            (-294.6, 5e-2),
            (-1.473, 5e-4),
            (0.814, 5e-4),
            (-24.149, 5e-4),
            (43.0, 0.0),
        ],
        _ => unreachable!(),
    };
    let got = [
        m.crystal.vdw,
        m.crystal.elec,
        m.crystal.polar,
        m.crystal.area,
        m.crystal.nonpolar,
        m.rmsd,
        m.minimised.vacuum(),
        m.residues as f64,
    ];
    let names = [
        "vdW", "elec", "polar", "ΔSASA", "nonpolar", "RMSD", "ΔE min", "residues",
    ];
    for ((g, (w, tol)), name) in got.iter().zip(committed).zip(names) {
        assert!(
            (g - w).abs() <= tol,
            "181L {cutoff} Å {name}: {g} against the committed {w}"
        );
    }
}

/// The statistics of one computed column against experiment, as a table row.
fn statistics(label: &str, exp: &[f64], calc: &[f64]) -> String {
    let n = exp.len();
    let r = pearson(exp, calc);
    let rho = spearman(exp, calc);
    let (rl, rh) = fisher_interval(r, n, 1.0);
    let (sl, sh) = fisher_interval(rho, n, 1.06);
    format!(
        "| {label} | {r:+.2} [{rl:+.2}, {rh:+.2}] | {rho:+.2} [{sl:+.2}, {sh:+.2}] | {:+.2} | \
         {:.2} | {:.2} |",
        slope(exp, calc),
        rms_after_offset(exp, calc),
        calc.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
            - calc.iter().cloned().fold(f64::INFINITY, f64::min)
    )
}

/// **The series, against experiment, reported and not asserted.** For each entry and each of 6 and
/// 8 Å: the rigid-pocket ΔE in vacuum, the polar and nonpolar desolvation and their total at the
/// crystal pose with its closest H–H contact; the vacuum ΔE with every hydrogen relaxed on the
/// crystal's heavy atoms ([`hydrogens_relaxed`]); and the crystal-pose quantities again after the
/// ligand is minimised in the frozen pocket, with the pose's RMSD. Asserted in [`measure`], and
/// listed there; and 181L's numbers as 2c-2 committed them. Then Pearson's r and Spearman's ρ of each column against ΔG°_exp, each with its 95%
/// Fisher interval, the slope of the column on ΔG°_exp, and its RMS about ΔG°_exp after the mean
/// offset is removed — printed, because nine points cannot carry a claim. Ignored: two cutoffs of
/// nine entries is eighteen pairs of QEq solves and twice as many minimisations, 100 s with
/// `--release`, and unoptimised it would be many times that.
#[test]
#[ignore = "nine entries at two cutoffs: QEq and a minimisation each; minutes with --release -- --ignored"]
fn the_congener_series_against_experiment() {
    let exp: Vec<f64> = ENTRIES.iter().map(|e| e.dg.0).collect();
    for cutoff in [6.0, 8.0] {
        eprintln!(
            "\n{cutoff} Å, kcal/mol:\n| entry | ligand | ΔG°exp | residues | atoms | pocket q | \
             vdW | elec | ΔE vac | polar | ΔSASA Å² | nonpolar | total | closest H–H Å | \
             H relaxed: vdW | elec | ΔE vac | closest H–H Å | → RMSD Å | steps | vdW | elec | \
             ΔE vac | polar | ΔSASA Å² | nonpolar | total |"
        );
        let mut rows = Vec::new();
        for e in &ENTRIES {
            let s = build(e);
            let m = measure(e, &s, cutoff);
            if e.code == "181L" {
                assert_181l_unchanged(cutoff, &m);
            }
            let (c, n) = (&m.crystal, &m.minimised);
            eprintln!(
                "| {} | {} | {:.2} ± {:.2} | {} | {} | {:+} | {:.2} | {:.2} | {:.2} | {:+.2} | {:.1} \
                 | {:.2} | {:.2} | {:.2} ({} to {}) | {:.2} | {:.2} | {:.2} | {:.2} | {:.3} | {} | \
                 {:.2} | {:.2} | {:.2} | {:+.2} | {:.1} | {:.2} | {:.2} |",
                e.code,
                e.name,
                e.dg.0,
                e.dg.1,
                m.residues,
                m.atoms,
                m.pocket_charge,
                c.vdw,
                c.elec,
                c.vacuum(),
                c.polar,
                c.area,
                c.nonpolar,
                c.total(),
                m.contact.0,
                m.contact.1,
                m.contact.2,
                m.relaxed.0,
                m.relaxed.1,
                m.relaxed.0 + m.relaxed.1,
                m.relaxed_contact,
                m.rmsd,
                m.steps,
                n.vdw,
                n.elec,
                n.vacuum(),
                n.polar,
                n.area,
                n.nonpolar,
                n.total()
            );
            rows.push(m);
        }
        let col = |f: &dyn Fn(&Measured) -> f64| rows.iter().map(f).collect::<Vec<f64>>();
        eprintln!(
            "\n{cutoff} Å against ΔG°exp, n = {}; 95% Fisher intervals in brackets:\n| column | \
             Pearson r | Spearman ρ | slope on ΔG°exp | RMS after offset | spread |",
            exp.len()
        );
        for (label, calc) in [
            ("ΔE vacuum, crystal", col(&|m| m.crystal.vacuum())),
            ("total, crystal", col(&|m| m.crystal.total())),
            (
                "ΔE vacuum, hydrogens relaxed",
                col(&|m| m.relaxed.0 + m.relaxed.1),
            ),
            ("ΔE vacuum, minimised", col(&|m| m.minimised.vacuum())),
            ("total, minimised", col(&|m| m.minimised.total())),
        ] {
            eprintln!("{}", statistics(label, &exp, &calc));
        }
    }
    let spread = exp.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        - exp.iter().cloned().fold(f64::INFINITY, f64::min);
    eprintln!(
        "\nΔG°exp spans {spread:.2} kcal/mol; for n = 9 a Pearson r is different from zero at the \
         two-sided 5% level only past {:.3}",
        {
            // t = r √(n − 2) / √(1 − r²) at t(0.975, 7) = 2.3646.
            let t: f64 = 2.3646;
            t / (t * t + 7.0).sqrt()
        }
    );
}
