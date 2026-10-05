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
//! **Each system's hydrogens are relaxed** ([`Binding::relaxing_hydrogens`], step A-1) as well as
//! taken as placed: at the crystal pose three of the nine clash hydrogen to hydrogen, and the
//! relaxed columns are the ones that measure binding rather than the placement.
//!
//! **What is asserted about the binding energies is only what must hold**: van der Waals
//! attraction, buried area, frozen atoms, convergence, no hydrogen–hydrogen contact under 1.8 Å
//! once relaxed, and 181L's numbers as 2c-2 committed them.
//! **No correlation with experiment is asserted.** The series is nine ligands of one size whose
//! measured ΔG spans 2.1 kcal/mol; a correlation measured on nine points has a 95% interval that
//! the series test computes and prints, and it is wide.
//!
//! **The polar term is split into its three parts** (step A-2) — the pocket desolvated by the
//! ligand's volume, the ligand's own, the screened cross terms — and given beside the empty-cavity
//! reference ([`ApoCavity::Empty`]), at 6 and 8 Å and, in a test of its own, at 20 Å, where it
//! has converged.
//!
//! The series takes minutes, so it is ignored by default:
//! `cargo test -p pantometry-forcefield --release --test the_congener_series -- --ignored`.

mod protein;

use pantometry_forcefield::minimise::KCAL_PER_MOL_ANGSTROM;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{
    ApoCavity, Binding, BindingError, Component, Element, Histidine, Part, Selection, Status,
    System,
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

/// The most steps one hydrogen relaxation may take: reaching it is a failure to converge.
const RELAX_STEPS: usize = 20000;

/// The tolerance hydrogens are relaxed to, and each ligand minimised to from a relaxed complex:
/// [`Binding::HYDROGEN_TOLERANCE`], 2e-3 kcal mol⁻¹ Å⁻¹. At 2c-2's 1e-4 the minimiser stalled four
/// times on these pockets, where their energies' rounding can hide a step's decrease; that
/// constant gives the evidence, which is empirical, and the 6 Å table this test prints compares
/// the two tolerances on all nine. The ligand minimised from the crystal keeps 2c-2's
/// [`TOLERANCE`].
const H_TOL: f64 = Binding::HYDROGEN_TOLERANCE;

/// The closest a ligand hydrogen may come to a pocket hydrogen once relaxed, Å: where UFF's H···H
/// pair is already +11 kcal/mol, as `benzene_in_its_pocket.rs` computes.
const H_H_CONTACT: f64 = 1.8;

/// The binding terms at one pose, kcal/mol and Å².
#[derive(Clone, Copy, Debug)]
struct Terms {
    vdw: f64,
    elec: f64,
    reorganisation: f64,
    /// The ligand's part of `reorganisation`; the rest is the pocket's.
    reorganisation_ligand: f64,
    polar: f64,
    /// The polar term's three parts, which sum to it: the pocket desolvated by the ligand's
    /// volume (the ligand's charges zero), the ligand's own (the pocket's zero), and the screened
    /// cross terms (the rest).
    pocket_part: f64,
    ligand_part: f64,
    cross_part: f64,
    /// The polar term with the apo cavity empty ([`ApoCavity::Empty`]).
    polar_empty: f64,
    area: f64,
    nonpolar: f64,
}

impl Terms {
    fn of(b: &Binding) -> Terms {
        let i = b.interaction();
        let d = b.desolvation(POINTS);
        let (n0, n) = (b.pocket_len(), b.positions().len());
        let zeroed = |r: std::ops::Range<usize>| {
            let mut q = b.charges().to_vec();
            q[r].fill(0.0);
            kcal(b.desolvation_with(&q, 1).polar)
        };
        let (pocket_part, ligand_part) = (zeroed(n0..n), zeroed(0..n0));
        Terms {
            vdw: kcal(i.van_der_waals),
            elec: kcal(i.electrostatic),
            reorganisation: kcal(i.reorganisation.total),
            reorganisation_ligand: kcal(
                b.ligand_force_field().energy(b.ligand_positions()).total
                    - b.ligand_force_field()
                        .energy(b.ligand_alone_positions())
                        .total,
            ),
            polar: kcal(d.polar),
            pocket_part,
            ligand_part,
            cross_part: kcal(d.polar) - pocket_part - ligand_part,
            polar_empty: kcal(
                b.desolvation_with_cavity(b.charges(), 1, ApoCavity::Empty)
                    .polar,
            ),
            area: d.buried_area / (ANGSTROM * ANGSTROM),
            nonpolar: kcal(d.nonpolar),
        }
    }

    fn vacuum(&self) -> f64 {
        self.vdw + self.elec + self.reorganisation
    }

    fn total(&self) -> f64 {
        self.vacuum() + self.polar + self.nonpolar
    }

    /// [`Terms::total`] with the apo cavity empty.
    fn total_empty(&self) -> f64 {
        self.vacuum() + self.polar_empty + self.nonpolar
    }
}

/// A contact: its distance in Å and the two atoms' names in the system.
type Contact = (f64, String, String);

/// The closest pair of a ligand hydrogen and a pocket hydrogen in `b`'s complex: the distance in Å
/// and the two atoms' names in the system.
fn closest_hydrogens(s: &System, b: &Binding) -> Contact {
    let (el, ids, at) = (b.elements(), b.system_atoms(), b.positions());
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

/// One entry at one cutoff: the crystal pose; each system's hydrogens relaxed
/// ([`Binding::relaxing_hydrogens`]), and the same with every pocket hydrogen free instead of the
/// rule's; the ligand minimised in the frozen pocket from the crystal, as 2c-2 and 2c-3a did; and
/// minimised from the relaxed complex, its free hydrogens with it.
struct Measured {
    residues: usize,
    atoms: usize,
    pocket_charge: i32,
    crystal: Terms,
    contact: Contact,
    relaxed: Terms,
    relaxed_contact: Contact,
    every_hydrogen_free: f64,
    minimised: Terms,
    rmsd: f64,
    relaxed_minimised: Terms,
    relaxed_rmsd: f64,
    relaxed_minimised_contact: f64,
    steps: usize,
    /// At 6 Å: the relaxed and the relaxed-minimised ΔE_bind and pose again at 2c-2's 1e-4, as
    /// `(|ΔE relaxed|, |ΔE minimised|, ligand RMS displacement Å, every run converged)`.
    tight: Option<(f64, f64, f64, bool)>,
}

/// Asserts that `b`'s three relaxations converged, and that every atom they hold is where `start`
/// has it, bit for bit, in all three systems — every pocket atom the mask holds always, and the
/// ligand's heavy atoms too unless `ligand_moved`, in which case the ligand alone's heavy atoms are
/// the complex's. **A stall is a failure**, every-hydrogen-free comparison included: a stall means
/// the force is still above the tolerance, so accepting one would accept it at any force.
fn assert_relaxed(
    (code, cutoff, pose): (&str, f64, &str),
    start: &Binding,
    b: &Binding,
    ligand_moved: bool,
) {
    let r = b.hydrogen_relaxation().expect("relaxed");
    for (sys, p) in [
        ("complex", r.complex),
        ("pocket", r.pocket),
        ("ligand", r.ligand),
    ] {
        assert_eq!(
            p.status,
            Status::Converged,
            "{code} {cutoff} Å {pose} {sys}: {p:?}"
        );
    }
    let n0 = b.pocket_len();
    let held = |i: usize| !r.free[i];
    for (sys, at, first) in [
        ("complex", b.positions(), 0),
        ("pocket", b.pocket_alone_positions(), 0),
        ("ligand", b.ligand_alone_positions(), n0),
    ] {
        for (k, p) in at.iter().enumerate() {
            let i = first + k;
            if !held(i) {
                continue;
            }
            let want = if ligand_moved && i >= n0 {
                if sys == "complex" {
                    continue;
                }
                b.positions()[i]
            } else {
                start.positions()[i]
            };
            assert_eq!(
                p.map(f64::to_bits),
                want.map(f64::to_bits),
                "{code} {cutoff} Å {sys}: held atom {i} moved"
            );
        }
    }
}

/// Builds the pocket and measures every pose, asserting the facts that must hold whatever the
/// energies are: ΔSASA negative at the two poses whose systems share positions (the crystal and
/// the unrelaxed minimum), every held atom unmoved to the bit, every relaxation and minimisation
/// converged, van der Waals attractive at every relaxed or minimised pose, and **no ligand
/// hydrogen within 1.8 Å of a pocket hydrogen once relaxed**. Not van der Waals at the crystal
/// pose, where the placed hydrogens clash: that is the defect the relaxation removes.
fn measure(e: &Entry, s: &System, cutoff: f64) -> Measured {
    let b = Binding::new(s, cutoff * ANGSTROM).unwrap_or_else(|err| panic!("{}: {err}", e.code));
    let crystal = Terms::of(&b);
    let contact = closest_hydrogens(s, &b);

    let r = b.relaxing_hydrogens(RELAX_STEPS, H_TOL);
    assert_relaxed((e.code, cutoff, "relaxed"), &b, &r, false);
    let relaxed = Terms::of(&r);
    let relaxed_contact = closest_hydrogens(s, &r);
    let all: Vec<bool> = b.elements().iter().map(|&x| x == Element::H).collect();
    let every = b.relaxing_hydrogens_of(&all, RELAX_STEPS, H_TOL);
    assert_relaxed((e.code, cutoff, "every H free"), &b, &every, false);
    let every_hydrogen_free = kcal(every.interaction().total());

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

    let mut rm = r.clone();
    let q = rm.minimise_ligand(20000, H_TOL);
    assert_eq!(q.status, Status::Converged, "{} {cutoff} Å: {q:?}", e.code);
    assert_relaxed((e.code, cutoff, "relaxed, minimised"), &b, &rm, true);
    let relaxed_minimised = Terms::of(&rm);
    let relaxed_minimised_contact = closest_hydrogens(s, &rm).0;

    for (pose, t) in [("crystal", &crystal), ("minimised", &minimised)] {
        assert!(
            t.area < 0.0,
            "{} {cutoff} Å {pose}: ΔSASA {}",
            e.code,
            t.area
        );
    }
    for (pose, d) in [
        ("hydrogens relaxed", relaxed_contact.0),
        ("relaxed and minimised", relaxed_minimised_contact),
    ] {
        assert!(
            d >= H_H_CONTACT,
            "{} {cutoff} Å {pose}: closest H–H {d} Å",
            e.code
        );
    }
    for (pose, vdw) in [
        ("minimised", minimised.vdw),
        ("hydrogens relaxed", relaxed.vdw),
        ("relaxed and minimised", relaxed_minimised.vdw),
    ] {
        assert!(vdw < 0.0, "{} {cutoff} Å {pose}: vdW {vdw}", e.code);
    }
    Measured {
        residues: b.residues().len(),
        atoms: b.positions().len(),
        pocket_charge: b.formal_charges().0,
        crystal,
        contact,
        relaxed,
        relaxed_contact,
        every_hydrogen_free,
        minimised,
        rmsd: m.ligand_rmsd() / ANGSTROM,
        relaxed_minimised,
        relaxed_rmsd: rm.ligand_rmsd() / ANGSTROM,
        relaxed_minimised_contact,
        steps: q.steps,
        tight: (cutoff == 6.0).then(|| {
            let t = b.relaxing_hydrogens(RELAX_STEPS, TOLERANCE);
            let mut tm = t.clone();
            let p = tm.minimise_ligand(20000, TOLERANCE);
            let h = [t.hydrogen_relaxation(), tm.hydrogen_relaxation()];
            let converged = p.status == Status::Converged
                && h.iter().flatten().all(|h| {
                    [h.complex, h.pocket, h.ligand]
                        .iter()
                        .all(|x| x.status == Status::Converged)
                });
            let lig = rm.ligand_range();
            let rms = (lig
                .clone()
                .map(|k| len(sub(rm.positions()[k], tm.positions()[k])).powi(2))
                .sum::<f64>()
                / lig.len() as f64)
                .sqrt();
            (
                (kcal(t.interaction().total()) - relaxed.vacuum()).abs(),
                (kcal(tm.interaction().total()) - relaxed_minimised.vacuum()).abs(),
                rms / ANGSTROM,
                converged,
            )
        }),
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

/// The header of [`polar_row`]'s table, the second pose named `pose`.
fn polar_header(pose: &str) -> String {
    format!(
        "| entry | ligand | crystal: polar | pocket by volume | ligand's own | cross | empty cavity \
         | {pose}: polar | pocket by volume | ligand's own | cross | empty cavity | total | total, \
         cavity empty |"
    )
}

/// One entry's polar term split into its three parts, at the crystal pose and relaxed and
/// minimised (or, in the converged pocket, minimised from the crystal), with the empty-cavity
/// reference beside each.
fn polar_row(e: &Entry, c: &Terms, q: &Terms) -> String {
    format!(
        "| {} | {} | {:+.2} | {:+.2} | {:+.2} | {:+.2} | {:+.2} | {:+.2} | {:+.2} | {:+.2} | {:+.2} \
         | {:+.2} | {:+.2} | {:+.2} |",
        e.code,
        e.name,
        c.polar,
        c.pocket_part,
        c.ligand_part,
        c.cross_part,
        c.polar_empty,
        q.polar,
        q.pocket_part,
        q.ligand_part,
        q.cross_part,
        q.polar_empty,
        q.total(),
        q.total_empty()
    )
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
/// 8 Å: the crystal pose with its closest H–H contact; each system's hydrogens relaxed, with the
/// reorganisation, the solvation at each system's own positions and the closest H–H after; the
/// vacuum ΔE with every pocket hydrogen free instead; the ligand minimised from the crystal (2c-3a's
/// column); and minimised from the relaxed complex, with its RMSD. Asserted in [`measure`] and
/// listed there; and 181L's numbers as 2c-2 committed them. Then Pearson's r and Spearman's ρ of
/// each column against ΔG°_exp, each with its 95% Fisher interval, the slope of the column on
/// ΔG°_exp, and its RMS about ΔG°_exp after the mean offset is removed — printed, because nine
/// points cannot carry a claim. Ignored: two cutoffs of nine entries is eighteen pairs of QEq solves
/// and many relaxations and minimisations, minutes with `--release`, and unoptimised it would be
/// many times that.
#[test]
#[ignore = "nine entries at two cutoffs: QEq, relaxations and minimisations each; minutes with --release -- --ignored"]
fn the_congener_series_against_experiment() {
    let exp: Vec<f64> = ENTRIES.iter().map(|e| e.dg.0).collect();
    for cutoff in [6.0, 8.0] {
        eprintln!(
            "\n{cutoff} Å, kcal/mol:\n| entry | ligand | ΔG°exp | residues | atoms | pocket q | \
             crystal: vdW | elec | ΔE vac | polar | total | closest H–H Å | relaxed: vdW | elec | \
             reorg (ligand's) | ΔE vac | polar | ΔSASA Å² | nonpolar | total | closest H–H Å | \
             every H free: \
             ΔE vac | minimised from crystal: RMSD Å | ΔE vac | total | relaxed, minimised: RMSD Å \
             | steps | vdW | elec | reorg (ligand's) | ΔE vac | polar | ΔSASA Å² | nonpolar | total | closest \
             H–H Å |"
        );
        let mut rows = Vec::new();
        for e in &ENTRIES {
            let s = build(e);
            let m = measure(e, &s, cutoff);
            if e.code == "181L" {
                assert_181l_unchanged(cutoff, &m);
            }
            let (c, r, n, q) = (&m.crystal, &m.relaxed, &m.minimised, &m.relaxed_minimised);
            eprintln!(
                "| {} | {} | {:.2} ± {:.2} | {} | {} | {:+} | {:.2} | {:.2} | {:.2} | {:+.2} | {:.2} \
                 | {:.2} ({} to {}) | {:.2} | {:.2} | {:+.2} ({:+.2}) | {:.2} | {:+.2} | {:.1} | {:.2} | \
                 {:.2} | {:.2} ({} to {}) | {:.2} | {:.3} | {:.2} | {:.2} | {:.3} | {} | {:.2} | \
                 {:.2} | {:+.2} ({:+.2}) | {:.2} | {:+.2} | {:.1} | {:.2} | {:.2} | {:.2} |",
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
                c.total(),
                m.contact.0,
                m.contact.1,
                m.contact.2,
                r.vdw,
                r.elec,
                r.reorganisation,
                r.reorganisation_ligand,
                r.vacuum(),
                r.polar,
                r.area,
                r.nonpolar,
                r.total(),
                m.relaxed_contact.0,
                m.relaxed_contact.1,
                m.relaxed_contact.2,
                m.every_hydrogen_free,
                m.rmsd,
                n.vacuum(),
                n.total(),
                m.relaxed_rmsd,
                m.steps,
                q.vdw,
                q.elec,
                q.reorganisation,
                q.reorganisation_ligand,
                q.vacuum(),
                q.polar,
                q.area,
                q.nonpolar,
                q.total(),
                m.relaxed_minimised_contact
            );
            rows.push(m);
        }
        eprintln!("\n{cutoff} Å, the polar term's parts and the empty cavity, kcal/mol:");
        eprintln!("{}", polar_header("relaxed, minimised"));
        for (e, m) in ENTRIES.iter().zip(&rows) {
            eprintln!("{}", polar_row(e, &m.crystal, &m.relaxed_minimised));
        }
        if cutoff == 6.0 {
            eprintln!(
                "\n6 Å at {:.0e} against 2c-2's 1e-4 kcal/mol/Å, where 1e-4 converges:\n| entry | \
                 |ΔE relaxed| | |ΔE relaxed, minimised| | ligand RMS apart Å | 1e-4 converged |",
                H_TOL / KCAL_PER_MOL_ANGSTROM
            );
            for (e, m) in ENTRIES.iter().zip(&rows) {
                let (a, b, c, ok) = m.tight.expect("measured at 6 Å");
                eprintln!("| {} | {a:.1e} | {b:.1e} | {c:.1e} | {ok} |", e.code);
            }
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
            ("ΔE vacuum, hydrogens relaxed", col(&|m| m.relaxed.vacuum())),
            ("total, hydrogens relaxed", col(&|m| m.relaxed.total())),
            (
                "ΔE vacuum, every hydrogen free",
                col(&|m| m.every_hydrogen_free),
            ),
            (
                "ΔE vacuum, minimised from crystal",
                col(&|m| m.minimised.vacuum()),
            ),
            (
                "total, minimised from crystal",
                col(&|m| m.minimised.total()),
            ),
            (
                "ΔE vacuum, relaxed and minimised",
                col(&|m| m.relaxed_minimised.vacuum()),
            ),
            (
                "total, relaxed and minimised",
                col(&|m| m.relaxed_minimised.total()),
            ),
            (
                "total, minimised from crystal, cavity empty",
                col(&|m| m.minimised.total_empty()),
            ),
            (
                "total, relaxed and minimised, cavity empty",
                col(&|m| m.relaxed_minimised.total_empty()),
            ),
            (
                "polar, relaxed and minimised",
                col(&|m| m.relaxed_minimised.polar),
            ),
            (
                "polar, relaxed and minimised, cavity empty",
                col(&|m| m.relaxed_minimised.polar_empty),
            ),
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

/// The cutoff at which the polar term has converged, Å: on benzene within 0.10 kcal/mol of the
/// whole protein (+14.86 against +14.96), measured by `the_polar_term_diagnosed.rs`, and 0.80 at
/// 15 Å. **Not the whole protein**, because QEq on 184L's whole protein does not settle: its
/// hydrogen iteration is still changing a charge by 1.8e-9 e after the default 100 solves, against
/// a tolerance of 1e-10.
const CONVERGED: f64 = 20.0;

/// The entries whose QEq fails at [`CONVERGED`] Å, which the test skips by name: none, measured
/// on all nine.
const EXPECTED_QEQ_FAILURES: [&str; 0] = [];

/// **The polar term where it has converged**, a [`CONVERGED`] Å pocket, for all nine: the
/// reference the 6 and 8 Å pockets fall short of. At the crystal pose, which is the same geometry
/// at every cutoff, and with the ligand minimised from it in the frozen protein (2c-3a's pose; the
/// hydrogens are not relaxed here, since the rule would free some thousand protein hydrogens in
/// three systems of about two thousand atoms). Each polar term split into its three parts, the
/// empty-cavity reference beside it, and the solvated totals' statistics against ΔG°exp —
/// reported, not asserted. Asserted: the minimisation converged. **An entry whose QEq fails is
/// named and left out of the statistics**, rather than stopping the measurement of the rest, and
/// the skipped entries must be [`EXPECTED_QEQ_FAILURES`]; any other error from building the
/// binding fails the test, and so do statistics over fewer than three entries. Ignored: QEq on
/// about two thousand atoms for each entry.
#[test]
#[ignore = "QEq on pockets of about two thousand atoms for each of nine entries; tens of minutes with --release -- --ignored"]
fn the_polar_term_where_it_has_converged() {
    eprintln!(
        "\n{CONVERGED} Å, kcal/mol:\n{}",
        polar_header("minimised from crystal")
    );
    let mut rows = Vec::new();
    let mut exp = Vec::new();
    let mut skipped: Vec<&str> = Vec::new();
    for e in &ENTRIES {
        let s = build(e);
        let b = match Binding::new(&s, CONVERGED * ANGSTROM) {
            Ok(b) => b,
            // Only QEq failing is a skip; any other refusal is a defect and stops the test.
            Err(err @ BindingError::Charges { .. }) => {
                eprintln!("| {} | {} | not measured: {err} |", e.code, e.name);
                skipped.push(e.code);
                continue;
            }
            Err(err) => panic!("{} {CONVERGED} Å: {err}", e.code),
        };
        let crystal = Terms::of(&b);
        let mut m = b.clone();
        let p = m.minimise_ligand(20000, TOLERANCE);
        assert_eq!(
            p.status,
            Status::Converged,
            "{} {CONVERGED} Å: {p:?}",
            e.code
        );
        let minimised = Terms::of(&m);
        eprintln!(
            "{} ({} residues, {} atoms)",
            polar_row(e, &crystal, &minimised),
            b.residues().len(),
            b.positions().len()
        );
        rows.push((crystal, minimised));
        exp.push(e.dg.0);
    }
    // The entries QEq is known to fail on at this cutoff: none, measured. A new failure, or one
    // that went away, changes what the statistics are over and must be said.
    assert_eq!(skipped, EXPECTED_QEQ_FAILURES, "entries skipped for QEq");
    // Statistics on fewer than three points are not a correlation; on none they are NaN.
    assert!(rows.len() >= 3, "{} entries measured", rows.len());
    let col = |f: &dyn Fn(&(Terms, Terms)) -> f64| rows.iter().map(f).collect::<Vec<f64>>();
    eprintln!(
        "\n{CONVERGED} Å against ΔG°exp, n = {} of {}; 95% Fisher intervals in brackets:\n| \
         column | Pearson r | Spearman ρ | slope on ΔG°exp | RMS after offset | spread |",
        exp.len(),
        ENTRIES.len()
    );
    for (label, calc) in [
        ("ΔE vacuum, minimised from crystal", col(&|m| m.1.vacuum())),
        ("total, minimised from crystal", col(&|m| m.1.total())),
        (
            "total, minimised from crystal, cavity empty",
            col(&|m| m.1.total_empty()),
        ),
        ("polar, crystal", col(&|m| m.0.polar)),
        ("polar, crystal, cavity empty", col(&|m| m.0.polar_empty)),
    ] {
        eprintln!("{}", statistics(label, &exp, &calc));
    }
}
