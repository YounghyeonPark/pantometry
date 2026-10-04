//! What the protein tests share: PDB 181L as fetched, the dictionary's twenty amino acids and
//! benzene, and the selection that keeps T4 lysozyme and its benzene. Everything is embedded with
//! `include_str!`, so the tests run under `wasm32`, which has no filesystem.

#![allow(dead_code)]

use pantometry_forcefield::{Component, Selection, System};

/// PDB entry 181L, T4 lysozyme L99A with benzene, byte for byte as RCSB serves it.
pub const PDB_181L: &str = include_str!("../../components/181L.pdb");

/// Every template code, the twenty amino acids then the ligand.
pub const CODES: [&str; 21] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "BNZ",
];

/// The dictionary entry for `code`, as fetched.
pub fn ccd(code: &str) -> &'static str {
    match code {
        "ALA" => include_str!("../../components/ALA.cif"),
        "ARG" => include_str!("../../components/ARG.cif"),
        "ASN" => include_str!("../../components/ASN.cif"),
        "ASP" => include_str!("../../components/ASP.cif"),
        "CYS" => include_str!("../../components/CYS.cif"),
        "GLN" => include_str!("../../components/GLN.cif"),
        "GLU" => include_str!("../../components/GLU.cif"),
        "GLY" => include_str!("../../components/GLY.cif"),
        "HIS" => include_str!("../../components/HIS.cif"),
        "ILE" => include_str!("../../components/ILE.cif"),
        "LEU" => include_str!("../../components/LEU.cif"),
        "LYS" => include_str!("../../components/LYS.cif"),
        "MET" => include_str!("../../components/MET.cif"),
        "PHE" => include_str!("../../components/PHE.cif"),
        "PRO" => include_str!("../../components/PRO.cif"),
        "SER" => include_str!("../../components/SER.cif"),
        "THR" => include_str!("../../components/THR.cif"),
        "TRP" => include_str!("../../components/TRP.cif"),
        "TYR" => include_str!("../../components/TYR.cif"),
        "VAL" => include_str!("../../components/VAL.cif"),
        "BNZ" => include_str!("../../components/BNZ.cif"),
        other => panic!("no fixture for {other}"),
    }
}

/// The template for `code`.
pub fn template(code: &str) -> Component {
    Component::from_ccd(ccd(code)).unwrap_or_else(|e| panic!("{code}: {e}"))
}

/// All twenty-one templates.
pub fn templates() -> Vec<Component> {
    CODES.iter().map(|c| template(c)).collect()
}

/// Benzene kept; the waters, the two chloride ions and the 2-hydroxyethyl disulfide (HED, a
/// crystallisation additive) dropped.
pub fn selection() -> Selection {
    Selection::new("BNZ").dropping(&["HOH", "CL", "HED"])
}

/// 181L as a system, or the test fails saying why.
pub fn system() -> System {
    System::from_pdb(PDB_181L, &templates(), &selection()).unwrap_or_else(|e| panic!("181L: {e}"))
}

/// Columns `from` to `to` of a fixed-column line, counting from one.
pub fn columns(line: &str, from: usize, to: usize) -> &str {
    line.get(from - 1..to.min(line.len())).unwrap_or("")
}

/// `SEQRES` for chain `chain`, read here with string operations and nothing from the crate.
pub fn seqres(text: &str, chain: char) -> Vec<String> {
    text.lines()
        .filter(|l| l.starts_with("SEQRES") && columns(l, 12, 12) == chain.to_string())
        .flat_map(|l| columns(l, 20, 80).split_whitespace().map(str::to_string))
        .collect()
}

/// The atom lines of `text` (`ATOM` and `HETATM`).
pub fn atom_lines(text: &str) -> impl Iterator<Item = &str> {
    text.lines()
        .filter(|l| l.starts_with("ATOM  ") || l.starts_with("HETATM"))
}

/// One fixed-column atom record, as the format lays it out.
#[allow(clippy::too_many_arguments)]
pub fn atom_line(
    het: bool,
    serial: usize,
    name: &str,
    residue: &str,
    chain: char,
    number: i32,
    at: [f64; 3],
    element: &str,
) -> String {
    let name4 = if name.len() < 4 {
        format!(" {name}")
    } else {
        name.to_string()
    };
    format!(
        "{:<6}{:>5} {:<4}{}{:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}",
        if het { "HETATM" } else { "ATOM" },
        serial,
        name4,
        ' ',
        residue,
        chain,
        number,
        ' ',
        at[0],
        at[1],
        at[2],
        1.0,
        20.0,
        element
    )
}

/// `text` with the one line that `pick` selects replaced by `with` (`None` deletes it). Panics
/// unless exactly one line is picked, so an edit that misses its anchor cannot pass silently.
pub fn edit_one(text: &str, pick: impl Fn(&str) -> bool, with: Option<&str>) -> String {
    let hits = text.lines().filter(|l| pick(l)).count();
    assert_eq!(hits, 1, "the edit's anchor must pick exactly one line");
    let mut out = String::new();
    for line in text.lines() {
        if pick(line) {
            if let Some(w) = with {
                out.push_str(w);
                out.push('\n');
            }
        } else {
            out.push_str(line);
            out.push('\n');
        }
    }
    out
}

/// Whether atom line `l` is atom `name` of residue `number`.
pub fn is_atom(l: &str, number: i32, name: &str) -> bool {
    (l.starts_with("ATOM  ") || l.starts_with("HETATM"))
        && columns(l, 23, 26).trim() == number.to_string()
        && columns(l, 13, 16).trim() == name
}

pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

pub fn len(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// The angle between two vectors, radians, by `atan2` of the cross and dot products — accurate at
/// every angle, unlike `acos` near 0 and π.
pub fn angle_between(a: [f64; 3], b: [f64; 3]) -> f64 {
    len(cross(a, b)).atan2(dot(a, b))
}
