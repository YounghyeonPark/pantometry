//! **FreeSolv's GAFF/AM1-BCC benzene**, for W4's charge perturbation: the `.mol2` its calculated
//! hydration free energy was run from, kept here byte for byte, a strict reader of it, and the map
//! of its atoms onto the crate's by element and bonding. Data, cited, not code.
//!
//! - **The file** is `mol2files_gaff/mobley_3053621.mol2` in `mol2files_gaff.tar.gz` of FreeSolv
//!   v0.52 (D. L. Mobley and J. P. Guthrie, *J. Comput.-Aided Mol. Des.* **28**, 711 (2014); the
//!   v0.52 values G. Duarte Ramos Matos et al., *J. Chem. Eng. Data* **62**, 1559 (2017)), at the
//!   tag's commit ([`SOURCE`]), CC BY 4.0. Its SHA-256 is [`SHA256`], and [`MOBLEY_3053621`] is
//!   its text unchanged. The same commit's `gromacs/mobley_3053621.top` (in `gromacs.tar.gz`,
//!   SHA-256 `a76082492ed4845a1f207fc07fdbed90e3027415b8c623215b846896eca87a34`) carries the same
//!   twelve charges to its eight printed digits, so the file's are the ones FreeSolv's GROMACS run
//!   used; FreeSolv's v0.31 notes name this compound among those whose files' charges were repaired
//!   to match the calculation.
//! - **What it says**: C1 −0.1302, the other five carbons −0.1301, every hydrogen +0.1301. **They
//!   sum to −0.0001 e, not zero**, and C1 is 0.0001 off its ring: AM1-BCC's charges rounded to four
//!   places, as FreeSolv ran them. Its database row is −0.90 ± 0.20 by experiment and −0.81 ± 0.02
//!   calculated.
//! - **The reader** ([`parse`]) takes one molecule's `MOLECULE`, `ATOM` and `BOND` sections and
//!   refuses by name a count that is not the records', an atom or bond numbered out of turn, a row
//!   short of its columns, a charge or coordinate that is not a finite number, a bond to an atom
//!   that is not there, a repeated bond, and a second molecule. **The element is read off the
//!   atom's name** (its letters before the first digit, matched to a symbol the crate knows), since
//!   a GAFF file's type column holds GAFF's types, not SYBYL's elements.
//! - **The map** ([`charges_onto`]) is a graph isomorphism, element and bonds: never the file's
//!   order. All of them are enumerated, so a symmetric molecule's charges that are not equal by
//!   symmetry show as a spread.

#![allow(dead_code)]

use pantometry_forcefield::{Component, Element};

/// Where [`MOBLEY_3053621`] came from.
pub const SOURCE: &str = "github.com/MobleyLab/FreeSolv v0.52, commit \
                          9d8a8097d587a6c084befbb2abde6dd2a0f31ba3, mol2files_gaff.tar.gz \
                          (SHA-256 15ece6114442a8eb5e720c0ad25bb48e4a6a35392da36f72b80a6f2f415fd133) \
                          → mol2files_gaff/mobley_3053621.mol2";

/// The SHA-256 of the file, which is [`MOBLEY_3053621`]'s bytes.
pub const SHA256: &str = "42a6f4acb8bc838070cb6c1be02543fdf29bd630eded5e891383d3b4cc56a41a";

/// FreeSolv v0.52's `mobley_3053621.mol2`, benzene, byte for byte ([`SOURCE`]).
pub const MOBLEY_3053621: &str = r#"@<TRIPOS>MOLECULE
benzene
   12    12     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          1.8850   -1.0360   -0.1120 ca        1 MOL     -0.130200
      2 C2          2.9210   -1.6310    0.6080 ca        1 MOL     -0.130100
      3 C3          2.8780   -1.6520    2.0020 ca        1 MOL     -0.130100
      4 C4          1.8000   -1.0790    2.6760 ca        1 MOL     -0.130100
      5 C5          0.7640   -0.4840    1.9550 ca        1 MOL     -0.130100
      6 C6          0.8070   -0.4630    0.5610 ca        1 MOL     -0.130100
      7 H1          1.9170   -1.0200   -1.1980 ha        1 MOL      0.130100
      8 H2          3.7610   -2.0770    0.0840 ha        1 MOL      0.130100
      9 H3          3.6850   -2.1150    2.5630 ha        1 MOL      0.130100
     10 H4          1.7660   -1.0950    3.7610 ha        1 MOL      0.130100
     11 H5         -0.0760   -0.0380    2.4790 ha        1 MOL      0.130100
     12 H6          0.0000    0.0000   -0.0000 ha        1 MOL      0.130100
@<TRIPOS>BOND
     1    1    6 ar  
     2    1    2 ar  
     3    2    3 ar  
     4    3    4 ar  
     5    4    5 ar  
     6    5    6 ar  
     7    1    7 1   
     8    2    8 1   
     9    3    9 1   
    10    4   10 1   
    11    5   11 1   
    12    6   12 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's experimental hydration free energy of benzene and its uncertainty, kcal/mol: the
/// database row of [`SOURCE`]'s commit.
pub const EXPERIMENT: (f64, f64) = (-0.90, 0.20);

/// FreeSolv's calculated GAFF/AM1-BCC hydration free energy of benzene in TIP3P and its
/// uncertainty, kcal/mol: the same row.
pub const CALCULATED: (f64, f64) = (-0.81, 0.02);

/// One atom of a `.mol2`.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2Atom {
    pub name: String,
    pub element: Element,
    /// The type column as written: GAFF's type in a GAFF file.
    pub kind: String,
    /// Ångström, as written.
    pub at: [f64; 3],
    /// Elementary charges, the last column.
    pub charge: f64,
}

/// One molecule of a `.mol2`: its name, atoms and bonds, the bonds as indices into the atoms.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2 {
    pub name: String,
    pub atoms: Vec<Mol2Atom>,
    pub bonds: Vec<[usize; 2]>,
}

/// An atom name's element: the letters before its first digit, which must be a symbol the crate
/// knows, `C`, `H`, `Cl`.
pub fn element_of(name: &str) -> Result<Element, String> {
    let letters: String = name
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    Element::ALL
        .iter()
        .copied()
        .find(|e| e.symbol() == letters)
        .ok_or_else(|| format!("atom {name:?}: {letters:?} is not an element this crate knows"))
}

/// **One molecule of a `.mol2`, read strictly**: see the module documentation.
pub fn parse(text: &str) -> Result<Mol2, String> {
    let mut section = "";
    let mut molecule: Vec<&str> = Vec::new();
    let mut atom_lines: Vec<&str> = Vec::new();
    let mut bond_lines: Vec<&str> = Vec::new();
    let mut molecules = 0;
    for line in text.lines() {
        if let Some(name) = line.trim().strip_prefix("@<TRIPOS>") {
            section = match name {
                "MOLECULE" => {
                    molecules += 1;
                    if molecules > 1 {
                        return Err("a second molecule".into());
                    }
                    "MOLECULE"
                }
                "ATOM" => "ATOM",
                "BOND" => "BOND",
                _ => "other",
            };
            continue;
        }
        match section {
            "MOLECULE" => molecule.push(line),
            "ATOM" if !line.trim().is_empty() => atom_lines.push(line),
            "BOND" if !line.trim().is_empty() => bond_lines.push(line),
            _ => {}
        }
    }
    if molecules == 0 {
        return Err("no MOLECULE section".into());
    }
    let name = molecule
        .first()
        .map(|s| s.trim().to_string())
        .ok_or_else(|| "a MOLECULE section without its name".to_string())?;
    let counts: Vec<usize> = molecule
        .get(1)
        .ok_or_else(|| "a MOLECULE section without its counts".to_string())?
        .split_whitespace()
        .take(2)
        .map(|w| w.parse().map_err(|_| format!("a count {w:?}")))
        .collect::<Result<_, _>>()?;
    if counts.len() != 2 {
        return Err("a MOLECULE section without its atom and bond counts".into());
    }
    if atom_lines.len() != counts[0] {
        return Err(format!(
            "{} atoms where the counts say {}",
            atom_lines.len(),
            counts[0]
        ));
    }
    if bond_lines.len() != counts[1] {
        return Err(format!(
            "{} bonds where the counts say {}",
            bond_lines.len(),
            counts[1]
        ));
    }
    let number = |w: &str, what: &str| -> Result<f64, String> {
        match w.parse::<f64>() {
            Ok(x) if x.is_finite() => Ok(x),
            _ => Err(format!("{what} {w:?} is not a finite number")),
        }
    };
    let mut atoms = Vec::with_capacity(atom_lines.len());
    for (k, line) in atom_lines.iter().enumerate() {
        let w: Vec<&str> = line.split_whitespace().collect();
        if w.len() < 9 {
            return Err(format!(
                "atom row {} has {} columns, not the nine to its charge: {line:?}",
                k + 1,
                w.len()
            ));
        }
        if w[0] != (k + 1).to_string() {
            return Err(format!("atom row {} is numbered {}", k + 1, w[0]));
        }
        atoms.push(Mol2Atom {
            name: w[1].to_string(),
            element: element_of(w[1])?,
            kind: w[5].to_string(),
            at: [
                number(w[2], "a coordinate")?,
                number(w[3], "a coordinate")?,
                number(w[4], "a coordinate")?,
            ],
            charge: number(w[8], &format!("atom {}'s charge", w[1]))?,
        });
    }
    let mut bonds: Vec<[usize; 2]> = Vec::with_capacity(bond_lines.len());
    for (k, line) in bond_lines.iter().enumerate() {
        let w: Vec<&str> = line.split_whitespace().collect();
        if w.len() < 4 {
            return Err(format!("bond row {} is short: {line:?}", k + 1));
        }
        if w[0] != (k + 1).to_string() {
            return Err(format!("bond row {} is numbered {}", k + 1, w[0]));
        }
        let atom = |s: &str| -> Result<usize, String> {
            match s.parse::<usize>() {
                Ok(i) if (1..=atoms.len()).contains(&i) => Ok(i - 1),
                _ => Err(format!(
                    "bond {} names atom {s:?}, which is not there",
                    k + 1
                )),
            }
        };
        let (a, b) = (atom(w[1])?, atom(w[2])?);
        if a == b {
            return Err(format!("bond {} joins atom {} to itself", k + 1, a + 1));
        }
        let pair = [a.min(b), a.max(b)];
        if bonds
            .iter()
            .any(|p| [p[0].min(p[1]), p[0].max(p[1])] == pair)
        {
            return Err(format!("bond {} repeats {}–{}", k + 1, a + 1, b + 1));
        }
        bonds.push([a, b]);
    }
    Ok(Mol2 { name, atoms, bonds })
}

/// Each atom's neighbours.
fn adjacency(n: usize, bonds: &[[usize; 2]]) -> Vec<Vec<usize>> {
    let mut adj = vec![Vec::new(); n];
    for &[a, b] in bonds {
        adj[a].push(b);
        adj[b].push(a);
    }
    adj
}

/// **Every isomorphism of graph `a` onto graph `b`** — element for element, bond for bond — as
/// `map[i]`, `a`'s atom `i`'s atom of `b`; in the order a search over `a`'s atoms in turn, each
/// trying `b`'s in turn, finds them. Empty when the two are not the same molecule.
pub fn isomorphisms(
    a: (&[Element], &[[usize; 2]]),
    b: (&[Element], &[[usize; 2]]),
) -> Vec<Vec<usize>> {
    let n = a.0.len();
    if b.0.len() != n || a.1.len() != b.1.len() {
        return Vec::new();
    }
    let mut search = Search {
        a: a.0,
        b: b.0,
        adj_a: adjacency(n, a.1),
        adj_b: adjacency(n, b.1),
        map: vec![usize::MAX; n],
        used: vec![false; n],
        out: Vec::new(),
    };
    search.extend(0);
    search.out
}

/// [`isomorphisms`]' backtracking search: the two graphs, the map so far, which of `b`'s atoms it
/// uses, and the whole maps found.
struct Search<'a> {
    a: &'a [Element],
    b: &'a [Element],
    adj_a: Vec<Vec<usize>>,
    adj_b: Vec<Vec<usize>>,
    map: Vec<usize>,
    used: Vec<bool>,
    out: Vec<Vec<usize>>,
}

impl Search<'_> {
    /// Every whole map that extends the one placing `a`'s atoms before `k`.
    fn extend(&mut self, k: usize) {
        if k == self.a.len() {
            self.out.push(self.map.clone());
            return;
        }
        for c in 0..self.b.len() {
            if self.used[c] || self.b[c] != self.a[k] || self.adj_b[c].len() != self.adj_a[k].len()
            {
                continue;
            }
            // Every bond of `k` to an atom already placed is a bond of `c`: with the degrees and
            // the bond counts equal, a whole map is then a bijection of the bonds.
            if self.adj_a[k]
                .iter()
                .filter(|&&j| j < k)
                .all(|&j| self.adj_b[c].contains(&self.map[j]))
            {
                self.map[k] = c;
                self.used[c] = true;
                self.extend(k + 1);
                self.used[c] = false;
                self.map[k] = usize::MAX;
            }
        }
    }
}

/// A `.mol2`'s charges on another description of the same molecule.
#[derive(Clone, Debug)]
pub struct Mapped {
    /// The charge of each atom of the target, in the target's order.
    pub charges: Vec<f64>,
    /// The map used, the first [`isomorphisms`] finds: the target's atom `i` is the file's
    /// `map[i]`.
    pub map: Vec<usize>,
    /// How many isomorphisms there are: the molecule's automorphisms.
    pub automorphisms: usize,
    /// The largest difference, over the target's atoms, between the charges two isomorphisms give
    /// one atom: zero when the file's charges are equal by symmetry.
    pub spread: f64,
}

/// **`mol2`'s charges on the molecule whose atoms are `elements` and bonds `bonds`**, by element
/// and bonding ([`isomorphisms`]); refused when they are not the same molecule.
pub fn charges_onto(
    mol2: &Mol2,
    elements: &[Element],
    bonds: &[[usize; 2]],
) -> Result<Mapped, String> {
    let file: Vec<Element> = mol2.atoms.iter().map(|a| a.element).collect();
    let all = isomorphisms((elements, bonds), (&file, &mol2.bonds));
    let map = all
        .first()
        .ok_or_else(|| format!("{} is not this molecule", mol2.name))?
        .clone();
    let mut spread: f64 = 0.0;
    for i in 0..elements.len() {
        let q = all.iter().map(|m| mol2.atoms[m[i]].charge);
        let (lo, hi) = q.fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), x| {
            (l.min(x), h.max(x))
        });
        spread = spread.max(hi - lo);
    }
    Ok(Mapped {
        charges: map.iter().map(|&j| mol2.atoms[j].charge).collect(),
        map,
        automorphisms: all.len(),
        spread,
    })
}

/// A component's elements and bonds, as [`charges_onto`] takes them.
pub fn component_graph(c: &Component) -> (Vec<Element>, Vec<[usize; 2]>) {
    (
        c.atoms().iter().map(|a| a.element).collect(),
        c.bonds().iter().map(|b| b.atoms).collect(),
    )
}

/// FreeSolv's benzene, read ([`MOBLEY_3053621`]).
pub fn benzene() -> Mol2 {
    parse(MOBLEY_3053621).expect("FreeSolv's mobley_3053621.mol2")
}
