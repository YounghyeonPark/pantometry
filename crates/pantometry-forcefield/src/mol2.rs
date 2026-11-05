//! One molecule of a Tripos `.mol2`, read strictly, and its atoms mapped onto a [`Component`]'s.
//!
//! A `.mol2` is how a parameterised small molecule usually leaves the tool that parameterised it:
//! antechamber writes one, FreeSolv distributes one per molecule, and each atom row carries a
//! **type** and a **partial charge** — GAFF's type and AM1-BCC's charge, in the files this crate has been
//! measured with. This module reads those two columns and puts them on the crate's own description
//! of the same molecule. The bonded terms stay UFF's, whatever the file says; see
//! [`crate::amber`] for what a model built from it is.
//!
//! # What is read, and what is refused
//!
//! One molecule's `MOLECULE`, `ATOM` and `BOND` sections; every other section is skipped. **Every
//! refusal is a [`Mol2Error`] naming what it refused**:
//! - a count in `MOLECULE` that is not the number of records;
//! - an atom or bond numbered out of turn;
//! - an atom row short of the nine columns that end in its charge;
//! - a charge or coordinate that is not a finite number;
//! - an atom name whose element is not one the crate knows;
//! - a bond to an atom that is not there, to itself, or given twice;
//! - a second molecule.
//!
//! **The element is read off the atom's name**, by [`element_of`]: its leading letters matched to
//! a symbol the crate knows, without regard to case. A GAFF file's type column holds GAFF's types,
//! not SYBYL's elements, so the type column cannot say what the atom is.
//!
//! # The map is a graph isomorphism, never the file's order
//!
//! [`Mol2::onto`] matches the file's atoms to a component's by element and bonds, and **finds
//! every match**. A file writes its atoms in whatever order its tool chose, so mapping by position
//! would put a hydrogen's charge on a carbon and look fine. A symmetric molecule has several
//! matches, one per automorphism. When the file's charges are not equal by that symmetry, the
//! difference between what two matches give one atom shows as [`Mapped::spread`], so a reader can
//! see the file's rounding and the choice of match is not hidden.
//!
//! The search is plain backtracking, by element and degree, and it visits every match without
//! keeping them. The number of matches is the molecule's symmetry: 12 for benzene, but
//! 2 · (3!)² · 2¹⁴ = 1 179 648 for hexadecane, whose fourteen CH₂ groups each swap their two
//! hydrogens. Hexadecane took 0.14 s in release on the machine this was written on, in a few
//! megabytes. Every further CH₂ doubles the time, and the memory does not grow. It is not meant for
//! a protein.
//!
//! # Why the parser takes text and not a path
//!
//! As [`crate::ccd`]: nothing in `crates/` opens a file.

use std::fmt;

use crate::ccd::{Component, Element};

/// One atom of a `.mol2`.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2Atom {
    /// The atom name, as written.
    pub name: String,
    /// Its element, read off the name.
    pub element: Element,
    /// The type column as written: GAFF's type in a GAFF file.
    pub kind: String,
    /// Ångström, as written.
    pub at: [f64; 3],
    /// Elementary charges, the ninth column.
    pub charge: f64,
}

/// One molecule of a `.mol2`: its name, atoms and bonds, the bonds as indices into the atoms in
/// the order the file writes them.
#[derive(Clone, Debug, PartialEq)]
pub struct Mol2 {
    /// The first line of the `MOLECULE` section.
    pub name: String,
    /// The atoms, in the file's order.
    pub atoms: Vec<Mol2Atom>,
    /// The bonds, each as the file writes it.
    pub bonds: Vec<[usize; 2]>,
}

/// Why a `.mol2` did not produce a molecule, or did not map onto one. Every variant names what it
/// refused.
#[derive(Clone, Debug, PartialEq)]
pub enum Mol2Error {
    /// No `@<TRIPOS>MOLECULE` section.
    NoMolecule,
    /// A second `@<TRIPOS>MOLECULE` section.
    SecondMolecule,
    /// The `MOLECULE` section has no name line, or no line with its atom and bond counts.
    NoCounts,
    /// A count that is not a number.
    Count {
        /// What it said.
        text: String,
    },
    /// The number of atom or bond records is not the section's count.
    Records {
        /// `"atoms"` or `"bonds"`.
        what: &'static str,
        /// How many rows there are.
        found: usize,
        /// How many the counts say.
        counted: usize,
    },
    /// An atom row has fewer columns than the nine that end in its charge.
    Columns {
        /// The row, counting from one.
        row: usize,
        /// How many columns it has.
        columns: usize,
        /// The line.
        line: String,
    },
    /// A bond row is shorter than its number, two atoms and an order.
    ShortBond {
        /// The row, counting from one.
        row: usize,
        /// The line.
        line: String,
    },
    /// A row numbered other than its place.
    Numbered {
        /// `"atom"` or `"bond"`.
        what: &'static str,
        /// The row, counting from one.
        row: usize,
        /// The number it carries.
        number: String,
    },
    /// An atom name whose leading letters are not an element the crate knows.
    Element {
        /// The atom name.
        name: String,
    },
    /// A coordinate or charge that is not a finite number.
    NotFinite {
        /// What it is: `"a coordinate"` or `"atom C1's charge"`.
        what: String,
        /// What it said.
        text: String,
    },
    /// A bond names an atom that is not there.
    NoSuchAtom {
        /// The bond, counting from one.
        bond: usize,
        /// What it said.
        text: String,
    },
    /// A bond joins an atom to itself.
    SelfBond {
        /// The bond, counting from one.
        bond: usize,
        /// The atom, counting from one.
        atom: usize,
    },
    /// A bond repeats an earlier one.
    RepeatedBond {
        /// The bond, counting from one.
        bond: usize,
        /// Its atoms, counting from one.
        atoms: [usize; 2],
    },
    /// [`Mol2::onto`]: no match of the file's molecule onto the target, by element and bonds.
    NotThisMolecule {
        /// The file's molecule name.
        name: String,
        /// The file's composition, `C6H6`.
        file: String,
        /// The target's.
        target: String,
    },
}

impl fmt::Display for Mol2Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Mol2Error::NoMolecule => write!(f, "no MOLECULE section"),
            Mol2Error::SecondMolecule => write!(f, "a second molecule"),
            Mol2Error::NoCounts => write!(
                f,
                "a MOLECULE section without its name and its atom and bond counts"
            ),
            Mol2Error::Count { text } => write!(f, "a count {text:?}"),
            Mol2Error::Records {
                what,
                found,
                counted,
            } => write!(f, "{found} {what} where the counts say {counted}"),
            Mol2Error::Columns { row, columns, line } => write!(
                f,
                "atom row {row} has {columns} columns, not the nine to its charge: {line:?}"
            ),
            Mol2Error::ShortBond { row, line } => write!(f, "bond row {row} is short: {line:?}"),
            Mol2Error::Numbered { what, row, number } => {
                write!(f, "{what} row {row} is numbered {number}")
            }
            Mol2Error::Element { name } => {
                write!(
                    f,
                    "atom {name:?}: its leading letters are not an element this crate knows, nor \
                     begin with one (name each atom by its element, as C1, Cl2 or HO)"
                )
            }
            Mol2Error::NotFinite { what, text } => {
                write!(f, "{what} {text:?} is not a finite number")
            }
            Mol2Error::NoSuchAtom { bond, text } => {
                write!(f, "bond {bond} names atom {text:?}, which is not there")
            }
            Mol2Error::SelfBond { bond, atom } => {
                write!(f, "bond {bond} joins atom {atom} to itself")
            }
            Mol2Error::RepeatedBond { bond, atoms } => {
                write!(f, "bond {bond} repeats {}–{}", atoms[0], atoms[1])
            }
            Mol2Error::NotThisMolecule { name, file, target } if file != target => write!(
                f,
                "{name} is not this molecule: the file is {file} and the molecule {target}"
            ),
            Mol2Error::NotThisMolecule { name, .. } => write!(
                f,
                "{name} is not this molecule: the same atoms, bonded differently"
            ),
        }
    }
}

impl std::error::Error for Mol2Error {}

/// An atom name's element, read off its letters: those after any leading digits and before the
/// next character that is not a letter, compared without regard to case. Taken in turn:
/// - all of them, if they are a symbol the crate knows: `C1`, `Cl1`, `CL1`, `c1`;
/// - otherwise their first two, if those are a two-letter symbol: `CLA1` is chlorine;
/// - otherwise their first one, if it is a one-letter symbol: `HO`, `CA`, `OAB`, `1H2`.
///
/// So `CA` is carbon, as a PDB-derived name means, never calcium, which this crate has no
/// parameters for. A name read as the wrong element is not silent: [`Mol2::onto`] then finds no
/// match.
///
/// # Errors
///
/// [`Mol2Error::Element`] when none of the three is a symbol the crate knows.
pub fn element_of(name: &str) -> Result<Element, Mol2Error> {
    let letters: String = name
        .chars()
        .skip_while(|c| c.is_ascii_digit())
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    let is = |e: &Element, s: &str| e.symbol().eq_ignore_ascii_case(s);
    let whole = Element::ALL.iter().find(|e| is(e, &letters));
    let two = || {
        letters.get(..2).and_then(|s| {
            Element::ALL
                .iter()
                .find(|e| e.symbol().len() == 2 && is(e, s))
        })
    };
    let one = || {
        letters.get(..1).and_then(|s| {
            Element::ALL
                .iter()
                .find(|e| e.symbol().len() == 1 && is(e, s))
        })
    };
    whole
        .or_else(two)
        .or_else(one)
        .copied()
        .ok_or_else(|| Mol2Error::Element { name: name.into() })
}

impl Mol2 {
    /// **One molecule of a `.mol2`, read strictly**: see the module documentation.
    ///
    /// # Errors
    ///
    /// [`Mol2Error`], naming what it refused.
    pub fn parse(text: &str) -> Result<Mol2, Mol2Error> {
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
                            return Err(Mol2Error::SecondMolecule);
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
            return Err(Mol2Error::NoMolecule);
        }
        let name = molecule
            .first()
            .map(|s| s.trim().to_string())
            .ok_or(Mol2Error::NoCounts)?;
        let counts: Vec<usize> = molecule
            .get(1)
            .ok_or(Mol2Error::NoCounts)?
            .split_whitespace()
            .take(2)
            .map(|w| w.parse().map_err(|_| Mol2Error::Count { text: w.into() }))
            .collect::<Result<_, _>>()?;
        if counts.len() != 2 {
            return Err(Mol2Error::NoCounts);
        }
        for (what, found, counted) in [
            ("atoms", atom_lines.len(), counts[0]),
            ("bonds", bond_lines.len(), counts[1]),
        ] {
            if found != counted {
                return Err(Mol2Error::Records {
                    what,
                    found,
                    counted,
                });
            }
        }
        let number = |w: &str, what: &str| -> Result<f64, Mol2Error> {
            match w.parse::<f64>() {
                Ok(x) if x.is_finite() => Ok(x),
                _ => Err(Mol2Error::NotFinite {
                    what: what.into(),
                    text: w.into(),
                }),
            }
        };
        let mut atoms = Vec::with_capacity(atom_lines.len());
        for (k, line) in atom_lines.iter().enumerate() {
            let w: Vec<&str> = line.split_whitespace().collect();
            if w.len() < 9 {
                return Err(Mol2Error::Columns {
                    row: k + 1,
                    columns: w.len(),
                    line: (*line).into(),
                });
            }
            if w[0] != (k + 1).to_string() {
                return Err(Mol2Error::Numbered {
                    what: "atom",
                    row: k + 1,
                    number: w[0].into(),
                });
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
                return Err(Mol2Error::ShortBond {
                    row: k + 1,
                    line: (*line).into(),
                });
            }
            if w[0] != (k + 1).to_string() {
                return Err(Mol2Error::Numbered {
                    what: "bond",
                    row: k + 1,
                    number: w[0].into(),
                });
            }
            let atom = |s: &str| -> Result<usize, Mol2Error> {
                match s.parse::<usize>() {
                    Ok(i) if (1..=atoms.len()).contains(&i) => Ok(i - 1),
                    _ => Err(Mol2Error::NoSuchAtom {
                        bond: k + 1,
                        text: s.into(),
                    }),
                }
            };
            let (a, b) = (atom(w[1])?, atom(w[2])?);
            if a == b {
                return Err(Mol2Error::SelfBond {
                    bond: k + 1,
                    atom: a + 1,
                });
            }
            let pair = [a.min(b), a.max(b)];
            if bonds
                .iter()
                .any(|p| [p[0].min(p[1]), p[0].max(p[1])] == pair)
            {
                return Err(Mol2Error::RepeatedBond {
                    bond: k + 1,
                    atoms: [a + 1, b + 1],
                });
            }
            bonds.push([a, b]);
        }
        Ok(Mol2 { name, atoms, bonds })
    }

    /// **This file's atoms on `c`'s**, by element and bonding: see the module documentation.
    ///
    /// # Errors
    ///
    /// [`Mol2Error::NotThisMolecule`] when `c` is not the file's molecule.
    pub fn onto(&self, c: &Component) -> Result<Mapped, Mol2Error> {
        let (elements, bonds) = graph(c);
        self.onto_graph(&elements, &bonds)
    }

    /// **This file's atoms on the molecule whose atoms are `elements` and bonds `bonds`**, as
    /// [`Mol2::onto`] does for a component: for a molecule that is part of something larger, such
    /// as a ligand in a complex.
    ///
    /// # Errors
    ///
    /// [`Mol2Error::NotThisMolecule`] when they are not the file's molecule.
    pub fn onto_graph(
        &self,
        elements: &[Element],
        bonds: &[[usize; 2]],
    ) -> Result<Mapped, Mol2Error> {
        let file: Vec<Element> = self.atoms.iter().map(|a| a.element).collect();
        // Each match is seen once and not kept: the first, a count, and each atom's lowest and
        // highest charge over all of them.
        let n = elements.len();
        let mut first: Option<Vec<usize>> = None;
        let mut count = 0usize;
        let mut lo = vec![f64::INFINITY; n];
        let mut hi = vec![f64::NEG_INFINITY; n];
        each_isomorphism(
            (elements, bonds),
            (&file, &self.bonds),
            &mut |m: &[usize]| {
                if first.is_none() {
                    first = Some(m.to_vec());
                }
                count += 1;
                for i in 0..n {
                    let q = self.atoms[m[i]].charge;
                    lo[i] = lo[i].min(q);
                    hi[i] = hi[i].max(q);
                }
            },
        );
        let map = first.ok_or_else(|| Mol2Error::NotThisMolecule {
            name: self.name.clone(),
            file: composition(&file),
            target: composition(elements),
        })?;
        let spread = lo.iter().zip(&hi).fold(0.0f64, |s, (l, h)| s.max(h - l));
        Ok(Mapped {
            charges: map.iter().map(|&j| self.atoms[j].charge).collect(),
            kinds: map.iter().map(|&j| self.atoms[j].kind.clone()).collect(),
            map,
            automorphisms: count,
            spread,
        })
    }
}

/// A `.mol2`'s atoms on another description of the same molecule ([`Mol2::onto`]).
#[derive(Clone, Debug, PartialEq)]
pub struct Mapped {
    /// The charge of each atom of the target, in the target's order.
    pub charges: Vec<f64>,
    /// The type column of each atom of the target, in the target's order.
    pub kinds: Vec<String>,
    /// The match used, the first [`isomorphisms`] finds: the target's atom `i` is the file's
    /// `map[i]`.
    pub map: Vec<usize>,
    /// How many matches there are: the molecule's automorphisms.
    pub automorphisms: usize,
    /// The largest difference, over the target's atoms, between the charges two matches give one
    /// atom: zero when the file's charges are equal by symmetry.
    pub spread: f64,
}

/// A component's elements and bonds, as [`Mol2::onto_graph`] takes them.
pub fn graph(c: &Component) -> (Vec<Element>, Vec<[usize; 2]>) {
    (
        c.atoms().iter().map(|a| a.element).collect(),
        c.bonds().iter().map(|b| b.atoms).collect(),
    )
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

/// **Every isomorphism of graph `a` onto graph `b`**, element for element and bond for bond, as
/// `map[i]`: `a`'s atom `i` is `b`'s atom `map[i]`. They come in the order a search finds them,
/// taking `a`'s atoms in turn and trying `b`'s in turn for each. Empty when the two are not the
/// same molecule.
///
/// # Panics
///
/// If a bond names an atom past the end of its graph's elements.
pub fn isomorphisms(
    a: (&[Element], &[[usize; 2]]),
    b: (&[Element], &[[usize; 2]]),
) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    each_isomorphism(a, b, &mut |m: &[usize]| out.push(m.to_vec()));
    out
}

/// [`isomorphisms`], each handed to `found` as the search finds it and none kept.
fn each_isomorphism(
    a: (&[Element], &[[usize; 2]]),
    b: (&[Element], &[[usize; 2]]),
    found: &mut dyn FnMut(&[usize]),
) {
    let n = a.0.len();
    if b.0.len() != n || a.1.len() != b.1.len() {
        return;
    }
    let mut search = Search {
        a: a.0,
        b: b.0,
        adj_a: adjacency(n, a.1),
        adj_b: adjacency(n, b.1),
        map: vec![usize::MAX; n],
        used: vec![false; n],
        found,
    };
    search.extend(0);
}

/// A molecule's composition in Hill's order, `C6H5Cl`: carbon, hydrogen, then the rest by symbol,
/// counts of one left out.
fn composition(elements: &[Element]) -> String {
    let mut order: Vec<Element> = Element::ALL
        .iter()
        .copied()
        .filter(|&e| e != Element::C && e != Element::H)
        .collect();
    order.sort_by_key(|e| e.symbol());
    order.insert(0, Element::H);
    order.insert(0, Element::C);
    let mut s = String::new();
    for e in order {
        match elements.iter().filter(|&&x| x == e).count() {
            0 => {}
            1 => s += e.symbol(),
            k => s += &format!("{}{k}", e.symbol()),
        }
    }
    s
}

/// [`isomorphisms`]' backtracking search: the two graphs, the map so far, which of `b`'s atoms it
/// uses, and what to do with each whole map.
struct Search<'a> {
    a: &'a [Element],
    b: &'a [Element],
    adj_a: Vec<Vec<usize>>,
    adj_b: Vec<Vec<usize>>,
    map: Vec<usize>,
    used: Vec<bool>,
    found: &'a mut dyn FnMut(&[usize]),
}

impl Search<'_> {
    /// Every whole map that extends the one placing `a`'s atoms before `k`.
    fn extend(&mut self, k: usize) {
        if k == self.a.len() {
            (self.found)(&self.map);
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
