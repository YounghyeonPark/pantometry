//! A protein's backbone, read out of a PDB file.
//!
//! One node per residue, at its alpha carbon. That is the whole model: no side chains, no
//! hydrogens, no solvent. A coarse-grained network model is not an approximation to an all-atom
//! one that got cheaper — it is a different claim, that the **collective** motions of a protein
//! are set by the shape of the fold and not by its chemistry, and the evidence for that claim is
//! that the motions come out right anyway. [`Structure::experimental_fluctuations`] is what makes
//! that checkable here rather than asserted.
//!
//! # Why the parser takes text and not a path
//!
//! Nothing in `crates/` opens a file. Every crate in this workspace compiles to `wasm32`, where
//! there is no filesystem, and the one that reads the scenes off a disk lives in `app/` for
//! exactly this reason. So [`Structure::from_pdb`] takes the contents and the caller finds them.
//!
//! # What the file format does that a reader has to get right
//!
//! The PDB format is fixed-column, which is a source of two silent failures rather than of
//! loud ones:
//!
//! - **`CA` is both an alpha carbon and a calcium ion.** Both write `CA` in the atom-name field.
//!   They are told apart by the record type — `ATOM` against `HETATM` — and by the element
//!   column, and a reader that keys on the name alone silently gains a node wherever a structure
//!   binds calcium.
//! - **Alternate locations and multiple models repeat the same residue.** A residue with two
//!   modelled side-chain positions appears twice with `altLoc` `A` and `B`, and an NMR entry
//!   repeats every atom once per model. A reader that takes them all builds a network with
//!   duplicate nodes at nearly the same place, whose Hessian has a near-zero eigenvalue that
//!   looks exactly like a soft collective mode.
//!
//! Both are handled, and both have a test that fails when the handling is removed.

use crate::spectrum::Symmetric;
use std::fmt;

/// Ångström, in metres. PDB coordinates and B-factors are in these; everything past the parser
/// is SI.
const ANGSTROM: f64 = 1e-10;

/// One residue: where its alpha carbon is, and what the crystallographer measured there.
#[derive(Clone, Debug, PartialEq)]
pub struct Residue {
    /// The three-letter name, `"ALA"`, `"GLY"`. Kept as read, uppercase, not validated against a
    /// table of twenty — a structure may legitimately contain a modified residue.
    pub name: String,
    /// The chain it belongs to, as one character. A complex of two chains is one structure with
    /// two values here, which is what makes an interface an interface.
    pub chain: char,
    /// The residue number from the file. Not an index: numbering starts where the
    /// crystallographer says and may have gaps where residues were not resolved.
    pub number: i32,
    /// Where the alpha carbon is, in metres.
    pub at: [f64; 3],
    /// The temperature factor, converted from `Å²` to `m²`.
    ///
    /// `B = 8π²⟨Δr²⟩/3`, so this is not itself a mean-square displacement — see
    /// [`Structure::experimental_fluctuations`], which does that division.
    pub b_factor: f64,
}

/// Why a PDB file did not produce a structure.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// No `ATOM` record with an alpha carbon anywhere in the text.
    ///
    /// The failure this exists for is the quiet one: a file of the wrong kind, or a file that is
    /// all `HETATM`, parses without complaint into an empty list, and an empty network has an
    /// empty spectrum and no error anywhere.
    NoAlphaCarbons,
    /// A line began `ATOM` and was too short to hold coordinates.
    ShortRecord {
        /// Which line, counting from one.
        line: usize,
    },
    /// A coordinate or B-factor column did not hold a number.
    NotANumber {
        /// Which line, counting from one.
        line: usize,
        /// What it said instead.
        text: String,
    },
    /// Two residues landed at the same place, to within a picometre.
    ///
    /// A duplicate node makes the network singular in a way that is indistinguishable from
    /// physics: the pair's antisymmetric motion costs nothing, so the spectrum gains a
    /// near-zero eigenvalue and the lowest "collective mode" is two atoms sliding through each
    /// other.
    Coincident {
        /// The first of the pair, by index.
        first: usize,
        /// The second.
        second: usize,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::NoAlphaCarbons => write!(
                f,
                "no ATOM record with an alpha carbon; a file with none parses to an empty \
                 structure, which has no spectrum and no error"
            ),
            ParseError::ShortRecord { line } => {
                write!(
                    f,
                    "line {line} begins ATOM and is too short to hold coordinates"
                )
            }
            ParseError::NotANumber { line, text } => {
                write!(f, "line {line}: {text:?} is not a number")
            }
            ParseError::Coincident { first, second } => write!(
                f,
                "residues {first} and {second} are at the same place, which makes the network \
                 singular in a way that reads as a soft mode"
            ),
        }
    }
}

impl std::error::Error for ParseError {}

/// One atom of a ligand: where it is, and which element the file says it is.
///
/// Not a [`Residue`], which is one alpha carbon standing for a whole amino acid. A ligand is
/// modelled atom by atom, and the two are different things however similar the fields look.
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// Where it is, in metres.
    pub at: [f64; 3],
    /// The element, upper case, as columns 77-78 of the record spell it: `C`, `N`, `O`, `P`.
    ///
    /// **Taken from that column and not from the atom's name**, because the name is ambiguous in
    /// exactly the cases that matter. `PA` in `1AKE` is the alpha *phosphorus* of a phosphate and
    /// not protactinium; `CA` is a calcium ion in one file and an alpha carbon in the next. The
    /// element column exists because the name cannot be read, and a reader that guesses from the
    /// name is a reader that will be wrong about a metal one day and not say so.
    ///
    /// Empty when the record is too short to have the column, which is a pre-1996 file or a
    /// hand-written one. A caller drawing by element gets to decide what that means; this will
    /// not invent one.
    pub element: String,
}

/// A protein as a list of alpha carbons.
#[derive(Clone, Debug, PartialEq)]
pub struct Structure {
    residues: Vec<Residue>,
}

impl Structure {
    /// The atoms of one heteroatom group, for a ligand that came in the same file.
    ///
    /// [`Structure::from_pdb`] refuses `HETATM` and says why: a calcium ion writes `CA` in the
    /// atom-name column and is not an alpha carbon. A bound ligand is in those same records, and
    /// [`Network::with_ligand`](crate::Network::with_ligand) wants its coordinates -- so this is
    /// the one way from a crystal structure to both halves of the complex it holds.
    ///
    /// `name` is the residue name in columns 18-20: `AP5` for the bi-substrate inhibitor in
    /// `1AKE`, `HOH` for water. **Every atom of every copy**, so a file with the complex twice in
    /// its asymmetric unit gives two ligands' worth -- see [`Structure::chain`] for the same
    /// problem on the protein side, and pass `chain` here rather than filtering afterwards.
    ///
    /// Hydrogens are not special-cased. A crystal structure at ordinary resolution has none, and
    /// one that does is a file whose author meant them to be there.
    ///
    /// Each atom carries its [`element`](Atom::element), read from columns 77-78. It used to
    /// return bare coordinates, and a caller drawing the molecule could then only draw every atom
    /// the same size — which for `AP5A` means five phosphorus atoms the size of a carbon, a
    /// picture asserting something about the molecule that the data does not.
    ///
    /// Empty when nothing matched, which is a real answer: an apo structure has no ligand, and
    /// [`Network::with_ligand`](crate::Network::with_ligand) of nothing is the network it started
    /// as. A caller that asked for a group it expected to find should say so itself -- this
    /// cannot tell a typo from an apo form.
    pub fn ligand_from_pdb(text: &str, name: &str, chain: Option<char>) -> Vec<Atom> {
        let mut atoms: Vec<Atom> = Vec::new();
        let mut seen_alt: Vec<(char, i32, String, char)> = Vec::new();
        for line in text.lines() {
            if line.starts_with("ENDMDL") {
                break;
            }
            if !line.starts_with("HETATM") || line.len() < 54 {
                continue;
            }
            if line[17..20].trim() != name {
                continue;
            }
            let at = line[21..22].chars().next().unwrap_or(' ');
            if chain.is_some_and(|c| c != at) {
                continue;
            }
            // The same alternate-location rule the protein side uses, for the same reason: a
            // second conformer of one atom is that atom modelled twice, not a second atom.
            let alt = line[16..17].chars().next().unwrap_or(' ');
            let seq: i32 = line[22..26].trim().parse().unwrap_or(0);
            let atom = line[12..16].trim().to_string();
            let key = (at, seq, atom, alt);
            let first = (key.0, key.1, key.2.clone(), ' ');
            if seen_alt.contains(&first) || seen_alt.contains(&key) {
                continue;
            }
            seen_alt.push(first);
            let parse = |from: usize, to: usize| line[from..to].trim().parse::<f64>().ok();
            if let (Some(x), Some(y), Some(z)) = (parse(30, 38), parse(38, 46), parse(46, 54)) {
                // Columns 77-78, and nothing else. A record too short to hold them gives an empty
                // string rather than a guess off the atom name -- see `Atom::element`.
                let element = line
                    .get(76..78)
                    .map(|e| e.trim().to_ascii_uppercase())
                    .unwrap_or_default();
                // Angstrom in the file, metres here, as everywhere in this crate.
                atoms.push(Atom {
                    at: [x * 1e-10, y * 1e-10, z * 1e-10],
                    element,
                });
            }
        }
        atoms
    }

    /// The alpha carbons of the first model in `text`.
    ///
    /// Takes `ATOM` records whose atom name is `CA` and whose element column is blank or `C`,
    /// keeping the first alternate location of each; stops at the first `ENDMDL`, so an NMR
    /// ensemble gives its first model rather than every model's atoms piled on top of each other.
    ///
    /// # Errors
    ///
    /// [`ParseError::NoAlphaCarbons`] if nothing matched — the case that would otherwise succeed
    /// quietly and produce an empty spectrum. [`ParseError::Coincident`] if two of them are at
    /// the same place. The two parse errors carry the line.
    pub fn from_pdb(text: &str) -> Result<Structure, ParseError> {
        let mut residues = Vec::new();
        let mut seen_alt: Vec<(char, i32, char)> = Vec::new();
        for (index, line) in text.lines().enumerate() {
            let number = index + 1;
            if line.starts_with("ENDMDL") {
                break;
            }
            // `HETATM` is deliberately not here: a calcium ion writes `CA` in the atom-name
            // column and is not an alpha carbon.
            if !line.starts_with("ATOM") {
                continue;
            }
            let bytes = line.as_bytes();
            if bytes.len() < 54 {
                return Err(ParseError::ShortRecord { line: number });
            }
            if line[12..16].trim() != "CA" {
                continue;
            }
            // The element column, where a file has one. Calcium says `CA` here and carbon says
            // `C`, which is the only column that tells them apart in an `ATOM` record.
            let element = if bytes.len() >= 78 {
                line[76..78].trim().to_string()
            } else {
                String::new()
            };
            if !element.is_empty() && element != "C" {
                continue;
            }

            let alt = line[16..17].chars().next().unwrap_or(' ');
            let chain = line[21..22].chars().next().unwrap_or(' ');
            let seq: i32 = match line[22..26].trim().parse() {
                Ok(v) => v,
                Err(_) => {
                    return Err(ParseError::NotANumber {
                        line: number,
                        text: line[22..26].to_string(),
                    })
                }
            };
            // One residue, one node. A second alternate location for a residue already taken is
            // the same residue modelled twice, not a second residue.
            if seen_alt.contains(&(chain, seq, ' ')) || seen_alt.contains(&(chain, seq, alt)) {
                continue;
            }
            seen_alt.push((chain, seq, ' '));

            let mut coordinate = [0.0; 3];
            for (k, span) in [30..38, 38..46, 46..54].into_iter().enumerate() {
                coordinate[k] = match line[span.clone()].trim().parse::<f64>() {
                    Ok(v) => v * ANGSTROM,
                    Err(_) => {
                        return Err(ParseError::NotANumber {
                            line: number,
                            text: line[span].to_string(),
                        })
                    }
                };
            }
            // A file without the B-factor column is a structure without experimental
            // fluctuations, not a parse failure: the model still runs, and
            // `experimental_fluctuations` is what refuses.
            let b_factor = if bytes.len() >= 66 {
                line[60..66].trim().parse::<f64>().unwrap_or(0.0) * ANGSTROM * ANGSTROM
            } else {
                0.0
            };

            residues.push(Residue {
                name: line[17..20].trim().to_string(),
                chain,
                number: seq,
                at: coordinate,
                b_factor,
            });
        }

        if residues.is_empty() {
            return Err(ParseError::NoAlphaCarbons);
        }
        // A picometre is a hundredth of the last digit a PDB file writes, so this fires on a
        // genuine duplicate and never on rounding.
        let close = 1e-12;
        for i in 0..residues.len() {
            for j in (i + 1)..residues.len() {
                let d: f64 = (0..3)
                    .map(|k| (residues[i].at[k] - residues[j].at[k]).powi(2))
                    .sum();
                if d < close * close {
                    return Err(ParseError::Coincident {
                        first: i,
                        second: j,
                    });
                }
            }
        }
        Ok(Structure { residues })
    }

    /// A structure from positions in metres, for a test that wants a shape rather than a protein.
    ///
    /// # Panics
    ///
    /// If `positions` is empty.
    pub fn from_positions(positions: Vec<[f64; 3]>) -> Structure {
        assert!(!positions.is_empty(), "a structure with no nodes");
        Structure {
            residues: positions
                .into_iter()
                .enumerate()
                .map(|(i, at)| Residue {
                    name: "UNK".to_string(),
                    chain: 'A',
                    number: i as i32 + 1,
                    at,
                    b_factor: 0.0,
                })
                .collect(),
        }
    }

    /// How many residues.
    pub fn len(&self) -> usize {
        self.residues.len()
    }

    /// Whether there are none. There never are — [`Structure::from_pdb`] refuses an empty
    /// structure and [`Structure::from_positions`] panics on one — so this is always false, and
    /// it is here because `len` without `is_empty` is a clippy lint and because a caller should
    /// not have to know that.
    pub fn is_empty(&self) -> bool {
        self.residues.is_empty()
    }

    /// The residues, in file order.
    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    /// Only the residues of one chain.
    ///
    /// A deposited file often holds more than one copy of the protein — both adenylate kinase
    /// entries in this crate's `structures/` are dimers — and the second copy is
    /// crystallography rather than biology. Building a network across both joins two molecules
    /// that are not joined, and everything downstream is about that dimer.
    ///
    /// `None` if no residue has that chain, rather than an empty structure.
    pub fn chain(&self, id: char) -> Option<Structure> {
        let residues: Vec<Residue> = self
            .residues
            .iter()
            .filter(|r| r.chain == id)
            .cloned()
            .collect();
        if residues.is_empty() {
            None
        } else {
            Some(Structure { residues })
        }
    }

    /// Root-mean-square distance between matching residues, with nothing fitted.
    ///
    /// `None` if the two are different lengths. This pairs residue `i` with residue `i` and does
    /// not align sequences: two structures compared this way have to be the same protein, read
    /// the same way.
    pub fn rmsd_to(&self, other: &Structure) -> Option<f64> {
        if self.residues.len() != other.residues.len() {
            return None;
        }
        Some(
            (self
                .residues
                .iter()
                .zip(&other.residues)
                .map(|(a, b)| (0..3).map(|k| (a.at[k] - b.at[k]).powi(2)).sum::<f64>())
                .sum::<f64>()
                / self.residues.len() as f64)
                .sqrt(),
        )
    }

    /// The same structure moved rigidly to sit on `reference` as closely as it can.
    ///
    /// Horn's quaternion method: centre both, build the `4 × 4` symmetric matrix of the
    /// correlation between them, and the eigenvector of its **largest** eigenvalue is the
    /// quaternion of the rotation that minimises the sum of squared distances. It reuses this
    /// crate's own eigensolver, which is the whole reason a superposition costs no new code.
    ///
    /// The rotation is always proper — a quaternion cannot describe a reflection — so a
    /// structure is never fitted onto its own mirror image, which is a fit that looks good and
    /// means nothing for a molecule that has a handedness.
    ///
    /// `None` if the two are different lengths.
    pub fn superposed_onto(&self, reference: &Structure) -> Option<Structure> {
        if self.residues.len() != reference.residues.len() {
            return None;
        }
        let centre = |r: &[Residue]| {
            let n = r.len() as f64;
            let mut c = [0.0; 3];
            for x in r {
                for (k, slot) in c.iter_mut().enumerate() {
                    *slot += x.at[k] / n;
                }
            }
            c
        };
        let (mine, theirs) = (centre(&self.residues), centre(&reference.residues));

        // The correlation between the two centred sets.
        let mut s = [[0.0f64; 3]; 3];
        for (a, b) in self.residues.iter().zip(&reference.residues) {
            for i in 0..3 {
                for j in 0..3 {
                    s[i][j] += (a.at[i] - mine[i]) * (b.at[j] - theirs[j]);
                }
            }
        }
        let mut n = Symmetric::zeros(4);
        n.set(0, 0, s[0][0] + s[1][1] + s[2][2]);
        n.set(1, 1, s[0][0] - s[1][1] - s[2][2]);
        n.set(2, 2, -s[0][0] + s[1][1] - s[2][2]);
        n.set(3, 3, -s[0][0] - s[1][1] + s[2][2]);
        n.set(0, 1, s[1][2] - s[2][1]);
        n.set(0, 2, s[2][0] - s[0][2]);
        n.set(0, 3, s[0][1] - s[1][0]);
        n.set(1, 2, s[0][1] + s[1][0]);
        n.set(1, 3, s[2][0] + s[0][2]);
        n.set(2, 3, s[1][2] + s[2][1]);
        let e = n.eigen();
        let q = e.vector(3);
        let (w, x, y, z) = (q[0], q[1], q[2], q[3]);
        let r = [
            [
                w * w + x * x - y * y - z * z,
                2.0 * (x * y - w * z),
                2.0 * (x * z + w * y),
            ],
            [
                2.0 * (x * y + w * z),
                w * w - x * x + y * y - z * z,
                2.0 * (y * z - w * x),
            ],
            [
                2.0 * (x * z - w * y),
                2.0 * (y * z + w * x),
                w * w - x * x - y * y + z * z,
            ],
        ];

        let residues = self
            .residues
            .iter()
            .map(|residue| {
                let v = [
                    residue.at[0] - mine[0],
                    residue.at[1] - mine[1],
                    residue.at[2] - mine[2],
                ];
                let mut at = [0.0; 3];
                for k in 0..3 {
                    at[k] = r[k][0] * v[0] + r[k][1] * v[1] + r[k][2] * v[2] + theirs[k];
                }
                Residue {
                    at,
                    ..residue.clone()
                }
            })
            .collect();
        Some(Structure { residues })
    }

    /// The unit `3n` vector pointing from this structure to `other`.
    ///
    /// What a protein *did*, ready to be compared against what a model says it can do. Superpose
    /// first: without it the vector is mostly a rigid motion, which every mode is orthogonal to
    /// by construction, so the overlap comes out small for a reason that has nothing to do with
    /// the physics.
    ///
    /// `None` if the two are different lengths, or if they are the same structure — there is no
    /// direction from a point to itself, and a zero vector normalised is `NaN`, which compares
    /// false against every threshold.
    pub fn direction_to(&self, other: &Structure) -> Option<Vec<f64>> {
        if self.residues.len() != other.residues.len() {
            return None;
        }
        let mut d = vec![0.0; 3 * self.residues.len()];
        for (i, (a, b)) in self.residues.iter().zip(&other.residues).enumerate() {
            for k in 0..3 {
                d[3 * i + k] = b.at[k] - a.at[k];
            }
        }
        let norm: f64 = d.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm == 0.0 {
            return None;
        }
        for x in d.iter_mut() {
            *x /= norm;
        }
        Some(d)
    }

    /// The mean-square displacement each residue's B-factor implies, in `m²`.
    ///
    /// `B = 8π²⟨Δr²⟩/3`, inverted. This is the **independent measurement** an elastic network
    /// model can be checked against: it is what the crystal said, not what another model
    /// computed, which is what convention 1 of this workspace asks for and what a comparison
    /// against a second implementation would not give.
    ///
    /// `None` when the file carried no B-factor column, rather than a vector of zeros that would
    /// correlate with nothing and report a correlation of `NaN`.
    pub fn experimental_fluctuations(&self) -> Option<Vec<f64>> {
        if self.residues.iter().all(|r| r.b_factor == 0.0) {
            return None;
        }
        let factor = 3.0 / (8.0 * std::f64::consts::PI * std::f64::consts::PI);
        Some(self.residues.iter().map(|r| factor * r.b_factor).collect())
    }
}
