//! A protein chain and its bound ligand, read out of a PDB entry and made into one typed system.
//!
//! An X-ray structure gives heavy atoms and nothing else: no hydrogens, no bond orders, no charges.
//! [`System::from_pdb`] supplies all three from the wwPDB Chemical Component Dictionary — one
//! [`Component`] per residue type, the twenty amino acids and the ligand, passed in by the caller —
//! and produces one [`Component`] that [`crate::uff::assign`], [`ForceField::new`], [`crate::qeq`]
//! and [`crate::solvation`] accept unchanged, with every atom marked [`Part::Protein`] or
//! [`Part::Ligand`].
//!
//! # The file format, and what is read
//!
//! **PDB format, not mmCIF**, because it is the simpler of the two to read strictly: every field of
//! an `ATOM` or `HETATM` record is at a fixed column, so a value cannot be misread by a tokeniser,
//! and the four record types read here — `HEADER` (the entry's code), `SEQRES` (the sequence),
//! `MODEL` (refused past one) and `ATOM`/`HETATM` — carry everything needed. The columns are the
//! format's (wwPDB, *Protein Data Bank Contents Guide*, version 3.30): name 13–16, alternate
//! location 17, residue 18–20, chain 22, number 23–26, insertion code 27, coordinates 31–54 and
//! element 77–78. The occupancy, B-factor and charge columns are not read: occupancy is decided
//! by the alternate-location rule below, and the charge is this module's to assign.
//!
//! # What is refused, by name
//!
//! Every one of these is a [`PdbError`] rather than a system, because each silent alternative is a
//! protein that looks complete and is not: more than one `MODEL`; an insertion code; a hydrogen in
//! the file (this module places every hydrogen, and two sources for one is ambiguous); a residue
//! with no template; an atom its template does not have, or of another element; **a residue
//! missing any heavy atom**, every missing name listed; a gap in a chain, by number or by a C–N
//! distance past [`PEPTIDE_BOND_LIMIT`]; a residue that is not the one `SEQRES` puts at its number;
//! a `HETATM` residue that is neither the named ligand nor named for dropping; and a ligand present
//! other than exactly once. The one heavy atom placed rather than read is the C-terminal `OXT`
//! (see *Termini*), and [`Placement`] says so for it.
//!
//! **Residues missing from the end of a chain** — 181L's Asn163 and Leu164, which its `REMARK 465`
//! lists as not located — are not an error: the chain is built from what was modelled, its last
//! modelled residue is the C-terminus, and [`System::unmodelled`] names what is not there. A
//! residue missing from the *middle* of a chain is a [`PdbError::ChainBreak`].
//!
//! # Hydrogens: placed by superposition, on each atom's own neighbourhood
//!
//! Each hydrogen is carried over from its template's ideal coordinates by a rigid superposition
//! ([`superpose`], Horn's closed-form quaternion, deterministic) — **not of the whole residue but
//! of the hydrogen's parent atom and its heavy neighbours** in the template onto the same atoms in
//! the crystal, with the neighbours' heavy neighbours added when that is fewer than three points.
//! Then `H = P + R (H_template − P_template)`: the hydrogen keeps its template bond length exactly
//! and its direction in the fitted frame. A whole-residue fit is reported too
//! ([`Residue::fit_rmsd`]), because it says how far the crystal's side chain is from the
//! dictionary's conformer; it is not used to place anything, because a lysine whose χ angles differ
//! from the ideal conformer's would put its ε-hydrogens an ångström from their carbon. The local
//! fragment never spans a torsion: its points are fixed by bond lengths and angles alone, so its
//! residual is the difference between the dictionary's bond geometry and the refinement's.
//!
//! **A rotor the crystal cannot fix** — a methyl, a hydroxyl, a lysine's NH₃⁺, the N-terminal
//! NH₃⁺: a heavy atom with hydrogens and one heavy neighbour, every bond of both single and not
//! aromatic — **is turned to its clearest step**: its hydrogens rotate together about the bond in
//! [`ROTOR_STEP`]s of 10° (over 120° for three hydrogens, 360° for one) from the dictionary
//! conformer's torsion, and stay at the step whose nearest heavy atom, other than the two on the
//! bond, is farthest; the first on a tie, so it repeats. **Why**: the dictionary's torsion alone
//! put Lys162's HZ3 1.495 Å from Asp159's OD1 in 181L, an overlap no structure has. The rule
//! clears contacts and knows nothing of hydrogen bonds — a hydroxyl may turn away from the acceptor
//! it donates to — and that is the minimiser's business in the next step, not a placement rule's.
//! A tyrosine's OH (on an aromatic carbon), an amide's NH₂ and arginine's NH₂ are planar by their
//! bonds and are not rotors. Rotation about the bond keeps every bond length and every angle at
//! `P`; [`System::turns`] says how far each hydrogen went. Contacts with other hydrogens are not
//! scored.
//!
//! **The backbone amide H is rebuilt from the peptide geometry**, not the template's amine:
//! on the external bisector of `C(i−1)–N–CA`, in that plane, at the template's N–H length. The
//! template's H belongs to a free amino acid's pyramidal NH₂, and a peptide nitrogen is planar.
//!
//! # Protonation at pH 7
//!
//! The dictionary's templates are free amino acids, and their own protonation is mixed: every one
//! carries a neutral NH₂ (`H`, `H2`) and a neutral COOH (`OXT`, `HXT`); **Asp and Glu are drawn
//! neutral** (`HD2`, `HE2`); **Lys, Arg and His are drawn charged** (`NZ` +1 with three H, `NH2`
//! +1, `ND1` +1 with both `HD1` and `HE2`). So at pH 7:
//!
//! | residue | template | change | charge |
//! | --- | --- | --- | --- |
//! | Asp | neutral, `HD2` | remove `HD2`, `OD2` −1 | −1 |
//! | Glu | neutral, `HE2` | remove `HE2`, `OE2` −1 | −1 |
//! | Lys | `NZ` +1 | none | +1 |
//! | Arg | `NH2` +1 | none | +1 |
//! | His in a salt bridge | `ND1` +1, `HD1` and `HE2` | none ([`Histidine::Both`]) | +1 |
//! | any other His | `ND1` +1, `HD1` and `HE2` | remove `HD1`, `ND1` 0 ([`Histidine::Epsilon`]) | 0 |
//! | Cys in a disulfide | `HG` | remove `HG` | 0 |
//! | every residue but the first | `H2` (Pro: `H`) | removed: the peptide N | — |
//! | every residue but the last | `OXT`, `HXT` | removed: the peptide C | — |
//! | N-terminus | NH₂ | `H3` added, `N` +1 | +1 |
//! | C-terminus | COOH | remove `HXT`, `OXT` −1 | −1 |
//!
//! **Histidine, by a stated rule: a histidine with either ring nitrogen within
//! [`HISTIDINE_ACCEPTOR_LIMIT`] of an Asp or Glu carboxylate oxygen is a salt bridge, doubly
//! protonated and +1 ([`Histidine::Both`]); every other histidine is neutral on Nε2**
//! ([`Histidine::Epsilon`]), the tautomer usually taken as the commoner in solution. A carboxylate
//! that close to a ring nitrogen is an ion pair, which raises the histidine's pKa — the case the
//! rule was built from is T4 lysozyme's His31–Asp70 (see [`System::from_pdb`]). So a sequence's
//! charge is `Arg + Lys − Asp − Glu` plus one per salt-bridged histidine. The rule never gives the
//! Nδ1 tautomer; [`Selection::histidine`] sets any of the three for one residue.
//!
//! # Termini
//!
//! The first modelled residue of a chain is NH₃⁺ and the last COO⁻, **whether or not the chain was
//! modelled to its end**: a chain whose last residues are disordered gets its carboxylate on the
//! last residue there is, and its `OXT` — absent from the crystal there — is placed by the same
//! superposition as a hydrogen, on `C`, `CA` and `O`. That is an artificial terminus, said so in
//! [`System::unmodelled`]; its charge is the real terminus's, so the total differs from the whole
//! sequence's only by the side chains of the unmodelled residues, which are listed. An N-terminal
//! proline is refused ([`PdbError::UnsupportedTerminus`]): its NH₂⁺ needs two hydrogens built where
//! the template has one.
//!
//! # Bonds
//!
//! Within a residue, the template's, with their orders and aromatic flags, less those to removed
//! atoms (and for a δ-protonated histidine, its Kekulé structure moved: `ND1–CE1` single,
//! `CE1=NE2`). Between residues, **one peptide bond `C(i)–N(i+1)` per consecutive pair**, single and
//! not aromatic — the amide typing rule of [`crate::uff::assign`] gives it order 1.41 — and **a
//! disulfide for every two cysteine `SG` within [`DISULFIDE_LIMIT`]**, single: a disulfide sulfur
//! has bond-order sum two, which is divalent `S_3+2` and not the hypervalent sulfur
//! [`ForceField::new`] refuses.
//!
//! # Waters, ions and alternate locations
//!
//! A `HETATM` residue is the ligand ([`Selection::new`]), or dropped because the caller named it
//! ([`Selection::dropping`]), or an error. [`System::dropped`] counts what went. **The first
//! alternate location is kept**: within a residue, the first non-blank alternate-location letter
//! met is the one read, every atom with another letter is dropped, and
//! [`System::alternates_dropped`] counts them. 181L has none.
//!
//! [`ForceField::new`]: crate::energy::ForceField::new

use crate::ccd::{Atom, Bond, BondOrder, Component, Element, ANGSTROM};
use std::collections::BTreeMap;
use std::fmt;
use std::ops::Range;

/// The twenty standard amino acids, by their dictionary codes. A protein residue (`ATOM` record)
/// named anything else is refused.
pub const AMINO_ACIDS: [&str; 20] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];

/// The longest `C(i)–N(i+1)` distance read as a peptide bond, in ångström. A peptide bond is
/// 1.33 Å; 2.0 Å is past any covalent C–N, so a pair farther apart is a break in the chain.
pub const PEPTIDE_BOND_LIMIT: f64 = 2.0;

/// The longest `SG–SG` distance read as a disulfide, in ångström. A disulfide is 2.04 Å and two
/// unbonded sulfurs touch at about 3.6 Å (twice Bondi's 1.80 Å radius), so 2.5 Å separates them.
pub const DISULFIDE_LIMIT: f64 = 2.5;

/// How near a carboxylate oxygen must be to a histidine ring nitrogen to make a salt bridge, in
/// ångström: the usual heavy-atom hydrogen-bond cutoff for N–H···O. See the module documentation.
pub const HISTIDINE_ACCEPTOR_LIMIT: f64 = 3.2;

/// Which part of the system an atom belongs to — the partition a binding energy needs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Part {
    /// A residue of a protein chain.
    Protein,
    /// The ligand.
    Ligand,
}

/// A histidine's protonation: which ring nitrogen carries the hydrogen.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Histidine {
    /// Neutral, H on Nδ1 (`HD1`; AMBER's HID).
    Delta,
    /// Neutral, H on Nε2 (`HE2`; AMBER's HIE). The rule's for a histidine with no carboxylate near.
    Epsilon,
    /// Charged, +1, H on both (AMBER's HIP) — the dictionary template as drawn, and the rule's for a salt bridge.
    Both,
}

/// How an atom's position was obtained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Placement {
    /// Read from the file.
    Crystal,
    /// Carried from the template by a superposition of `points` atoms around its parent, which
    /// fitted with residual `rmsd` (metres). See the module documentation.
    Template {
        /// The RMS distance of the fitted template atoms from the crystal's, in metres.
        rmsd: f64,
        /// How many atoms the fit used.
        points: usize,
    },
    /// A backbone amide H, on the external bisector of `C(i−1)–N–CA` at the template's N–H length.
    Backbone,
    /// The third hydrogen of an N-terminal NH₃⁺, opposite the sum of the three other bonds'
    /// directions, at the template's N–H length.
    Completion,
}

/// One residue of the system, protein or ligand.
#[derive(Clone, Debug, PartialEq)]
pub struct Residue {
    /// The dictionary code, `ASP`, `BNZ`.
    pub name: String,
    /// The chain identifier.
    pub chain: char,
    /// The residue number in the file.
    pub number: i32,
    /// Protein or ligand.
    pub part: Part,
    /// Its atoms, as a range of indices into the system's [`Component`].
    pub atoms: Range<usize>,
    /// The first modelled residue of its chain: NH₃⁺.
    pub n_terminal: bool,
    /// The last modelled residue of its chain: COO⁻.
    pub c_terminal: bool,
    /// For a histidine, the protonation it was given.
    pub histidine: Option<Histidine>,
    /// Whether its `SG` is in a disulfide.
    pub disulfide: bool,
    /// The sum of its atoms' formal charges.
    pub charge: i32,
    /// The RMS distance, in metres, after superposing the template's ideal heavy atoms on the
    /// crystal's — the whole residue at once. A measure of the conformer's difference from the
    /// dictionary's, reported and not used to place anything.
    pub fit_rmsd: f64,
}

impl Residue {
    /// `A:ASP70`, as errors and atom names write it.
    pub fn label(&self) -> String {
        format!("{}:{}{}", self.chain, self.name, self.number)
    }
}

/// A residue `SEQRES` lists and the file has no coordinates for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Unmodelled {
    /// The chain.
    pub chain: char,
    /// Its number: its position in `SEQRES`, counting from one.
    pub number: i32,
    /// Its dictionary code.
    pub name: String,
}

/// A `HETATM` residue type the selection dropped, and how much of it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Dropped {
    /// The dictionary code, `HOH`.
    pub name: String,
    /// How many residues of it.
    pub residues: usize,
    /// How many atom records.
    pub atoms: usize,
}

/// What to keep from a PDB entry besides its protein chains.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Selection {
    ligand: String,
    drop: Vec<String>,
    histidines: Vec<(char, i32, Histidine)>,
}

impl Selection {
    /// Keep the ligand named `ligand` (its dictionary code, `BNZ`), which must occur exactly once,
    /// and drop no other `HETATM` residue — any there are is an error until named in
    /// [`Selection::dropping`].
    pub fn new(ligand: impl Into<String>) -> Selection {
        Selection {
            ligand: ligand.into(),
            drop: Vec::new(),
            histidines: Vec::new(),
        }
    }

    /// The same selection dropping every `HETATM` residue named in `names` — waters, ions, buffer
    /// components.
    pub fn dropping(mut self, names: &[&str]) -> Selection {
        self.drop.extend(names.iter().map(|s| s.to_string()));
        self
    }

    /// The same selection with histidine `number` of `chain` given `state`, instead of the rule's.
    pub fn histidine(mut self, chain: char, number: i32, state: Histidine) -> Selection {
        self.histidines.push((chain, number, state));
        self
    }
}

/// Why a PDB entry did not produce a system. Every variant names what it refused.
#[derive(Clone, Debug, PartialEq)]
pub enum PdbError {
    /// An `ATOM`, `HETATM` or `SEQRES` record whose fixed-column field could not be read.
    Record {
        /// Which line, counting from one.
        line: usize,
        /// The field.
        field: &'static str,
        /// What it said.
        text: String,
    },
    /// More than one `MODEL`: an ensemble, and taking the first would be reading one of several.
    Models {
        /// How many.
        count: usize,
    },
    /// A residue with an insertion code, which this reader does not number.
    InsertionCode {
        /// The residue, `A:GLY52A`.
        residue: String,
    },
    /// A hydrogen in the file. This module places every hydrogen; two sources for one is
    /// ambiguous.
    HydrogenInFile {
        /// The residue.
        residue: String,
        /// The atom.
        atom: String,
    },
    /// Two template components with one code.
    DuplicateTemplate {
        /// The code.
        template: String,
    },
    /// A residue with no template among those given, or a protein residue that is not one of the
    /// twenty amino acids.
    NoTemplate {
        /// The residue.
        residue: String,
    },
    /// A `HETATM` residue that is neither the ligand nor named for dropping.
    Unexpected {
        /// The residue.
        residue: String,
    },
    /// The ligand occurs other than exactly once.
    LigandCount {
        /// Its code.
        ligand: String,
        /// How many times it occurs.
        count: usize,
    },
    /// One residue number used for two residue names.
    MixedResidue {
        /// The chain and number, `A:70`.
        residue: String,
        /// The two names.
        names: [String; 2],
    },
    /// An atom whose name the residue's template does not have (after protonation: an `OXT` on a
    /// residue that is not a C-terminus is one).
    UnknownAtom {
        /// The residue.
        residue: String,
        /// The atom.
        atom: String,
    },
    /// An atom whose element is not its template's.
    ElementMismatch {
        /// The residue.
        residue: String,
        /// The atom.
        atom: String,
        /// The element column as the file gives it.
        file: String,
        /// The template's element.
        template: Element,
    },
    /// One atom name twice in a residue, after the alternate-location rule.
    DuplicateAtom {
        /// The residue.
        residue: String,
        /// The atom.
        atom: String,
    },
    /// A residue without every heavy atom its template has.
    MissingAtoms {
        /// The residue.
        residue: String,
        /// Every missing atom, in the template's order.
        atoms: Vec<String>,
    },
    /// A protein chain with no `SEQRES`.
    NoSequence {
        /// The chain.
        chain: char,
    },
    /// A residue that is not the one `SEQRES` puts at its number, or whose number is outside it.
    SequenceMismatch {
        /// The residue.
        residue: String,
        /// What `SEQRES` has there, if anything.
        seqres: Option<String>,
    },
    /// Two consecutive modelled residues of one chain that are not bonded: numbers not
    /// consecutive, or `C–N` past [`PEPTIDE_BOND_LIMIT`].
    ChainBreak {
        /// The residue before.
        after: String,
        /// The residue after.
        before: String,
        /// The `C–N` distance, ångström.
        distance: f64,
    },
    /// A terminus this module does not build: an N-terminal proline.
    UnsupportedTerminus {
        /// The residue.
        residue: String,
    },
    /// A template without an atom or bond the protonation rules name — not the dictionary's entry
    /// for that code.
    TemplateMismatch {
        /// The template's code.
        template: String,
        /// What is missing.
        what: String,
    },
    /// An atom whose parent's neighbourhood has fewer than three atoms with positions, so no
    /// superposition can place it.
    Unplaceable {
        /// The residue.
        residue: String,
        /// The atom.
        atom: String,
    },
    /// A cysteine `SG` within [`DISULFIDE_LIMIT`] of more than one other.
    SulfurPartners {
        /// The residue.
        residue: String,
    },
    /// A [`Selection::histidine`] naming a residue that is not a histidine.
    NoSuchHistidine {
        /// The chain and number, `A:31`.
        residue: String,
    },
    /// No protein residue at all.
    NoProtein,
}

impl fmt::Display for PdbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PdbError::Record { line, field, text } => {
                write!(
                    f,
                    "line {line}: the {field} field is {text:?}, which it cannot be"
                )
            }
            PdbError::Models { count } => write!(
                f,
                "{count} models; an ensemble is not one structure, and taking the first would \
                 read one of several"
            ),
            PdbError::InsertionCode { residue } => {
                write!(
                    f,
                    "{residue} has an insertion code, which this reader does not number"
                )
            }
            PdbError::HydrogenInFile { residue, atom } => write!(
                f,
                "{residue} {atom} is a hydrogen in the file; every hydrogen is placed here, and \
                 two sources for one is ambiguous"
            ),
            PdbError::DuplicateTemplate { template } => {
                write!(f, "two templates are named {template}")
            }
            PdbError::NoTemplate { residue } => write!(f, "{residue} has no template"),
            PdbError::Unexpected { residue } => write!(
                f,
                "{residue} is a HETATM residue that is neither the ligand nor named for dropping"
            ),
            PdbError::LigandCount { ligand, count } => write!(
                f,
                "the ligand {ligand} occurs {count} times; it must occur exactly once"
            ),
            PdbError::MixedResidue { residue, names } => {
                write!(f, "{residue} is both {} and {}", names[0], names[1])
            }
            PdbError::UnknownAtom { residue, atom } => {
                write!(f, "{residue} has an atom {atom} that its template does not")
            }
            PdbError::ElementMismatch {
                residue,
                atom,
                file,
                template,
            } => write!(
                f,
                "{residue} {atom} is element {file:?} in the file and {template} in its template"
            ),
            PdbError::DuplicateAtom { residue, atom } => {
                write!(f, "{residue} has two atoms named {atom}")
            }
            PdbError::MissingAtoms { residue, atoms } => {
                write!(f, "{residue} is missing heavy atoms {}", atoms.join(", "))
            }
            PdbError::NoSequence { chain } => write!(f, "chain {chain} has no SEQRES"),
            PdbError::SequenceMismatch { residue, seqres } => match seqres {
                Some(s) => write!(f, "{residue} is {s} in SEQRES"),
                None => write!(f, "{residue} is outside SEQRES"),
            },
            PdbError::ChainBreak {
                after,
                before,
                distance,
            } => write!(
                f,
                "the chain breaks between {after} and {before} (C-N {distance:.2} Å)"
            ),
            PdbError::UnsupportedTerminus { residue } => write!(
                f,
                "{residue} is an N-terminal proline, whose NH2+ this module does not build"
            ),
            PdbError::TemplateMismatch { template, what } => write!(
                f,
                "the template {template} has no {what}, which the protonation rules need"
            ),
            PdbError::Unplaceable { residue, atom } => write!(
                f,
                "{residue} {atom} cannot be placed: its parent has fewer than three positioned \
                 neighbours"
            ),
            PdbError::SulfurPartners { residue } => write!(
                f,
                "{residue} SG is within {DISULFIDE_LIMIT} Å of more than one other SG"
            ),
            PdbError::NoSuchHistidine { residue } => {
                write!(f, "{residue} is named as a histidine and is not one")
            }
            PdbError::NoProtein => f.write_str("the file has no protein residue"),
        }
    }
}

impl std::error::Error for PdbError {}

/// A rigid superposition: `to + R (p − from)` takes a point of the moving set onto the target.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fit {
    /// The rotation, row-major.
    pub rotation: [[f64; 3]; 3],
    /// The moving set's centroid.
    pub from: [f64; 3],
    /// The target set's centroid.
    pub to: [f64; 3],
    /// The RMS distance between the moved set and the target, in their unit.
    pub rmsd: f64,
}

impl Fit {
    /// `R v`: a direction moved.
    pub fn rotate(&self, v: [f64; 3]) -> [f64; 3] {
        let r = &self.rotation;
        [
            r[0][0] * v[0] + r[0][1] * v[1] + r[0][2] * v[2],
            r[1][0] * v[0] + r[1][1] * v[1] + r[1][2] * v[2],
            r[2][0] * v[0] + r[2][1] * v[1] + r[2][2] * v[2],
        ]
    }

    /// `to + R (p − from)`: a point moved.
    pub fn apply(&self, p: [f64; 3]) -> [f64; 3] {
        add(self.to, self.rotate(sub(p, self.from)))
    }
}

/// The proper rotation and translation that best superpose `moving` on `target` in the
/// least-squares sense, by Horn's closed form: the rotation is the unit quaternion that is the
/// eigenvector of the largest eigenvalue of a 4×4 symmetric matrix built from the two sets'
/// cross-covariance (B. K. P. Horn, *J. Opt. Soc. Am. A* **4**, 629 (1987)). The eigenvector is
/// found by cyclic Jacobi rotations in a fixed order, so the result repeats bit for bit. A
/// quaternion is always a proper rotation: a mirror image is not superposed by reflecting it.
///
/// # Panics
///
/// If the two sets differ in length or are empty.
pub fn superpose(moving: &[[f64; 3]], target: &[[f64; 3]]) -> Fit {
    assert_eq!(moving.len(), target.len(), "one target per moving point");
    assert!(!moving.is_empty(), "nothing to superpose");
    let n = moving.len() as f64;
    let centroid = |s: &[[f64; 3]]| {
        let mut c = [0.0; 3];
        for p in s {
            for k in 0..3 {
                c[k] += p[k];
            }
        }
        c.map(|x| x / n)
    };
    let (from, to) = (centroid(moving), centroid(target));
    let mut s = [[0.0f64; 3]; 3];
    for (a, b) in moving.iter().zip(target) {
        let (a, b) = (sub(*a, from), sub(*b, to));
        for i in 0..3 {
            for j in 0..3 {
                s[i][j] += a[i] * b[j];
            }
        }
    }
    let [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = s;
    let m = [
        [xx + yy + zz, yz - zy, zx - xz, xy - yx],
        [yz - zy, xx - yy - zz, xy + yx, zx + xz],
        [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
        [xy - yx, zx + xz, yz + zy, -xx - yy + zz],
    ];
    let q = largest_eigenvector(m);
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    let [a, b, c, d] = q.map(|x| x / norm);
    let rotation = [
        [
            a * a + b * b - c * c - d * d,
            2.0 * (b * c - a * d),
            2.0 * (b * d + a * c),
        ],
        [
            2.0 * (b * c + a * d),
            a * a - b * b + c * c - d * d,
            2.0 * (c * d - a * b),
        ],
        [
            2.0 * (b * d - a * c),
            2.0 * (c * d + a * b),
            a * a - b * b - c * c + d * d,
        ],
    ];
    let mut fit = Fit {
        rotation,
        from,
        to,
        rmsd: 0.0,
    };
    let sq: f64 = moving
        .iter()
        .zip(target)
        .map(|(a, b)| norm2(sub(fit.apply(*a), *b)))
        .sum();
    fit.rmsd = (sq / n).sqrt();
    fit
}

/// The eigenvector of the largest eigenvalue of a symmetric 4×4 matrix, by cyclic Jacobi.
#[allow(clippy::needless_range_loop)]
fn largest_eigenvector(mut a: [[f64; 4]; 4]) -> [f64; 4] {
    let mut v = [[0.0f64; 4]; 4];
    for (k, row) in v.iter_mut().enumerate() {
        row[k] = 1.0;
    }
    let scale: f64 = a.iter().flatten().map(|x| x * x).sum();
    for _sweep in 0..64 {
        let mut off = 0.0;
        for p in 0..4 {
            for q in p + 1..4 {
                off += a[p][q] * a[p][q];
            }
        }
        if off <= 1e-32 * scale {
            break;
        }
        for p in 0..4 {
            for q in p + 1..4 {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = if theta >= 0.0 {
                    1.0 / (theta + (theta * theta + 1.0).sqrt())
                } else {
                    -1.0 / (-theta + (theta * theta + 1.0).sqrt())
                };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..4 {
                    let (kp, kq) = (a[k][p], a[k][q]);
                    a[k][p] = c * kp - s * kq;
                    a[k][q] = s * kp + c * kq;
                }
                for k in 0..4 {
                    let (pk, qk) = (a[p][k], a[q][k]);
                    a[p][k] = c * pk - s * qk;
                    a[q][k] = s * pk + c * qk;
                }
                for k in 0..4 {
                    let (kp, kq) = (v[k][p], v[k][q]);
                    v[k][p] = c * kp - s * kq;
                    v[k][q] = s * kp + c * kq;
                }
            }
        }
    }
    let mut best = 0;
    for k in 1..4 {
        if a[k][k] > a[best][best] {
            best = k;
        }
    }
    [v[0][best], v[1][best], v[2][best], v[3][best]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn norm2(a: [f64; 3]) -> f64 {
    a[0] * a[0] + a[1] * a[1] + a[2] * a[2]
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    scale(a, 1.0 / norm2(a).sqrt())
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm2(sub(a, b)).sqrt()
}

/// One `ATOM` or `HETATM` record, positions in ångström.
#[derive(Clone, Debug)]
struct Record {
    het: bool,
    name: String,
    alt: char,
    residue: String,
    element: String,
    at: [f64; 3],
}

/// One residue's records, in file order.
#[derive(Clone, Debug)]
struct Raw {
    name: String,
    chain: char,
    number: i32,
    het: bool,
    records: Vec<Record>,
    alternates_dropped: usize,
}

impl Raw {
    fn label(&self) -> String {
        format!("{}:{}{}", self.chain, self.name, self.number)
    }

    fn position(&self, atom: &str) -> Option<[f64; 3]> {
        self.records.iter().find(|r| r.name == atom).map(|r| r.at)
    }
}

/// A fixed-column field: columns `from` to `to`, counting from one as the format does.
fn columns(line: &str, from: usize, to: usize) -> &str {
    let end = to.min(line.len());
    line.get(from - 1..end).unwrap_or("")
}

fn record(line: &str, number: usize) -> Result<Record, PdbError> {
    let bad = |field: &'static str, text: &str| PdbError::Record {
        line: number,
        field,
        text: text.into(),
    };
    if !line.is_ascii() {
        return Err(bad("record", "non-ASCII text"));
    }
    if line.len() < 54 {
        return Err(bad("record", line));
    }
    let coordinate = |from: usize, field: &'static str| {
        let text = columns(line, from, from + 7).trim();
        text.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| bad(field, text))
    };
    let name = columns(line, 13, 16).trim().to_string();
    if name.is_empty() {
        return Err(bad("atom name", ""));
    }
    Ok(Record {
        het: line.starts_with("HETATM"),
        name,
        alt: columns(line, 17, 17).chars().next().unwrap_or(' '),
        residue: columns(line, 18, 20).trim().to_string(),
        element: columns(line, 77, 78).trim().to_string(),
        at: [
            coordinate(31, "x")?,
            coordinate(39, "y")?,
            coordinate(47, "z")?,
        ],
    })
}

/// Everything the file says that is read: the code, the sequences, and the residues.
struct Entry {
    code: Option<String>,
    sequences: BTreeMap<char, Vec<String>>,
    residues: Vec<Raw>,
}

fn read_entry(text: &str) -> Result<Entry, PdbError> {
    let mut code = None;
    let mut sequences: BTreeMap<char, Vec<String>> = BTreeMap::new();
    let mut residues: Vec<Raw> = Vec::new();
    let mut index: BTreeMap<(char, i32), usize> = BTreeMap::new();
    let mut models = 0;
    for (k, line) in text.lines().enumerate() {
        let number = k + 1;
        if line.starts_with("HEADER") {
            let c = columns(line, 63, 66).trim();
            if !c.is_empty() {
                code = Some(c.to_string());
            }
        } else if line.starts_with("MODEL ") {
            models += 1;
        } else if line.starts_with("SEQRES") {
            let chain = columns(line, 12, 12).chars().next().unwrap_or(' ');
            let names = columns(line, 20, 80).split_whitespace().map(str::to_string);
            sequences.entry(chain).or_default().extend(names);
        } else if line.starts_with("ATOM  ") || line.starts_with("HETATM") {
            let r = record(line, number)?;
            let chain = columns(line, 22, 22).chars().next().unwrap_or(' ');
            let n_text = columns(line, 23, 26).trim();
            let n: i32 = n_text.parse().map_err(|_| PdbError::Record {
                line: number,
                field: "residue number",
                text: n_text.into(),
            })?;
            let label = || format!("{chain}:{}{n}", r.residue);
            let insertion = columns(line, 27, 27).trim();
            if !insertion.is_empty() {
                return Err(PdbError::InsertionCode {
                    residue: format!("{}{insertion}", label()),
                });
            }
            let at = match index.get(&(chain, n)) {
                Some(&at) => at,
                None => {
                    index.insert((chain, n), residues.len());
                    residues.push(Raw {
                        name: r.residue.clone(),
                        chain,
                        number: n,
                        het: r.het,
                        records: Vec::new(),
                        alternates_dropped: 0,
                    });
                    residues.len() - 1
                }
            };
            let raw = &mut residues[at];
            if raw.name != r.residue {
                return Err(PdbError::MixedResidue {
                    residue: format!("{chain}:{n}"),
                    names: [raw.name.clone(), r.residue.clone()],
                });
            }
            // The first alternate location met in a residue is the one kept.
            if r.alt != ' ' {
                let first = raw.records.iter().map(|x| x.alt).find(|&a| a != ' ');
                if first.is_some_and(|f| f != r.alt) {
                    raw.alternates_dropped += 1;
                    continue;
                }
            }
            if raw.records.iter().any(|x| x.name == r.name) {
                return Err(PdbError::DuplicateAtom {
                    residue: label(),
                    atom: r.name,
                });
            }
            raw.records.push(r);
        }
    }
    if models > 1 {
        return Err(PdbError::Models { count: models });
    }
    Ok(Entry {
        code,
        sequences,
        residues,
    })
}

/// A residue's role in its chain and the choices made for it.
#[derive(Clone, Copy, Debug)]
struct Role {
    part: Part,
    n_terminal: bool,
    c_terminal: bool,
    histidine: Option<Histidine>,
    disulfide: bool,
}

/// The system as it is assembled.
#[derive(Default)]
struct Builder {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    local: Vec<String>,
    residue_of: Vec<usize>,
    placements: Vec<Placement>,
    turns: Vec<f64>,
    residues: Vec<Residue>,
}

fn template_index(t: &Component, name: &str) -> Option<usize> {
    t.atoms().iter().position(|a| a.name == name)
}

fn need(t: &Component, name: &str) -> Result<usize, PdbError> {
    template_index(t, name).ok_or_else(|| PdbError::TemplateMismatch {
        template: t.id().into(),
        what: format!("atom {name}"),
    })
}

fn angstrom(p: [f64; 3]) -> [f64; 3] {
    p.map(|x| x / ANGSTROM)
}

/// Template atoms bonded to `i` that are heavy, kept, and positioned, excluding `not`.
fn heavy_neighbours(
    t: &Component,
    i: usize,
    keep: &[bool],
    at: &[Option<[f64; 3]>],
    not: usize,
) -> Vec<usize> {
    t.neighbours(i)
        .map(|(j, _)| j)
        .filter(|&j| j != not && keep[j] && t.atoms()[j].element != Element::H && at[j].is_some())
        .collect()
}

/// Place template atom `i` by superposing its parent's neighbourhood: see the module
/// documentation. Returns the position (Å) and the fit's residual (Å) and size.
fn place_on_parent(
    t: &Component,
    i: usize,
    keep: &[bool],
    at: &[Option<[f64; 3]>],
    label: &str,
) -> Result<([f64; 3], f64, usize), PdbError> {
    let unplaceable = || PdbError::Unplaceable {
        residue: label.into(),
        atom: t.atoms()[i].name.clone(),
    };
    let parent = t
        .neighbours(i)
        .map(|(j, _)| j)
        .find(|&j| keep[j] && t.atoms()[j].element != Element::H && at[j].is_some())
        .ok_or_else(unplaceable)?;
    let mut fragment = vec![parent];
    fragment.extend(heavy_neighbours(t, parent, keep, at, i));
    if fragment.len() < 3 {
        let shell: Vec<usize> = fragment[1..].to_vec();
        for j in shell {
            for k in heavy_neighbours(t, j, keep, at, i) {
                if !fragment.contains(&k) {
                    fragment.push(k);
                }
            }
        }
    }
    if fragment.len() < 3 {
        return Err(unplaceable());
    }
    let moving: Vec<[f64; 3]> = fragment
        .iter()
        .map(|&k| angstrom(t.atoms()[k].at))
        .collect();
    let target: Vec<[f64; 3]> = fragment
        .iter()
        .map(|&k| at[k].expect("filtered to positioned atoms"))
        .collect();
    let fit = superpose(&moving, &target);
    let bond = sub(angstrom(t.atoms()[i].at), angstrom(t.atoms()[parent].at));
    let p = at[parent].expect("filtered to a positioned parent");
    Ok((add(p, fit.rotate(bond)), fit.rmsd, fragment.len()))
}

/// The removals and charge changes the protonation rules make to template `t` in role `role`.
#[allow(clippy::type_complexity)]
fn protonation(t: &Component, role: Role) -> (Vec<&'static str>, Vec<(&'static str, i32)>) {
    let mut remove: Vec<&'static str> = Vec::new();
    let mut charge: Vec<(&'static str, i32)> = Vec::new();
    if role.part == Part::Ligand {
        return (remove, charge);
    }
    if role.n_terminal {
        charge.push(("N", 1));
    } else if t.id() == "PRO" {
        remove.push("H");
    } else {
        remove.push("H2");
    }
    if role.c_terminal {
        remove.push("HXT");
        charge.push(("OXT", -1));
    } else {
        remove.extend(["OXT", "HXT"]);
    }
    match t.id() {
        "ASP" => {
            remove.push("HD2");
            charge.push(("OD2", -1));
        }
        "GLU" => {
            remove.push("HE2");
            charge.push(("OE2", -1));
        }
        "HIS" => match role.histidine.unwrap_or(Histidine::Epsilon) {
            Histidine::Delta => {
                remove.push("HE2");
                charge.push(("ND1", 0));
            }
            Histidine::Epsilon => {
                remove.push("HD1");
                charge.push(("ND1", 0));
            }
            Histidine::Both => charge.push(("ND1", 1)),
        },
        "CYS" if role.disulfide => remove.push("HG"),
        _ => {}
    }
    (remove, charge)
}

/// Build one residue from its template into `out`.
fn build(
    raw: &Raw,
    t: &Component,
    role: Role,
    previous_c: Option<[f64; 3]>,
    out: &mut Builder,
) -> Result<(), PdbError> {
    let label = raw.label();
    if role.n_terminal && t.id() == "PRO" {
        return Err(PdbError::UnsupportedTerminus { residue: label });
    }
    let (remove, charges) = protonation(t, role);
    let n = t.atoms().len();
    let mut keep = vec![true; n];
    for name in &remove {
        keep[need(t, name)?] = false;
    }
    for (name, _) in &charges {
        need(t, name)?;
    }
    // Crystal atoms onto template atoms.
    let mut at: Vec<Option<[f64; 3]>> = vec![None; n];
    let mut placement = vec![Placement::Crystal; n];
    for r in &raw.records {
        if r.element.eq_ignore_ascii_case("H") || r.element.eq_ignore_ascii_case("D") {
            return Err(PdbError::HydrogenInFile {
                residue: label,
                atom: r.name.clone(),
            });
        }
        let i = template_index(t, &r.name)
            .filter(|&i| keep[i])
            .ok_or_else(|| PdbError::UnknownAtom {
                residue: label.clone(),
                atom: r.name.clone(),
            })?;
        let element = t.atoms()[i].element;
        if Element::from_symbol(&r.element) != Some(element) {
            return Err(PdbError::ElementMismatch {
                residue: label,
                atom: r.name.clone(),
                file: r.element.clone(),
                template: element,
            });
        }
        at[i] = Some(r.at);
    }
    let oxt = template_index(t, "OXT");
    let missing: Vec<String> = (0..n)
        .filter(|&i| keep[i] && t.atoms()[i].element != Element::H && at[i].is_none())
        .filter(|&i| !(role.c_terminal && Some(i) == oxt))
        .map(|i| t.atoms()[i].name.clone())
        .collect();
    if !missing.is_empty() {
        return Err(PdbError::MissingAtoms {
            residue: label,
            atoms: missing,
        });
    }
    // The whole-residue fit, reported.
    let fitted: Vec<usize> = (0..n).filter(|&i| at[i].is_some()).collect();
    let fit_rmsd = if fitted.len() >= 3 {
        let moving: Vec<[f64; 3]> = fitted.iter().map(|&i| angstrom(t.atoms()[i].at)).collect();
        let target: Vec<[f64; 3]> = fitted.iter().map(|&i| at[i].expect("fitted")).collect();
        superpose(&moving, &target).rmsd * ANGSTROM
    } else {
        0.0
    };
    // A C-terminal OXT the crystal lacks, then every hydrogen.
    if let Some(o) = oxt.filter(|&o| keep[o] && at[o].is_none()) {
        let (p, rmsd, points) = place_on_parent(t, o, &keep, &at, &label)?;
        at[o] = Some(p);
        placement[o] = Placement::Template {
            rmsd: rmsd * ANGSTROM,
            points,
        };
    }
    let template_nh = |t: &Component| -> Result<f64, PdbError> {
        let (nn, h) = (need(t, "N")?, need(t, "H")?);
        Ok(distance(
            angstrom(t.atoms()[nn].at),
            angstrom(t.atoms()[h].at),
        ))
    };
    for i in 0..n {
        if !keep[i] || t.atoms()[i].element != Element::H {
            continue;
        }
        let backbone = role.part == Part::Protein && !role.n_terminal && t.atoms()[i].name == "H";
        if backbone {
            let c = previous_c.expect("every residue but a chain's first has a previous C");
            let nn = at[need(t, "N")?].expect("N is a checked heavy atom");
            let ca = at[need(t, "CA")?].expect("CA is a checked heavy atom");
            let d = unit(scale(add(unit(sub(c, nn)), unit(sub(ca, nn))), -1.0));
            at[i] = Some(add(nn, scale(d, template_nh(t)?)));
            placement[i] = Placement::Backbone;
        } else {
            let (p, rmsd, points) = place_on_parent(t, i, &keep, &at, &label)?;
            at[i] = Some(p);
            placement[i] = Placement::Template {
                rmsd: rmsd * ANGSTROM,
                points,
            };
        }
    }
    // Append the residue.
    let start = out.atoms.len();
    let residue_index = out.residues.len();
    let mut local = vec![usize::MAX; n];
    for i in 0..n {
        if !keep[i] {
            continue;
        }
        let a = &t.atoms()[i];
        let charge = charges
            .iter()
            .find(|(name, _)| *name == a.name)
            .map_or(a.charge, |&(_, q)| q);
        local[i] = out.atoms.len();
        out.atoms.push(Atom {
            name: format!("{label}:{}", a.name),
            element: a.element,
            charge,
            aromatic: a.aromatic,
            at: at[i].expect("every kept atom placed").map(|x| x * ANGSTROM),
        });
        out.local.push(a.name.clone());
        out.residue_of.push(residue_index);
        out.placements.push(placement[i]);
        out.turns.push(0.0);
    }
    let delta = t.id() == "HIS" && role.histidine == Some(Histidine::Delta);
    let mut moved = 0;
    for b in t.bonds() {
        let [i, j] = b.atoms;
        if !keep[i] || !keep[j] {
            continue;
        }
        let mut order = b.order;
        if delta {
            let names = [t.atoms()[i].name.as_str(), t.atoms()[j].name.as_str()];
            let is = |x: &str, y: &str| names == [x, y] || names == [y, x];
            if is("ND1", "CE1") {
                order = BondOrder::Single;
                moved += 1;
            } else if is("CE1", "NE2") {
                order = BondOrder::Double;
                moved += 1;
            }
        }
        out.bonds.push(Bond {
            atoms: [local[i], local[j]],
            order,
            aromatic: b.aromatic,
        });
    }
    if delta && moved != 2 {
        return Err(PdbError::TemplateMismatch {
            template: t.id().into(),
            what: "the ND1-CE1 and CE1-NE2 bonds".into(),
        });
    }
    if role.part == Part::Protein && role.n_terminal {
        if template_index(t, "H3").is_some() {
            return Err(PdbError::TemplateMismatch {
                template: t.id().into(),
                what: "free name H3 for the third N-terminal hydrogen".into(),
            });
        }
        let nn = local[need(t, "N")?];
        let n_at = angstrom(out.atoms[nn].at);
        let around: Vec<[f64; 3]> = out
            .bonds
            .iter()
            .filter_map(|b| match b.atoms {
                [x, o] | [o, x] if x == nn => Some(o),
                _ => None,
            })
            .map(|o| unit(sub(angstrom(out.atoms[o].at), n_at)))
            .collect();
        let sum = around.iter().fold([0.0; 3], |s, &u| add(s, u));
        let h = add(n_at, scale(unit(scale(sum, -1.0)), template_nh(t)?));
        let index = out.atoms.len();
        out.atoms.push(Atom {
            name: format!("{label}:H3"),
            element: Element::H,
            charge: 0,
            aromatic: false,
            at: h.map(|x| x * ANGSTROM),
        });
        out.local.push("H3".into());
        out.residue_of.push(residue_index);
        out.placements.push(Placement::Completion);
        out.turns.push(0.0);
        out.bonds.push(Bond {
            atoms: [nn, index],
            order: BondOrder::Single,
            aromatic: false,
        });
    }
    let end = out.atoms.len();
    out.residues.push(Residue {
        name: raw.name.clone(),
        chain: raw.chain,
        number: raw.number,
        part: role.part,
        atoms: start..end,
        n_terminal: role.n_terminal,
        c_terminal: role.c_terminal,
        histidine: role.histidine,
        disulfide: role.disulfide,
        charge: out.atoms[start..end].iter().map(|a| a.charge).sum(),
        fit_rmsd,
    });
    Ok(())
}

/// The step a rotor is turned in, radians: 10°.
pub const ROTOR_STEP: f64 = std::f64::consts::PI / 18.0;

/// Turn every rotor to its clearest step: see the module documentation.
///
/// A rotor is a heavy atom `P` with hydrogens and exactly one heavy neighbour `A`, every bond of
/// both single and none aromatic — a methyl, a hydroxyl, a thiol, an NH₃⁺. Its hydrogens are
/// turned together about `A → P` in [`ROTOR_STEP`]s, over 120° for three hydrogens and 360°
/// otherwise, from the dictionary conformer's torsion, and kept at the step whose nearest heavy
/// atom other than `P` and `A` is farthest; the first such step on a tie, so the result repeats.
fn turn_rotors(out: &mut Builder) {
    let n = out.atoms.len();
    let mut adjacent: Vec<Vec<usize>> = vec![Vec::new(); n];
    let mut plain = vec![true; n];
    for b in &out.bonds {
        let [i, j] = b.atoms;
        adjacent[i].push(j);
        adjacent[j].push(i);
        if b.order != BondOrder::Single || b.aromatic {
            plain[i] = false;
            plain[j] = false;
        }
    }
    let hydrogen: Vec<bool> = out.atoms.iter().map(|a| a.element == Element::H).collect();
    let is_h = |i: usize| hydrogen[i];
    let heavy: Vec<usize> = (0..n).filter(|&i| !is_h(i)).collect();
    let at: Vec<[f64; 3]> = out.atoms.iter().map(|a| angstrom(a.at)).collect();
    let at = &at;
    let heavy = &heavy;
    for &p in heavy {
        let heavy_nb: Vec<usize> = adjacent[p].iter().copied().filter(|&j| !is_h(j)).collect();
        let hydrogens: Vec<usize> = adjacent[p].iter().copied().filter(|&j| is_h(j)).collect();
        let [a] = heavy_nb[..] else { continue };
        if hydrogens.is_empty() || !plain[p] || !plain[a] {
            continue;
        }
        let axis = unit(sub(at[p], at[a]));
        let span = if hydrogens.len() == 3 {
            2.0 * std::f64::consts::PI / 3.0
        } else {
            2.0 * std::f64::consts::PI
        };
        let steps = (span / ROTOR_STEP).round() as usize;
        let turned = |k: usize| -> Vec<[f64; 3]> {
            let phi = k as f64 * ROTOR_STEP;
            hydrogens
                .iter()
                .map(|&h| add(at[p], rotate_about(sub(at[h], at[p]), axis, phi)))
                .collect()
        };
        let clearance = |hs: &[[f64; 3]]| {
            hs.iter()
                .flat_map(|&h| {
                    heavy
                        .iter()
                        .filter(|&&x| x != p && x != a)
                        .map(move |&x| norm2(sub(h, at[x])))
                })
                .fold(f64::INFINITY, f64::min)
        };
        let (mut best, mut best_clearance) = (0, clearance(&turned(0)));
        for k in 1..steps {
            let c = clearance(&turned(k));
            if c > best_clearance {
                best = k;
                best_clearance = c;
            }
        }
        if best == 0 {
            continue;
        }
        let turn = best as f64 * ROTOR_STEP;
        for (&h, p_new) in hydrogens.iter().zip(turned(best)) {
            out.atoms[h].at = p_new.map(|x| x * ANGSTROM);
            out.turns[h] = turn;
        }
    }
}

/// `v` turned by `phi` about the unit vector `u` (Rodrigues).
fn rotate_about(v: [f64; 3], u: [f64; 3], phi: f64) -> [f64; 3] {
    let (s, c) = phi.sin_cos();
    let along = scale(u, (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) * (1.0 - c));
    let cross = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    add(add(scale(v, c), scale(cross, s)), along)
}

/// A protein and its ligand as one typed system: see the module documentation.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    component: Component,
    residues: Vec<Residue>,
    residue_of: Vec<usize>,
    local: Vec<String>,
    placements: Vec<Placement>,
    turns: Vec<f64>,
    peptide_bonds: Vec<usize>,
    disulfides: Vec<usize>,
    unmodelled: Vec<Unmodelled>,
    dropped: Vec<Dropped>,
    alternates_dropped: usize,
}

impl System {
    /// Read the protein chains of the PDB-format entry `text`, with the one ligand `selection`
    /// names, give them hydrogens, bonds and formal charges from `templates` — the dictionary
    /// entries of every amino acid present and of the ligand — and assemble one [`Component`].
    ///
    /// Protein residues come first, chain by chain in the file's order, then the ligand; within a
    /// residue the atoms are in the template's order, less those protonation removes, with an
    /// N-terminal `H3` last. Every atom is named `chain:RESnumber:atom`, `A:ASP70:OD2`.
    ///
    /// **T4 lysozyme has one histidine, His31**, whose Nδ1 is 2.66 Å from Asp70's Oδ2 in 181L, so
    /// the rule makes it a salt bridge, doubly protonated, +1 — and T4 lysozyme's total +9, the
    /// sequence's `Arg + Lys − Asp − Glu` of +8 and His31. The pair is the subject of
    /// Anderson, Becktel and Dahlquist, "pH-Induced Denaturation of Proteins: A Single
    /// Salt Bridge Contributes 3–5 kcal/mol to the Free Energy of Folding of T4 Lysozyme",
    /// *Biochemistry* **29**, 2403 (1990) — its title and citation verified, its body not read.
    /// [`Selection::histidine`] gives either neutral state instead.
    ///
    /// # Errors
    ///
    /// Any [`PdbError`], each naming what it refused; the module documentation lists them.
    pub fn from_pdb(
        text: &str,
        templates: &[Component],
        selection: &Selection,
    ) -> Result<System, PdbError> {
        let mut by_code: BTreeMap<&str, &Component> = BTreeMap::new();
        for t in templates {
            if by_code.insert(t.id(), t).is_some() {
                return Err(PdbError::DuplicateTemplate {
                    template: t.id().into(),
                });
            }
        }
        let entry = read_entry(text)?;

        // Sort the residues into protein, ligand and dropped.
        let mut protein: Vec<&Raw> = Vec::new();
        let mut ligands: Vec<&Raw> = Vec::new();
        let mut dropped: Vec<Dropped> = Vec::new();
        let mut alternates_dropped = 0;
        for raw in &entry.residues {
            if raw.name == selection.ligand {
                ligands.push(raw);
                alternates_dropped += raw.alternates_dropped;
            } else if selection.drop.contains(&raw.name) {
                match dropped.iter_mut().find(|d| d.name == raw.name) {
                    Some(d) => {
                        d.residues += 1;
                        d.atoms += raw.records.len() + raw.alternates_dropped;
                    }
                    None => dropped.push(Dropped {
                        name: raw.name.clone(),
                        residues: 1,
                        atoms: raw.records.len() + raw.alternates_dropped,
                    }),
                }
            } else if raw.het {
                return Err(PdbError::Unexpected {
                    residue: raw.label(),
                });
            } else {
                if !AMINO_ACIDS.contains(&raw.name.as_str()) {
                    return Err(PdbError::NoTemplate {
                        residue: raw.label(),
                    });
                }
                protein.push(raw);
                alternates_dropped += raw.alternates_dropped;
            }
        }
        if protein.is_empty() {
            return Err(PdbError::NoProtein);
        }
        if ligands.len() != 1 {
            return Err(PdbError::LigandCount {
                ligand: selection.ligand.clone(),
                count: ligands.len(),
            });
        }
        let template = |raw: &Raw| {
            by_code
                .get(raw.name.as_str())
                .copied()
                .ok_or_else(|| PdbError::NoTemplate {
                    residue: raw.label(),
                })
        };

        // Chains, in the order they first appear.
        let mut chains: Vec<char> = Vec::new();
        for r in &protein {
            if !chains.contains(&r.chain) {
                chains.push(r.chain);
            }
        }
        let mut unmodelled = Vec::new();
        for &chain in &chains {
            let residues: Vec<&Raw> = protein
                .iter()
                .copied()
                .filter(|r| r.chain == chain)
                .collect();
            let seqres = entry
                .sequences
                .get(&chain)
                .ok_or(PdbError::NoSequence { chain })?;
            for r in &residues {
                let at = usize::try_from(r.number - 1)
                    .ok()
                    .and_then(|k| seqres.get(k));
                if at != Some(&r.name) {
                    return Err(PdbError::SequenceMismatch {
                        residue: r.label(),
                        seqres: at.cloned(),
                    });
                }
            }
            for pair in residues.windows(2) {
                let c = pair[0].position("C");
                let n = pair[1].position("N");
                let d = match (c, n) {
                    (Some(c), Some(n)) => distance(c, n),
                    _ => f64::NAN,
                };
                // A missing C or N leaves the distance NaN, refused here by name.
                if pair[1].number != pair[0].number + 1 || d.is_nan() || d > PEPTIDE_BOND_LIMIT {
                    if c.is_none() || n.is_none() {
                        let missing = if c.is_none() { pair[0] } else { pair[1] };
                        return Err(PdbError::MissingAtoms {
                            residue: missing.label(),
                            atoms: vec![if c.is_none() { "C" } else { "N" }.into()],
                        });
                    }
                    return Err(PdbError::ChainBreak {
                        after: pair[0].label(),
                        before: pair[1].label(),
                        distance: d,
                    });
                }
            }
            let (first, last) = (residues[0].number, residues[residues.len() - 1].number);
            for (k, name) in seqres.iter().enumerate() {
                let number = k as i32 + 1;
                if number < first || number > last {
                    unmodelled.push(Unmodelled {
                        chain,
                        number,
                        name: name.clone(),
                    });
                }
            }
        }

        // Disulfides, from the crystal's SG positions.
        let sulfurs: Vec<(usize, [f64; 3])> = protein
            .iter()
            .enumerate()
            .filter(|(_, r)| r.name == "CYS")
            .filter_map(|(k, r)| r.position("SG").map(|p| (k, p)))
            .collect();
        let mut partner: Vec<Option<usize>> = vec![None; protein.len()];
        let mut bridges: Vec<(usize, usize)> = Vec::new();
        for (a, &(i, pi)) in sulfurs.iter().enumerate() {
            for &(j, pj) in &sulfurs[a + 1..] {
                if distance(pi, pj) <= DISULFIDE_LIMIT {
                    for k in [i, j] {
                        if partner[k].is_some() {
                            return Err(PdbError::SulfurPartners {
                                residue: protein[k].label(),
                            });
                        }
                    }
                    partner[i] = Some(j);
                    partner[j] = Some(i);
                    bridges.push((i, j));
                }
            }
        }

        // Histidines: the rule, then the selection's overrides.
        let carboxylates: Vec<[f64; 3]> = protein
            .iter()
            .flat_map(|r| {
                let names: &[&str] = match r.name.as_str() {
                    "ASP" => &["OD1", "OD2"],
                    "GLU" => &["OE1", "OE2"],
                    _ => &[],
                };
                names
                    .iter()
                    .filter_map(|n| r.position(n))
                    .collect::<Vec<_>>()
            })
            .collect();
        let nearest = |p: Option<[f64; 3]>| {
            p.map_or(f64::INFINITY, |p| {
                carboxylates
                    .iter()
                    .map(|&o| distance(p, o))
                    .fold(f64::INFINITY, f64::min)
            })
        };
        let mut histidine: Vec<Option<Histidine>> = protein
            .iter()
            .map(|r| {
                (r.name == "HIS").then(|| {
                    let (d, e) = (nearest(r.position("ND1")), nearest(r.position("NE2")));
                    if d.min(e) <= HISTIDINE_ACCEPTOR_LIMIT {
                        Histidine::Both
                    } else {
                        Histidine::Epsilon
                    }
                })
            })
            .collect();
        for &(chain, number, state) in &selection.histidines {
            let k = protein
                .iter()
                .position(|r| r.chain == chain && r.number == number && r.name == "HIS")
                .ok_or_else(|| PdbError::NoSuchHistidine {
                    residue: format!("{chain}:{number}"),
                })?;
            histidine[k] = Some(state);
        }

        // Assemble.
        let mut out = Builder::default();
        let mut peptide_bonds = Vec::new();
        let mut start_of = vec![0; protein.len()];
        for &chain in &chains {
            let members: Vec<usize> = (0..protein.len())
                .filter(|&k| protein[k].chain == chain)
                .collect();
            let mut previous: Option<(usize, [f64; 3])> = None;
            for (m, &k) in members.iter().enumerate() {
                let raw = protein[k];
                let role = Role {
                    part: Part::Protein,
                    n_terminal: m == 0,
                    c_terminal: m + 1 == members.len(),
                    histidine: histidine[k],
                    disulfide: partner[k].is_some(),
                };
                let start = out.atoms.len();
                build(raw, template(raw)?, role, previous.map(|p| p.1), &mut out)?;
                let local = |name: &str| {
                    (start..out.atoms.len())
                        .find(|&i| out.local[i] == name)
                        .expect("a built residue has its backbone")
                };
                if let Some((c, _)) = previous {
                    peptide_bonds.push(out.bonds.len());
                    out.bonds.push(Bond {
                        atoms: [c, local("N")],
                        order: BondOrder::Single,
                        aromatic: false,
                    });
                }
                let c = local("C");
                previous = Some((c, angstrom(out.atoms[c].at)));
                start_of[k] = start;
            }
        }
        let mut disulfides = Vec::new();
        for (i, j) in bridges {
            let sg = |k: usize| {
                (start_of[k]..out.atoms.len())
                    .find(|&a| out.local[a] == "SG")
                    .expect("a built cysteine has its SG")
            };
            disulfides.push(out.bonds.len());
            out.bonds.push(Bond {
                atoms: [sg(i), sg(j)],
                order: BondOrder::Single,
                aromatic: false,
            });
        }
        let ligand = ligands[0];
        let role = Role {
            part: Part::Ligand,
            n_terminal: false,
            c_terminal: false,
            histidine: None,
            disulfide: false,
        };
        build(ligand, template(ligand)?, role, None, &mut out)?;
        turn_rotors(&mut out);

        let component = Component::assembled(
            entry.code.unwrap_or_else(|| "structure".into()),
            out.atoms,
            out.bonds,
        );
        Ok(System {
            component,
            residues: out.residues,
            residue_of: out.residue_of,
            local: out.local,
            placements: out.placements,
            turns: out.turns,
            peptide_bonds,
            disulfides,
            unmodelled,
            dropped,
            alternates_dropped,
        })
    }

    /// The whole system as one component, which [`crate::uff::assign`], [`ForceField::new`]
    /// and [`crate::qeq::charges`] take as they take a dictionary entry.
    ///
    /// [`ForceField::new`]: crate::energy::ForceField::new
    pub fn component(&self) -> &Component {
        &self.component
    }

    /// Every residue, protein chains first and the ligand last.
    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }

    /// The residue atom `i` belongs to.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn residue_of(&self, i: usize) -> &Residue {
        &self.residues[self.residue_of[i]]
    }

    /// Atom `i`'s name within its residue, `OD2`, which is its template's name for it.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn atom_name(&self, i: usize) -> &str {
        &self.local[i]
    }

    /// Which part atom `i` is in.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn part(&self, i: usize) -> Part {
        self.residue_of(i).part
    }

    /// The indices of every atom in `part`, ascending.
    pub fn atoms_in(&self, part: Part) -> Vec<usize> {
        (0..self.local.len())
            .filter(|&i| self.part(i) == part)
            .collect()
    }

    /// How each atom's position was obtained, in atom order.
    pub fn placements(&self) -> &[Placement] {
        &self.placements
    }

    /// How far each atom was turned about its rotor's bond from the dictionary conformer's
    /// torsion, radians, in atom order: zero for every atom not on a rotor, and for a rotor whose
    /// clearest step was the dictionary's. See the module documentation.
    pub fn turns(&self) -> &[f64] {
        &self.turns
    }

    /// The peptide bonds, as indices into the component's bonds.
    pub fn peptide_bonds(&self) -> &[usize] {
        &self.peptide_bonds
    }

    /// The disulfides, as indices into the component's bonds.
    pub fn disulfides(&self) -> &[usize] {
        &self.disulfides
    }

    /// The residues `SEQRES` lists that the file did not model, every chain's, in order.
    pub fn unmodelled(&self) -> &[Unmodelled] {
        &self.unmodelled
    }

    /// The `HETATM` residue types dropped, in the order first met.
    pub fn dropped(&self) -> &[Dropped] {
        &self.dropped
    }

    /// How many atom records of kept residues were dropped by the alternate-location rule.
    pub fn alternates_dropped(&self) -> usize {
        self.alternates_dropped
    }

    /// The sum of every atom's formal charge.
    pub fn formal_charge(&self) -> i32 {
        self.component.atoms().iter().map(|a| a.charge).sum()
    }
}

/// The same placement in a coordinate-free word, for a report.
impl fmt::Display for Placement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Placement::Crystal => f.write_str("crystal"),
            Placement::Template { rmsd, points } => {
                write!(f, "template, {points}-atom fit to {:.3} Å", rmsd / ANGSTROM)
            }
            Placement::Backbone => f.write_str("backbone bisector"),
            Placement::Completion => f.write_str("tetrahedral completion"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rotation_about(axis: [f64; 3], angle: f64) -> [[f64; 3]; 3] {
        let [x, y, z] = unit(axis);
        let (s, c) = angle.sin_cos();
        let t = 1.0 - c;
        [
            [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
            [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
            [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
        ]
    }

    fn turn(r: &[[f64; 3]; 3], p: [f64; 3]) -> [f64; 3] {
        [
            r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2],
            r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2],
            r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2],
        ]
    }

    const SET: [[f64; 3]; 5] = [
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [2.0, 1.4, 0.0],
        [3.4, 1.6, 0.7],
        [-0.6, 0.9, -1.2],
    ];

    /// A set moved by a known rotation and translation is superposed back exactly, the rotation
    /// recovered element by element — for a half turn too, where a quaternion's real part is zero.
    #[test]
    #[allow(clippy::needless_range_loop)]
    fn a_rigid_motion_is_recovered() {
        for (axis, angle) in [
            ([1.0, 2.0, 3.0], 0.7),
            ([0.0, 0.0, 1.0], std::f64::consts::PI),
            ([1.0, -1.0, 0.5], 2.9),
            ([0.3, 0.1, -0.2], 1e-3),
        ] {
            let r = rotation_about(axis, angle);
            let shift = [4.0, -2.0, 7.5];
            let target: Vec<[f64; 3]> = SET.iter().map(|&p| add(turn(&r, p), shift)).collect();
            let fit = superpose(&SET, &target);
            assert!(fit.rmsd < 1e-12, "{axis:?} {angle}: rmsd {}", fit.rmsd);
            for i in 0..3 {
                for j in 0..3 {
                    assert!(
                        (fit.rotation[i][j] - r[i][j]).abs() < 1e-12,
                        "{axis:?} {angle}: R[{i}][{j}] {} against {}",
                        fit.rotation[i][j],
                        r[i][j]
                    );
                }
            }
        }
    }

    /// A mirror image is not superposed: the fit stays a proper rotation, so a chiral set and its
    /// reflection keep a residual, and the rotation's determinant is +1.
    #[test]
    fn a_mirror_image_is_not_superposed() {
        let mirror: Vec<[f64; 3]> = SET.iter().map(|p| [p[0], p[1], -p[2]]).collect();
        let fit = superpose(&SET, &mirror);
        let r = fit.rotation;
        let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])
            - r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])
            + r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
        assert!((det - 1.0).abs() < 1e-12, "determinant {det}");
        assert!(fit.rmsd > 0.1, "rmsd {}", fit.rmsd);
    }
}
