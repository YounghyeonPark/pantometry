//! One chemical component, read out of the wwPDB Chemical Component Dictionary.
//!
//! The Chemical Component Dictionary (CCD) is the wwPDB's reference for every small molecule that
//! appears in the Protein Data Bank: one mmCIF data block per component, keyed by a code of up to
//! five characters — `AIN` is aspirin, `ATP` is ATP. Each block lists the atoms with their
//! elements, formal charges and aromaticity, an idealised geometry, and every bond with its order.
//! **Hydrogens are explicit**, which is what makes the file enough to type a molecule for a force
//! field without guessing.
//!
//! # What is read, and what is refused
//!
//! Two categories: `_chem_comp_atom` (atom name, element, formal charge, aromatic flag,
//! coordinates) and `_chem_comp_bond` (the two atom names, the order, the aromatic flag), plus
//! `_chem_comp.id`, `.name`, `.formula` and `.pdbx_formal_charge` when they are there.
//!
//! **The reader is strict, and every refusal names what it refused.** A missing column, an element
//! outside the ten this crate has parameters for, a bond order other than `SING`, `DOUB` or
//! `TRIP`, a bond naming an atom that is not in the atom list, two atoms with one name, a bond
//! given twice: each is a [`CcdError`] rather than a molecule. The silent alternatives are all
//! worse than a refusal — a bond to nowhere dropped, a `QUAD` read as single, an unknown element
//! typed as carbon — because each produces a molecule that looks fine and is typed wrong.
//!
//! # Coordinates: ideal, unless the file has none
//!
//! A CCD entry carries two geometries. `pdbx_model_Cartn_*_ideal` is a computed, idealised
//! conformation centred near the origin; `model_Cartn_*` is the component as it was first
//! observed in some deposited structure, in *that* crystal's frame. **The two are in different
//! frames** — for aspirin, about twenty ångström apart — so mixing them atom by atom would tear a
//! molecule in two without any single value looking wrong. So the choice is made for the whole component: ideal if every atom has
//! all three ideal coordinates, otherwise model for every atom, and [`Component::coordinates`]
//! says which. A component with neither complete is refused.
//!
//! # Why the parser takes text and not a path
//!
//! Nothing in `crates/` opens a file — every crate compiles to `wasm32`, where there is no
//! filesystem. [`Component::from_ccd`] takes the contents and the caller finds them.
//!
//! # The subset of CIF this reads
//!
//! Enough of CIF 1.1 to read any CCD file and to refuse anything else rather than misread it:
//! `data_` blocks (exactly one), `loop_` tables, tag–value pairs (which is how mmCIF writes a
//! category with a single row, so a one-atom component's atoms arrive that way), `#` comments,
//! single- and double-quoted values — whose closing quote is one *followed by whitespace*, so
//! `"C1'"` is the atom name `C1'` — and `;`-delimited text fields. `save_` frames and `global_`
//! are dictionary constructs no component file uses, and they are refused by name.

use std::collections::BTreeMap;
use std::fmt;

/// Ångström, in metres. CCD coordinates are in these; everything past the parser is SI.
pub const ANGSTROM: f64 = 1e-10;

/// The elements this crate can type and has parameters for.
///
/// Ten, because the Universal Force Field types below are written for them and nothing else. An
/// element outside this list is refused at parse time by name, rather than typed as something it
/// is not — a metal centre read as carbon would be a molecule that minimised to a confident wrong
/// answer.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Element {
    /// Hydrogen.
    H,
    /// Carbon.
    C,
    /// Nitrogen.
    N,
    /// Oxygen.
    O,
    /// Fluorine.
    F,
    /// Phosphorus.
    P,
    /// Sulfur.
    S,
    /// Chlorine.
    Cl,
    /// Bromine.
    Br,
    /// Iodine.
    I,
}

impl Element {
    /// Every element this crate knows, in atomic-number order.
    pub const ALL: [Element; 10] = [
        Element::H,
        Element::C,
        Element::N,
        Element::O,
        Element::F,
        Element::P,
        Element::S,
        Element::Cl,
        Element::Br,
        Element::I,
    ];

    /// The symbol as chemists write it: `C`, `Cl`, `Br`.
    pub fn symbol(self) -> &'static str {
        match self {
            Element::H => "H",
            Element::C => "C",
            Element::N => "N",
            Element::O => "O",
            Element::F => "F",
            Element::P => "P",
            Element::S => "S",
            Element::Cl => "Cl",
            Element::Br => "Br",
            Element::I => "I",
        }
    }

    /// The atomic number.
    pub fn atomic_number(self) -> u32 {
        match self {
            Element::H => 1,
            Element::C => 6,
            Element::N => 7,
            Element::O => 8,
            Element::F => 9,
            Element::P => 15,
            Element::S => 16,
            Element::Cl => 17,
            Element::Br => 35,
            Element::I => 53,
        }
    }

    /// The row of the periodic table, one for hydrogen to five for iodine.
    ///
    /// The Universal Force Field's sp² torsion barrier is tabulated by this and nothing else.
    pub fn period(self) -> u32 {
        match self {
            Element::H => 1,
            Element::C | Element::N | Element::O | Element::F => 2,
            Element::P | Element::S | Element::Cl => 3,
            Element::Br => 4,
            Element::I => 5,
        }
    }

    /// Read a CCD `type_symbol`, which the dictionary writes in upper case (`CL`, `BR`).
    ///
    /// Case-insensitive, because the case is a convention of the file and not of the element.
    pub fn from_symbol(symbol: &str) -> Option<Element> {
        Element::ALL
            .into_iter()
            .find(|e| e.symbol().eq_ignore_ascii_case(symbol))
    }
}

impl fmt::Display for Element {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.symbol())
    }
}

/// A bond's order as the dictionary states it.
///
/// **An aromatic ring arrives Kekulé.** The CCD writes benzene as alternating `SING` and `DOUB`
/// with every ring bond's aromatic flag set, rather than as six `AROM` bonds — so the order here
/// is one of the two resonance structures and [`Bond::aromatic`] is what says the ring is
/// delocalised. The typing rules read the flag first for exactly this reason.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    /// `SING`.
    Single,
    /// `DOUB`.
    Double,
    /// `TRIP`.
    Triple,
}

impl BondOrder {
    /// One, two or three.
    pub fn count(self) -> u32 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
        }
    }
}

/// Which of the file's two geometries the positions came from. See the module documentation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Coordinates {
    /// `pdbx_model_Cartn_*_ideal`: the dictionary's idealised conformation. The default.
    Ideal,
    /// `model_Cartn_*`: the component as observed in a deposited structure, used only because at
    /// least one atom had no ideal coordinate.
    Model,
}

/// One atom of a component.
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// The atom's name in the dictionary, `C7`, `HO1`, `C1'` — unique within a component, and
    /// what every bond refers to it by.
    pub name: String,
    /// Its element, from `type_symbol`.
    pub element: Element,
    /// The formal charge, in elementary charges.
    pub charge: i32,
    /// The atom's own aromatic flag. [`Bond::aromatic`] carries the same information per bond.
    pub aromatic: bool,
    /// Where it is, in metres.
    pub at: [f64; 3],
}

/// One bond of a component.
#[derive(Clone, Debug, PartialEq)]
pub struct Bond {
    /// The two atoms, as indices into [`Component::atoms`], in the order the file names them.
    pub atoms: [usize; 2],
    /// The order the file states. For an aromatic bond this is a Kekulé order — see [`BondOrder`].
    pub order: BondOrder,
    /// Whether the dictionary marks this bond aromatic.
    pub aromatic: bool,
}

/// Why a CCD file did not produce a component. Every variant names what it refused.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CcdError {
    /// A quoted value with no closing quote on its line.
    UnterminatedQuote {
        /// Which line, counting from one.
        line: usize,
    },
    /// A `;` text field with no closing `;` line.
    UnterminatedText {
        /// The line it opened on, counting from one.
        line: usize,
    },
    /// A value where a tag, `loop_` or `data_` was expected.
    StrayValue {
        /// Which line, counting from one.
        line: usize,
        /// What it said.
        text: String,
    },
    /// A tag with no value after it.
    MissingValue {
        /// The tag.
        tag: String,
    },
    /// The same tag given twice, as two pairs or as a pair and a loop column.
    DuplicateTag {
        /// The tag.
        tag: String,
    },
    /// A construct this reader does not implement — `save_`, `global_`, or a stop word.
    Unsupported {
        /// Which line, counting from one.
        line: usize,
        /// The keyword.
        text: String,
    },
    /// More than one `data_` block, or none. A component file holds one component; the
    /// dictionary as a whole holds tens of thousands, and taking the first of them silently would
    /// be reading a different molecule from the one asked for.
    DataBlocks {
        /// How many there were.
        count: usize,
    },
    /// A `loop_` whose values do not fill a whole number of rows, or that has no columns.
    LoopShape {
        /// The category, `_chem_comp_atom`.
        category: String,
        /// How many columns the loop declared.
        columns: usize,
        /// How many values followed.
        values: usize,
    },
    /// A category this reader needs is not in the file at all.
    MissingCategory {
        /// The category, `_chem_comp_atom`.
        category: String,
    },
    /// A category is there and one of the columns it needs is not.
    MissingColumn {
        /// The category, `_chem_comp_atom`.
        category: String,
        /// The column, `type_symbol`.
        column: String,
    },
    /// An element this crate has no Universal Force Field type for.
    UnknownElement {
        /// The atom, by name.
        atom: String,
        /// The symbol the file gave.
        symbol: String,
    },
    /// A bond order other than `SING`, `DOUB` or `TRIP` — `AROM`, `QUAD`, `DELO`, `POLY` or `PI`.
    UnknownBondOrder {
        /// The first atom, by name.
        first: String,
        /// The second atom, by name.
        second: String,
        /// The order the file gave.
        order: String,
    },
    /// A bond names an atom that is not in `_chem_comp_atom`.
    UnknownAtom {
        /// The name the bond used.
        atom: String,
        /// The first atom of that bond, as the file names it.
        first: String,
        /// The second atom of that bond, as the file names it.
        second: String,
    },
    /// Two atoms with one name. Every bond refers to atoms by name, so a repeated name makes
    /// every bond to it ambiguous.
    DuplicateAtom {
        /// The name.
        atom: String,
    },
    /// The same pair of atoms bonded twice.
    DuplicateBond {
        /// The first atom, by name.
        first: String,
        /// The second atom, by name.
        second: String,
    },
    /// A bond from an atom to itself.
    SelfBond {
        /// The atom, by name.
        atom: String,
    },
    /// A value that is not what its column holds: a coordinate that is not a number, a charge
    /// that is not an integer, a flag that is not `Y` or `N`.
    BadValue {
        /// The atom or bond it belongs to, as the file names it.
        at: String,
        /// The column.
        column: String,
        /// What it said.
        text: String,
    },
    /// An atom with neither a complete ideal coordinate nor a complete model one.
    NoCoordinates {
        /// The atom, by name.
        atom: String,
    },
    /// More than one atom and no `_chem_comp_bond` at all. Every atom would type as though it
    /// were bonded to nothing, which is a molecule that looks fine and is typed wrong everywhere.
    NoBonds {
        /// How many atoms there were.
        atoms: usize,
    },
}

impl fmt::Display for CcdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            CcdError::UnterminatedQuote { line } => {
                write!(f, "line {line}: a quoted value is never closed")
            }
            CcdError::UnterminatedText { line } => {
                write!(f, "line {line}: a ';' text field is never closed")
            }
            CcdError::StrayValue { line, text } => {
                write!(f, "line {line}: {text:?} is a value with no tag before it")
            }
            CcdError::MissingValue { tag } => write!(f, "{tag} has no value"),
            CcdError::DuplicateTag { tag } => write!(f, "{tag} is given twice"),
            CcdError::Unsupported { line, text } => write!(
                f,
                "line {line}: {text:?} is a CIF construct no component file uses, and this \
                 reader refuses it rather than guess"
            ),
            CcdError::DataBlocks { count } => write!(
                f,
                "{count} data blocks; a component file holds exactly one, and taking the first of \
                 several would read a molecule nobody asked for"
            ),
            CcdError::LoopShape {
                category,
                columns,
                values,
            } => write!(
                f,
                "the {category} loop has {columns} columns and {values} values, which is not a \
                 whole number of rows"
            ),
            CcdError::MissingCategory { category } => {
                write!(f, "the file has no {category} category")
            }
            CcdError::MissingColumn { category, column } => {
                write!(f, "{category} has no {column} column")
            }
            CcdError::UnknownElement { atom, symbol } => write!(
                f,
                "atom {atom} is element {symbol:?}, which has no Universal Force Field type here \
                 (H, C, N, O, F, P, S, Cl, Br and I do)"
            ),
            CcdError::UnknownBondOrder {
                first,
                second,
                order,
            } => write!(
                f,
                "the bond {first}-{second} has order {order:?}; only SING, DOUB and TRIP are read"
            ),
            CcdError::UnknownAtom {
                atom,
                first,
                second,
            } => write!(
                f,
                "the bond {first}-{second} names atom {atom}, which is not in _chem_comp_atom"
            ),
            CcdError::DuplicateAtom { atom } => {
                write!(
                    f,
                    "two atoms are named {atom}, so a bond to it is ambiguous"
                )
            }
            CcdError::DuplicateBond { first, second } => {
                write!(f, "the bond {first}-{second} is given twice")
            }
            CcdError::SelfBond { atom } => write!(f, "atom {atom} is bonded to itself"),
            CcdError::BadValue { at, column, text } => {
                write!(f, "{at}: {column} is {text:?}, which it cannot be")
            }
            CcdError::NoCoordinates { atom } => write!(
                f,
                "atom {atom} has neither a complete ideal coordinate nor a complete model one"
            ),
            CcdError::NoBonds { atoms } => write!(
                f,
                "{atoms} atoms and no _chem_comp_bond; every atom would be typed as bonded to \
                 nothing"
            ),
        }
    }
}

impl std::error::Error for CcdError {}

/// A chemical component: its atoms, its bonds, and what the file says about it as a whole.
#[derive(Clone, Debug, PartialEq)]
pub struct Component {
    id: String,
    name: Option<String>,
    formula: Option<String>,
    formal_charge: Option<i32>,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    coordinates: Coordinates,
}

impl Component {
    /// Read one component from the text of a CCD mmCIF file.
    ///
    /// # Errors
    ///
    /// Any [`CcdError`], each naming what was refused. See the module documentation for the list
    /// and for why none of them is recovered from.
    pub fn from_ccd(text: &str) -> Result<Component, CcdError> {
        let tables = read_tables(text)?;

        let header = tables
            .get("_chem_comp")
            .ok_or_else(|| CcdError::MissingCategory {
                category: "_chem_comp".into(),
            })?;
        let id =
            header
                .single("id")
                .map(str::to_string)
                .ok_or_else(|| CcdError::MissingColumn {
                    category: "_chem_comp".into(),
                    column: "id".into(),
                })?;
        let name = header.single("name").map(str::to_string);
        let formula = header.single("formula").map(str::to_string);
        let formal_charge = match header.single("pdbx_formal_charge") {
            None => None,
            Some(text) => Some(text.parse::<i32>().map_err(|_| CcdError::BadValue {
                at: id.clone(),
                column: "pdbx_formal_charge".into(),
                text: text.into(),
            })?),
        };

        let (atoms, coordinates) = read_atoms(tables.get("_chem_comp_atom"))?;
        let bonds = match tables.get("_chem_comp_bond") {
            Some(table) => read_bonds(table, &atoms)?,
            None if atoms.len() > 1 => return Err(CcdError::NoBonds { atoms: atoms.len() }),
            None => Vec::new(),
        };

        Ok(Component {
            id,
            name,
            formula,
            formal_charge,
            atoms,
            bonds,
            coordinates,
        })
    }

    /// The component's code in the dictionary, `AIN`.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// `_chem_comp.name`, when the file gives one: `2-(ACETYLOXY)BENZOIC ACID`.
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    /// `_chem_comp.formula` as the file writes it, `C9 H8 O4`. **Not computed from the atoms**:
    /// it is a second statement in the same file, which is what makes it worth comparing against
    /// [`Component::composition`].
    pub fn formula(&self) -> Option<&str> {
        self.formula.as_deref()
    }

    /// `_chem_comp.pdbx_formal_charge`, the file's statement of the net charge.
    pub fn formal_charge(&self) -> Option<i32> {
        self.formal_charge
    }

    /// The atoms, in the file's order.
    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    /// The bonds, in the file's order.
    pub fn bonds(&self) -> &[Bond] {
        &self.bonds
    }

    /// Which geometry the positions are.
    pub fn coordinates(&self) -> Coordinates {
        self.coordinates
    }

    /// How many atoms of each element, counted from the atoms.
    pub fn composition(&self) -> BTreeMap<Element, usize> {
        let mut count = BTreeMap::new();
        for atom in &self.atoms {
            *count.entry(atom.element).or_insert(0) += 1;
        }
        count
    }

    /// The bonds atom `i` takes part in, as `(the other atom, the bond)`.
    pub fn neighbours(&self, i: usize) -> impl Iterator<Item = (usize, &Bond)> {
        self.bonds.iter().filter_map(move |b| match b.atoms {
            [a, o] if a == i => Some((o, b)),
            [o, a] if a == i => Some((o, b)),
            _ => None,
        })
    }
}

/// One token of CIF, and whether it was quoted — `?` quoted is the string `?` and not "unknown".
#[derive(Clone, Debug)]
struct Token {
    text: String,
    quoted: bool,
    line: usize,
}

impl Token {
    fn is_tag(&self) -> bool {
        !self.quoted && self.text.starts_with('_')
    }

    fn keyword(&self) -> Option<String> {
        if self.quoted {
            return None;
        }
        let lower = self.text.to_ascii_lowercase();
        ["data_", "loop_", "save_", "global_", "stop_"]
            .iter()
            .any(|k| lower.starts_with(k))
            .then_some(lower)
    }

    /// `?` (unknown) and `.` (inapplicable), unquoted.
    fn is_null(&self) -> bool {
        !self.quoted && (self.text == "?" || self.text == ".")
    }
}

fn tokenise(text: &str) -> Result<Vec<Token>, CcdError> {
    let mut out = Vec::new();
    let lines: Vec<&str> = text.lines().collect();
    let mut n = 0;
    while n < lines.len() {
        let line = lines[n];
        let number = n + 1;
        if let Some(first) = line.strip_prefix(';') {
            // A text field: everything up to a line that begins with `;`.
            let mut body = vec![first];
            let mut end = None;
            for (k, l) in lines.iter().enumerate().skip(n + 1) {
                if l.starts_with(';') {
                    end = Some(k);
                    break;
                }
                body.push(l);
            }
            let Some(end) = end else {
                return Err(CcdError::UnterminatedText { line: number });
            };
            out.push(Token {
                text: body.join("\n").trim().to_string(),
                quoted: true,
                line: number,
            });
            // Whatever follows the closing `;` on its line is tokenised as ordinary text.
            words(&lines[end][1..], end + 1, &mut out)?;
            n = end + 1;
            continue;
        }
        words(line, number, &mut out)?;
        n += 1;
    }
    Ok(out)
}

/// The whitespace-separated tokens of one line, honouring quotes and comments.
fn words(line: &str, number: usize, out: &mut Vec<Token>) -> Result<(), CcdError> {
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        if c.is_whitespace() {
            i += 1;
            continue;
        }
        if c == '#' {
            return Ok(());
        }
        if c == '\'' || c == '"' {
            // A closing quote is the same character followed by whitespace or the end of the
            // line, so `"C1'"` holds an apostrophe and `'it's'` is one value.
            let start = i + 1;
            let mut j = start;
            let close = loop {
                if j >= chars.len() {
                    return Err(CcdError::UnterminatedQuote { line: number });
                }
                if chars[j] == c && (j + 1 == chars.len() || chars[j + 1].is_whitespace()) {
                    break j;
                }
                j += 1;
            };
            out.push(Token {
                text: chars[start..close].iter().collect(),
                quoted: true,
                line: number,
            });
            i = close + 1;
            continue;
        }
        let start = i;
        while i < chars.len() && !chars[i].is_whitespace() {
            i += 1;
        }
        out.push(Token {
            text: chars[start..i].iter().collect(),
            quoted: false,
            line: number,
        });
    }
    Ok(())
}

/// One category's values: column names (lower case, without the category) and rows.
#[derive(Debug)]
struct Table {
    columns: Vec<String>,
    rows: Vec<Vec<Token>>,
}

impl Table {
    fn column(&self, category: &str, name: &str) -> Result<usize, CcdError> {
        self.find(name).ok_or_else(|| CcdError::MissingColumn {
            category: category.into(),
            column: name.into(),
        })
    }

    fn find(&self, name: &str) -> Option<usize> {
        let name = name.to_ascii_lowercase();
        self.columns.iter().position(|c| *c == name)
    }

    /// A value from a one-row category, or `None` when the column is absent or null.
    fn single(&self, name: &str) -> Option<&str> {
        let k = self.find(name)?;
        let token = self.rows.first()?.get(k)?;
        (!token.is_null()).then_some(token.text.as_str())
    }
}

/// Split `_chem_comp_atom.type_symbol` into the category and the column, both lower case.
fn split_tag(tag: &str) -> (String, String) {
    let lower = tag.to_ascii_lowercase();
    match lower.split_once('.') {
        Some((cat, col)) => (cat.to_string(), col.to_string()),
        None => (lower.clone(), String::new()),
    }
}

fn read_tables(text: &str) -> Result<BTreeMap<String, Table>, CcdError> {
    let tokens = tokenise(text)?;
    let mut tables: BTreeMap<String, Table> = BTreeMap::new();
    let mut blocks = 0;
    let mut i = 0;
    while i < tokens.len() {
        let token = &tokens[i];
        if let Some(word) = token.keyword() {
            if word.starts_with("data_") {
                blocks += 1;
                i += 1;
                continue;
            }
            if word != "loop_" {
                return Err(CcdError::Unsupported {
                    line: token.line,
                    text: token.text.clone(),
                });
            }
            // The loop's columns, then its values up to the next tag or keyword.
            i += 1;
            let mut tags = Vec::new();
            while i < tokens.len() && tokens[i].is_tag() {
                tags.push(split_tag(&tokens[i].text));
                i += 1;
            }
            let mut values = Vec::new();
            while i < tokens.len() && !tokens[i].is_tag() && tokens[i].keyword().is_none() {
                values.push(tokens[i].clone());
                i += 1;
            }
            let category = tags.first().map(|t| t.0.clone()).unwrap_or_default();
            if tags.is_empty() || values.len() % tags.len() != 0 {
                return Err(CcdError::LoopShape {
                    category,
                    columns: tags.len(),
                    values: values.len(),
                });
            }
            // A loop holds one category; a column from another would be a malformed file.
            if let Some((cat, col)) = tags.iter().find(|t| t.0 != category) {
                return Err(CcdError::DuplicateTag {
                    tag: format!("{cat}.{col}"),
                });
            }
            if tables.contains_key(&category) {
                return Err(CcdError::DuplicateTag { tag: category });
            }
            let columns: Vec<String> = tags.into_iter().map(|t| t.1).collect();
            for (k, c) in columns.iter().enumerate() {
                if columns[..k].contains(c) {
                    return Err(CcdError::DuplicateTag {
                        tag: format!("{category}.{c}"),
                    });
                }
            }
            let rows = values
                .chunks(columns.len())
                .map(<[Token]>::to_vec)
                .collect();
            tables.insert(category, Table { columns, rows });
            continue;
        }
        if token.is_tag() {
            let value = tokens
                .get(i + 1)
                .filter(|v| !v.is_tag() && v.keyword().is_none())
                .ok_or_else(|| CcdError::MissingValue {
                    tag: token.text.clone(),
                })?;
            let (category, column) = split_tag(&token.text);
            let table = tables.entry(category.clone()).or_insert(Table {
                columns: Vec::new(),
                rows: vec![Vec::new()],
            });
            // A category already read as a loop has more than one row; a pair cannot extend it.
            if table.rows.len() != 1 || table.find(&column).is_some() {
                return Err(CcdError::DuplicateTag {
                    tag: token.text.clone(),
                });
            }
            table.columns.push(column);
            table.rows[0].push(value.clone());
            i += 2;
            continue;
        }
        return Err(CcdError::StrayValue {
            line: token.line,
            text: token.text.clone(),
        });
    }
    if blocks != 1 {
        return Err(CcdError::DataBlocks { count: blocks });
    }
    Ok(tables)
}

fn flag(token: &Token, at: &str, column: &str) -> Result<bool, CcdError> {
    match token.text.as_str() {
        "Y" | "y" => Ok(true),
        "N" | "n" => Ok(false),
        other => Err(CcdError::BadValue {
            at: at.into(),
            column: column.into(),
            text: other.into(),
        }),
    }
}

/// Three coordinates, in metres, or `None` if any of them is null.
fn point(
    row: &[Token],
    cols: [usize; 3],
    at: &str,
    prefix: &str,
) -> Result<Option<[f64; 3]>, CcdError> {
    let mut p = [0.0; 3];
    for (k, &c) in cols.iter().enumerate() {
        let token = &row[c];
        if token.is_null() {
            return Ok(None);
        }
        let v: f64 = token
            .text
            .parse()
            .ok()
            .filter(|v: &f64| v.is_finite())
            .ok_or_else(|| CcdError::BadValue {
                at: at.into(),
                column: format!("{prefix}{}", ["x", "y", "z"][k]),
                text: token.text.clone(),
            })?;
        p[k] = v * ANGSTROM;
    }
    Ok(Some(p))
}

fn read_atoms(table: Option<&Table>) -> Result<(Vec<Atom>, Coordinates), CcdError> {
    const CAT: &str = "_chem_comp_atom";
    let table = table.ok_or_else(|| CcdError::MissingCategory {
        category: CAT.into(),
    })?;
    let id = table.column(CAT, "atom_id")?;
    let symbol = table.column(CAT, "type_symbol")?;
    let charge = table.column(CAT, "charge")?;
    let aromatic = table.column(CAT, "pdbx_aromatic_flag")?;
    let ideal = [
        table.column(CAT, "pdbx_model_Cartn_x_ideal")?,
        table.column(CAT, "pdbx_model_Cartn_y_ideal")?,
        table.column(CAT, "pdbx_model_Cartn_z_ideal")?,
    ];

    let mut atoms = Vec::with_capacity(table.rows.len());
    let mut ideal_points = Vec::with_capacity(table.rows.len());
    for row in &table.rows {
        let name = row[id].text.clone();
        if row[id].is_null() {
            return Err(CcdError::BadValue {
                at: CAT.into(),
                column: "atom_id".into(),
                text: name,
            });
        }
        if atoms.iter().any(|a: &Atom| a.name == name) {
            return Err(CcdError::DuplicateAtom { atom: name });
        }
        let element =
            Element::from_symbol(&row[symbol].text).ok_or_else(|| CcdError::UnknownElement {
                atom: name.clone(),
                symbol: row[symbol].text.clone(),
            })?;
        let charge = row[charge]
            .text
            .parse::<i32>()
            .map_err(|_| CcdError::BadValue {
                at: name.clone(),
                column: "charge".into(),
                text: row[charge].text.clone(),
            })?;
        let aromatic = flag(&row[aromatic], &name, "pdbx_aromatic_flag")?;
        ideal_points.push(point(row, ideal, &name, "pdbx_model_Cartn_")?);
        atoms.push(Atom {
            name,
            element,
            charge,
            aromatic,
            at: [0.0; 3],
        });
    }

    // Ideal for every atom, or model for every atom: the two are in different frames.
    let coordinates = if ideal_points.iter().all(Option::is_some) {
        for (atom, p) in atoms.iter_mut().zip(&ideal_points) {
            atom.at = p.expect("checked just above");
        }
        Coordinates::Ideal
    } else {
        let model = [
            table.column(CAT, "model_Cartn_x")?,
            table.column(CAT, "model_Cartn_y")?,
            table.column(CAT, "model_Cartn_z")?,
        ];
        for (atom, row) in atoms.iter_mut().zip(&table.rows) {
            atom.at = point(row, model, &atom.name, "model_Cartn_")?.ok_or_else(|| {
                CcdError::NoCoordinates {
                    atom: atom.name.clone(),
                }
            })?;
        }
        Coordinates::Model
    };
    Ok((atoms, coordinates))
}

fn read_bonds(table: &Table, atoms: &[Atom]) -> Result<Vec<Bond>, CcdError> {
    const CAT: &str = "_chem_comp_bond";
    let first = table.column(CAT, "atom_id_1")?;
    let second = table.column(CAT, "atom_id_2")?;
    let order = table.column(CAT, "value_order")?;
    let aromatic = table.column(CAT, "pdbx_aromatic_flag")?;

    let mut bonds: Vec<Bond> = Vec::with_capacity(table.rows.len());
    for row in &table.rows {
        let (a, b) = (row[first].text.clone(), row[second].text.clone());
        let find = |name: &str| {
            atoms
                .iter()
                .position(|atom| atom.name == name)
                .ok_or_else(|| CcdError::UnknownAtom {
                    atom: name.into(),
                    first: a.clone(),
                    second: b.clone(),
                })
        };
        let (i, j) = (find(&a)?, find(&b)?);
        if i == j {
            return Err(CcdError::SelfBond { atom: a });
        }
        if bonds.iter().any(|x| x.atoms == [i, j] || x.atoms == [j, i]) {
            return Err(CcdError::DuplicateBond {
                first: a,
                second: b,
            });
        }
        let value = match row[order].text.to_ascii_uppercase().as_str() {
            "SING" => BondOrder::Single,
            "DOUB" => BondOrder::Double,
            "TRIP" => BondOrder::Triple,
            _ => {
                return Err(CcdError::UnknownBondOrder {
                    first: a,
                    second: b,
                    order: row[order].text.clone(),
                })
            }
        };
        let label = format!("bond {a}-{b}");
        bonds.push(Bond {
            atoms: [i, j],
            order: value,
            aromatic: flag(&row[aromatic], &label, "pdbx_aromatic_flag")?,
        });
    }
    Ok(bonds)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texts(line: &str) -> Vec<String> {
        let mut out = Vec::new();
        words(line, 1, &mut out).unwrap();
        out.into_iter().map(|t| t.text).collect()
    }

    /// The closing quote is one followed by whitespace, so a prime inside a name survives.
    #[test]
    fn a_quote_inside_a_quoted_value_is_part_of_it() {
        assert_eq!(texts(r#"AIN "C1'" 'it's' x"#), ["AIN", "C1'", "it's", "x"]);
        assert_eq!(texts("a b # c d"), ["a", "b"]);
        let mut out = Vec::new();
        assert_eq!(
            words("x 'never closed", 7, &mut out),
            Err(CcdError::UnterminatedQuote { line: 7 })
        );
    }

    /// A `;` text field is one value, however many lines it spans.
    #[test]
    fn a_text_field_is_one_value() {
        let tokens = tokenise("_a.b\n;first\nsecond\n;\n_a.c x\n").unwrap();
        let t: Vec<&str> = tokens.iter().map(|t| t.text.as_str()).collect();
        assert_eq!(t, ["_a.b", "first\nsecond", "_a.c", "x"]);
        assert_eq!(
            tokenise("_a.b\n;never\nclosed\n").unwrap_err(),
            CcdError::UnterminatedText { line: 2 }
        );
    }
}
