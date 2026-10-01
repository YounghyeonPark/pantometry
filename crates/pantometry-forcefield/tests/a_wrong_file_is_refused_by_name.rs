//! **A file that is wrong is refused, and the refusal names what is wrong.**
//!
//! Each test changes one thing in aspirin's real entry and asserts the specific error — not merely
//! that parsing failed, because a parser that refused everything would pass a test that only asked
//! for an error. And each asserts the message names the thing, because "parse error" is not
//! something a person can act on.
//!
//! Every edit checks its anchor occurs exactly once first: a `replace` whose anchor did not match
//! returns the file unchanged, the parse succeeds, and the test fails for the wrong reason — or,
//! written carelessly, passes.

use pantometry_forcefield::{CcdError, Component};

const AIN: &str = include_str!("../components/AIN.cif");

fn edit(from: &str, to: &str) -> String {
    assert_eq!(
        AIN.matches(from).count(),
        1,
        "anchor {from:?} must occur once"
    );
    AIN.replacen(from, to, 1)
}

fn refused(text: &str) -> CcdError {
    Component::from_ccd(text).expect_err("this file should be refused")
}

fn names(e: &CcdError, what: &str) {
    let message = e.to_string();
    assert!(message.contains(what), "{message:?} does not name {what:?}");
}

/// **A missing column**: `type_symbol` taken out of the header *and* out of every row, so the
/// table is well-formed and simply lacks it.
#[test]
fn a_missing_column_is_named() {
    let mut out = Vec::new();
    let mut in_atoms = false;
    let mut rows = 0;
    for line in AIN.lines() {
        if line.trim() == "_chem_comp_atom.type_symbol" {
            in_atoms = true;
            continue;
        }
        if in_atoms && line.starts_with("AIN ") {
            let mut fields: Vec<&str> = line.split_whitespace().collect();
            assert_eq!(fields[3].len(), 1, "the fourth field is the element");
            fields.remove(3);
            out.push(fields.join(" "));
            rows += 1;
            continue;
        }
        if in_atoms && line.starts_with('#') {
            in_atoms = false;
        }
        out.push(line.to_string());
    }
    assert_eq!(rows, 21, "every atom row was edited");
    let e = refused(&out.join("\n"));
    assert_eq!(
        e,
        CcdError::MissingColumn {
            category: "_chem_comp_atom".into(),
            column: "type_symbol".into()
        }
    );
    names(&e, "type_symbol");
}

/// **A bond to an atom that is not there.** Dropped silently, this would leave `O1` one bond short
/// and type it as nothing worse than it already is — a molecule that looks fine.
#[test]
fn a_bond_to_an_unknown_atom_is_named() {
    let e = refused(&edit("AIN O1 HO1 SING", "AIN O1 HX1 SING"));
    assert_eq!(
        e,
        CcdError::UnknownAtom {
            atom: "HX1".into(),
            first: "O1".into(),
            second: "HX1".into()
        }
    );
    names(&e, "HX1");
}

/// **A bond order this crate does not read.** `QUAD` read as single would be the silent version.
#[test]
fn an_unknown_bond_order_is_named() {
    let e = refused(&edit("AIN C7 O2  DOUB", "AIN C7 O2  QUAD"));
    assert_eq!(
        e,
        CcdError::UnknownBondOrder {
            first: "C7".into(),
            second: "O2".into(),
            order: "QUAD".into()
        }
    );
    names(&e, "QUAD");
    // `AROM` too: the dictionary does not use it, and a file that did would be typed from an
    // order this crate has no rule for.
    let e = refused(&edit("AIN C7 O2  DOUB", "AIN C7 O2  AROM"));
    names(&e, "AROM");
}

/// **An element with no UFF type here.** Typed as carbon, an iron would minimise to a confident
/// wrong answer.
#[test]
fn an_unknown_element_is_named() {
    let e = refused(&edit("AIN O1  O1  O 0", "AIN O1  O1  FE 0"));
    assert_eq!(
        e,
        CcdError::UnknownElement {
            atom: "O1".into(),
            symbol: "FE".into()
        }
    );
    names(&e, "FE");
}

/// **Two atoms with one name** make every bond to that name ambiguous.
#[test]
fn a_repeated_atom_name_is_named() {
    let e = refused(&edit("AIN H5  H5  H", "AIN H4  H5  H"));
    assert_eq!(e, CcdError::DuplicateAtom { atom: "H4".into() });
}

/// **A bond given twice**, and **a bond from an atom to itself.**
#[test]
fn a_repeated_or_reflexive_bond_is_named() {
    let e = refused(&edit("AIN C9 H93 SING", "AIN H91 C9 SING"));
    assert_eq!(
        e,
        CcdError::DuplicateBond {
            first: "H91".into(),
            second: "C9".into()
        }
    );
    let e = refused(&edit("AIN C9 H93 SING", "AIN C9 C9 SING"));
    assert_eq!(e, CcdError::SelfBond { atom: "C9".into() });
}

/// **No coordinates at all for one atom** — neither ideal nor model — is refused rather than
/// placed at the origin.
#[test]
fn an_atom_with_no_coordinates_is_named() {
    let e = refused(&edit(
        "13.394 16.144 -0.175 2.659  0.080  -3.183",
        "?      ?      ?      ?      ?      ?",
    ));
    assert_eq!(e, CcdError::NoCoordinates { atom: "HO1".into() });
}

/// **A row one value short** makes the whole loop the wrong shape, rather than shifting every
/// value after it into the wrong column.
#[test]
fn a_short_row_is_refused() {
    let e = refused(&edit("AIN C9 H93 SING N N 21", "AIN C9 H93 SING N N"));
    assert_eq!(
        e,
        CcdError::LoopShape {
            category: "_chem_comp_bond".into(),
            columns: 7,
            values: 7 * 21 - 1
        }
    );
}

/// **A charge that is not an integer** and **a flag that is not `Y` or `N`.**
#[test]
fn a_value_that_is_not_what_its_column_holds_is_named() {
    let e = refused(&edit("AIN O1  O1  O 0", "AIN O1  O1  O ?"));
    assert_eq!(
        e,
        CcdError::BadValue {
            at: "O1".into(),
            column: "charge".into(),
            text: "?".into()
        }
    );
    let e = refused(&edit("AIN C3 C4  SING Y", "AIN C3 C4  SING maybe"));
    names(&e, "maybe");
}

/// **No bonds at all** for twenty-one atoms is refused: every atom would type as bonded to
/// nothing — every carbon `C_3` — and nothing would look wrong.
#[test]
fn atoms_without_a_bond_table_are_refused() {
    let start = AIN.find("loop_\n_chem_comp_bond.").expect("the bond loop");
    let end = start + AIN[start..].find("\n# \n").expect("its end");
    let text = format!("{}{}", &AIN[..start], &AIN[end..]);
    let e = refused(&text);
    assert_eq!(e, CcdError::NoBonds { atoms: 21 });
}

/// **No atoms at all**, and **two components in one file.** The dictionary as a whole is one file
/// of tens of thousands of blocks; taking the first would read a molecule nobody asked for.
#[test]
fn a_file_with_no_atoms_or_two_components_is_refused() {
    let e = refused("data_X\n_chem_comp.id X\n");
    assert_eq!(
        e,
        CcdError::MissingCategory {
            category: "_chem_comp_atom".into()
        }
    );
    let e = refused(&format!("{AIN}\ndata_AIN2\n"));
    assert_eq!(e, CcdError::DataBlocks { count: 2 });
    names(&e, "2 data blocks");
}
