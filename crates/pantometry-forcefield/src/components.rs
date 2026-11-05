//! The wwPDB Chemical Component Dictionary entries this crate ships, by their identifiers.
//!
//! ```
//! use pantometry_forcefield::{components, Component};
//!
//! let benzene = Component::from_ccd(components::text("BNZ").unwrap()).unwrap();
//! assert_eq!(benzene.atoms().len(), 12);
//! assert!(components::text("XYZ").is_none());
//! ```
//!
//! Each is the entry as RCSB serves it, byte for byte, and `components/README.md` in the
//! repository gives each one's source, date and SHA-256. They are what this crate's tests are
//! checked against, so a consumer can run the same molecules. The texts are compiled in, and a
//! program that never calls [`text`] or reads [`ALL`] does not keep them.
//!
//! **Licence.** wwPDB data files are under the CC0 1.0 Universal public domain dedication
//! (wwpdb.org/about/usage-policies), so shipping them puts no condition on this crate's own
//! MIT OR Apache-2.0.
//!
//! **The list is the directory.** `every_shipped_component_is_listed` reads `components/` and
//! fails if a `.cif` there is missing here, or listed here and missing there.

/// Every entry, as `(identifier, text)`, in the identifiers' byte order.
pub const ALL: [(&str, &str); 52] = [
    ("16R", include_str!("../components/16R.cif")),
    ("2F2", include_str!("../components/2F2.cif")),
    ("2ME", include_str!("../components/2ME.cif")),
    ("61G", include_str!("../components/61G.cif")),
    ("6AC", include_str!("../components/6AC.cif")),
    ("A1JFW", include_str!("../components/A1JFW.cif")),
    ("ACE", include_str!("../components/ACE.cif")),
    ("ACM", include_str!("../components/ACM.cif")),
    ("ACN", include_str!("../components/ACN.cif")),
    ("AIN", include_str!("../components/AIN.cif")),
    ("ALA", include_str!("../components/ALA.cif")),
    ("ARG", include_str!("../components/ARG.cif")),
    ("ASN", include_str!("../components/ASN.cif")),
    ("ASP", include_str!("../components/ASP.cif")),
    ("BNZ", include_str!("../components/BNZ.cif")),
    ("BZF", include_str!("../components/BZF.cif")),
    ("CCN", include_str!("../components/CCN.cif")),
    ("CYS", include_str!("../components/CYS.cif")),
    ("D10", include_str!("../components/D10.cif")),
    ("DEN", include_str!("../components/DEN.cif")),
    ("DMF", include_str!("../components/DMF.cif")),
    ("DMN", include_str!("../components/DMN.cif")),
    ("GLN", include_str!("../components/GLN.cif")),
    ("GLU", include_str!("../components/GLU.cif")),
    ("GLY", include_str!("../components/GLY.cif")),
    ("HIS", include_str!("../components/HIS.cif")),
    ("I4B", include_str!("../components/I4B.cif")),
    ("ILE", include_str!("../components/ILE.cif")),
    ("IND", include_str!("../components/IND.cif")),
    ("KEN", include_str!("../components/KEN.cif")),
    ("LEU", include_str!("../components/LEU.cif")),
    ("LYS", include_str!("../components/LYS.cif")),
    ("MEE", include_str!("../components/MEE.cif")),
    ("MET", include_str!("../components/MET.cif")),
    ("MOH", include_str!("../components/MOH.cif")),
    ("N4B", include_str!("../components/N4B.cif")),
    ("NME", include_str!("../components/NME.cif")),
    ("OC9", include_str!("../components/OC9.cif")),
    ("OXE", include_str!("../components/OXE.cif")),
    ("PEO", include_str!("../components/PEO.cif")),
    ("PHE", include_str!("../components/PHE.cif")),
    ("PRO", include_str!("../components/PRO.cif")),
    ("PXY", include_str!("../components/PXY.cif")),
    ("PYJ", include_str!("../components/PYJ.cif")),
    ("R16", include_str!("../components/R16.cif")),
    ("S2H", include_str!("../components/S2H.cif")),
    ("SER", include_str!("../components/SER.cif")),
    ("THR", include_str!("../components/THR.cif")),
    ("TME", include_str!("../components/TME.cif")),
    ("TRP", include_str!("../components/TRP.cif")),
    ("TYR", include_str!("../components/TYR.cif")),
    ("VAL", include_str!("../components/VAL.cif")),
];

/// The text of entry `id`, as the dictionary writes the identifier: `"BNZ"`, `"AIN"`, `"16R"`.
/// Exact, not case-folded, because the dictionary's identifiers are upper case and a lower-case
/// one is more likely a different name than the same one.
pub fn text(id: &str) -> Option<&'static str> {
    ALL.iter().find(|(i, _)| *i == id).map(|(_, t)| *t)
}
