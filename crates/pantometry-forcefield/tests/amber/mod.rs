//! **ff14SB's protein Lennard-Jones, by AMBER atom type**, for W4's decomposition: which type each
//! of 181L's atoms is, from its residue's template, and that type's `R*/2` and `ε`. Data, cited,
//! not code: the type columns of the residue libraries and the van der Waals lines of the
//! parameter files, transcribed; only the lookup is written here.
//!
//! - **The types** are the `str type` column of `amino12.lib`, `aminont12.lib` (for the
//!   N-terminal `NMET`) and `aminoct12.lib` (for the C-terminal `CLYS`), the three protein
//!   libraries `leaprc.ff14SB` loads, in ambermini at commit `f7c421c4` ([`SOURCES`] has every path
//!   and SHA-256 in full). Only the templates 181L needs: nineteen amino acids, His31 as `HIP` (the
//!   crate's [`Histidine::Both`]), and the two termini. **The dictionary's atom names are AMBER's**
//!   (PDB version 3) but for the N-terminal hydrogens, which the crate names `H`, `H2`, `H3` and
//!   `NMET` names `H1`, `H2`, `H3`; all three are type `H`, and the table is written in the crate's
//!   names. The 356 atoms of these 21 templates were each read back against OpenMM's
//!   `amber14/protein.ff14SB.xml` by a script outside the repository, and none differs; the atom
//!   names of every residue of 181L are exactly the template's.
//! - **The parameters** are `parm10.dat`'s `MOD4 RE` block and `frcmod.ff14SB`'s `NONB` block at
//!   the same commit, with `parm10.dat`'s two equivalence lines (`N NA N2 …` and `C* CA CB CC …`):
//!   a type with no line of its own takes its equivalent's ([`EQUIVALENT`]). All thirty types here
//!   agree with OpenMM's `protein.ff14SB.xml` through `R*/2 = 2^⅙ σ/2` and 4.184 kJ/kcal, to its
//!   printed digits (`HO`'s ε is zero in both).
//!
//! ff14SB splits parm99's `CT` into `CX` (every α carbon), `2C`, `3C` and `C8` (side-chain sp³
//! carbons with two or three heavy neighbours, and Arg's and Lys's), and gives the carboxylate
//! carbon `CO`. **Every one of them carries `CT`'s or `C`'s Lennard-Jones**, so the split moves no
//! energy here, and the decomposition by type keeps them apart only as names.

#![allow(dead_code)]

use pantometry_forcefield::{Histidine, Part, Residue, System};

/// Where each table came from, with its SHA-256 in full.
pub const SOURCES: [(&str, &str); 6] = [
    (
        "github.com/choderalab/ambermini f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
         share/amber/dat/leap/lib/amino12.lib",
        "dcd92afbb7ee4c8a737f1ad47eb714618b6c186759c902cd5c06c22a5378df77",
    ),
    (
        "github.com/choderalab/ambermini f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
         share/amber/dat/leap/lib/aminont12.lib",
        "e1ab30755b34cb6d0d38f9ff13c35b581c39f92a5cc071b307cbb1630897c2ee",
    ),
    (
        "github.com/choderalab/ambermini f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
         share/amber/dat/leap/lib/aminoct12.lib",
        "14c9bdf01190406089ebad30f9b97b8e2130b88d118c73e45685d13ed222a500",
    ),
    (
        "github.com/choderalab/ambermini f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
         share/amber/dat/leap/parm/parm10.dat",
        "1ca292d5679e35450ea111a8d1b9442b5ece8518e824a8e43ed290dfff1a2ae6",
    ),
    (
        "github.com/choderalab/ambermini f7c421c46dca9f0e0f87b10a9c87985b446d8ff9 \
         share/amber/dat/leap/parm/frcmod.ff14SB",
        "dd2bdcf9bb1f89598b997d3d6bfc8317db71cc80063136448469de17c4b63398",
    ),
    (
        "github.com/openmm/openmm 3a2326d0dd091229eaa4c0d7b5ad1001d5fe376e \
         wrappers/python/openmm/app/data/amber14/protein.ff14SB.xml (the cross-check)",
        "d9f9779c09d67cd5f8bc657692f174ffab14c469dfd06d560ac1899fa7e976b8",
    ),
];

/// Each template 181L needs: its AMBER name, and `atom=type` for every atom in the crate's names.
pub const TEMPLATES: [(&str, &str); 21] = [
    (
        "ALA",
        "N=N CA=CX C=C O=O CB=CT H=H HA=H1 HB1=HC HB2=HC HB3=HC",
    ),
    (
        "ARG",
        "N=N CA=CX C=C O=O CB=C8 CG=C8 CD=C8 NE=N2 CZ=CA NH1=N2 NH2=N2 H=H HA=H1 \
         HB2=HC HB3=HC HG2=HC HG3=HC HD2=H1 HD3=H1 HE=H HH11=H HH12=H HH21=H HH22=H",
    ),
    (
        "ASN",
        "N=N CA=CX C=C O=O CB=2C CG=C OD1=O ND2=N H=H HA=H1 HB2=HC HB3=HC HD21=H \
         HD22=H",
    ),
    (
        "ASP",
        "N=N CA=CX C=C O=O CB=2C CG=CO OD1=O2 OD2=O2 H=H HA=H1 HB2=HC HB3=HC",
    ),
    (
        "CLYS",
        "N=N CA=CX C=C O=O2 CB=C8 CG=C8 CD=C8 CE=C8 NZ=N3 OXT=O2 H=H HA=H1 HB2=HC \
         HB3=HC HG2=HC HG3=HC HD2=HC HD3=HC HE2=HP HE3=HP HZ1=H HZ2=H HZ3=H",
    ),
    (
        "GLN",
        "N=N CA=CX C=C O=O CB=2C CG=2C CD=C OE1=O NE2=N H=H HA=H1 HB2=HC HB3=HC \
         HG2=HC HG3=HC HE21=H HE22=H",
    ),
    (
        "GLU",
        "N=N CA=CX C=C O=O CB=2C CG=2C CD=CO OE1=O2 OE2=O2 H=H HA=H1 HB2=HC HB3=HC \
         HG2=HC HG3=HC",
    ),
    ("GLY", "N=N CA=CX C=C O=O H=H HA2=H1 HA3=H1"),
    (
        "HIP",
        "N=N CA=CX C=C O=O CB=CT CG=CC ND1=NA CD2=CW CE1=CR NE2=NA H=H HA=H1 HB2=HC \
         HB3=HC HD1=H HD2=H4 HE1=H5 HE2=H",
    ),
    (
        "ILE",
        "N=N CA=CX C=C O=O CB=3C CG1=2C CG2=CT CD1=CT H=H HA=H1 HB=HC HG12=HC HG13=HC \
         HG21=HC HG22=HC HG23=HC HD11=HC HD12=HC HD13=HC",
    ),
    (
        "LEU",
        "N=N CA=CX C=C O=O CB=2C CG=3C CD1=CT CD2=CT H=H HA=H1 HB2=HC HB3=HC HG=HC \
         HD11=HC HD12=HC HD13=HC HD21=HC HD22=HC HD23=HC",
    ),
    (
        "LYS",
        "N=N CA=CX C=C O=O CB=C8 CG=C8 CD=C8 CE=C8 NZ=N3 H=H HA=H1 HB2=HC HB3=HC \
         HG2=HC HG3=HC HD2=HC HD3=HC HE2=HP HE3=HP HZ1=H HZ2=H HZ3=H",
    ),
    (
        "MET",
        "N=N CA=CX C=C O=O CB=2C CG=2C SD=S CE=CT H=H HA=H1 HB2=HC HB3=HC HG2=H1 \
         HG3=H1 HE1=H1 HE2=H1 HE3=H1",
    ),
    (
        "NMET",
        "N=N3 CA=CX C=C O=O CB=2C CG=2C SD=S CE=CT H=H H2=H HA=HP HB2=HC HB3=HC \
         HG2=H1 HG3=H1 HE1=H1 HE2=H1 HE3=H1 H3=H",
    ),
    (
        "PHE",
        "N=N CA=CX C=C O=O CB=CT CG=CA CD1=CA CD2=CA CE1=CA CE2=CA CZ=CA H=H HA=H1 \
         HB2=HC HB3=HC HD1=HA HD2=HA HE1=HA HE2=HA HZ=HA",
    ),
    (
        "PRO",
        "N=N CA=CX C=C O=O CB=CT CG=CT CD=CT HA=H1 HB2=HC HB3=HC HG2=HC HG3=HC HD2=H1 \
         HD3=H1",
    ),
    (
        "SER",
        "N=N CA=CX C=C O=O CB=2C OG=OH H=H HA=H1 HB2=H1 HB3=H1 HG=HO",
    ),
    (
        "THR",
        "N=N CA=CX C=C O=O CB=3C OG1=OH CG2=CT H=H HA=H1 HB=H1 HG1=HO HG21=HC HG22=HC \
         HG23=HC",
    ),
    (
        "TRP",
        "N=N CA=CX C=C O=O CB=CT CG=C* CD1=CW CD2=CB NE1=NA CE2=CN CE3=CA CZ2=CA \
         CZ3=CA CH2=CA H=H HA=H1 HB2=HC HB3=HC HD1=H4 HE1=H HE3=HA HZ2=HA HZ3=HA \
         HH2=HA",
    ),
    (
        "TYR",
        "N=N CA=CX C=C O=O CB=CT CG=CA CD1=CA CD2=CA CE1=CA CE2=CA CZ=C OH=OH H=H \
         HA=H1 HB2=HC HB3=HC HD1=HA HD2=HA HE1=HA HE2=HA HH=HO",
    ),
    (
        "VAL",
        "N=N CA=CX C=C O=O CB=3C CG1=CT CG2=CT H=H HA=H1 HB=HC HG11=HC HG12=HC \
         HG13=HC HG21=HC HG22=HC HG23=HC",
    ),
];

/// `parm10.dat`'s equivalences, for the types here with no line of their own: `CA CB CC CN CR CW`
/// take `C*`'s, and `N2 NA` take `N`'s.
pub const EQUIVALENT: [(&str, &str); 8] = [
    ("CA", "C*"),
    ("CB", "C*"),
    ("CC", "C*"),
    ("CN", "C*"),
    ("CR", "C*"),
    ("CW", "C*"),
    ("N2", "N"),
    ("NA", "N"),
];

/// `R*/2` (ångström) and `ε` (kcal/mol) of every line the types here take: `parm10.dat`'s `MOD4 RE`
/// and, for `2C`, `3C`, `C8` and `CO`, `frcmod.ff14SB`'s `NONB`.
pub const LENNARD_JONES: [(&str, f64, f64); 22] = [
    ("H", 0.6000, 0.0157),
    ("HO", 0.0000, 0.0000),
    ("HC", 1.4870, 0.0157),
    ("H1", 1.3870, 0.0157),
    ("HP", 1.1000, 0.0157),
    ("HA", 1.4590, 0.0150),
    ("H4", 1.4090, 0.0150),
    ("H5", 1.3590, 0.0150),
    ("O", 1.6612, 0.2100),
    ("O2", 1.6612, 0.2100),
    ("OH", 1.7210, 0.2104),
    ("C*", 1.9080, 0.0860),
    ("CT", 1.9080, 0.1094),
    ("CX", 1.9080, 0.1094),
    ("C", 1.9080, 0.0860),
    ("N", 1.8240, 0.1700),
    ("N3", 1.8240, 0.1700),
    ("S", 2.0000, 0.2500),
    ("2C", 1.9080, 0.1094),
    ("3C", 1.9080, 0.1094),
    ("C8", 1.9080, 0.1094),
    ("CO", 1.9080, 0.0860),
];

/// The template a residue takes: `HIP`, `HID` or `HIE` for a histidine by its protonation, and `N`
/// or `C` before the name at a terminus. **Refused by name** for a residue the table does not have.
pub fn template_name(r: &Residue) -> String {
    let base = match (r.name.as_str(), r.histidine) {
        ("HIS", Some(Histidine::Both)) => "HIP".to_string(),
        ("HIS", Some(Histidine::Delta)) => "HID".to_string(),
        ("HIS", Some(Histidine::Epsilon)) => "HIE".to_string(),
        (name, _) => name.to_string(),
    };
    assert!(
        !(r.n_terminal && r.c_terminal),
        "{}: a residue that is both termini has no template here",
        r.label()
    );
    let name = if r.n_terminal {
        format!("N{base}")
    } else if r.c_terminal {
        format!("C{base}")
    } else {
        base
    };
    assert!(
        TEMPLATES.iter().any(|(t, _)| *t == name),
        "{}: no template {name} in the table",
        r.label()
    );
    name
}

/// The AMBER type of `system`'s atom `i`, from its residue's template. **Refused by name** for a
/// ligand's atom, or an atom its template does not have.
pub fn atom_type(system: &System, i: usize) -> &'static str {
    let r = system.residue_of(i);
    assert_eq!(r.part, Part::Protein, "atom {i} is not the protein's");
    let template = template_name(r);
    let atoms = TEMPLATES
        .iter()
        .find(|(t, _)| *t == template)
        .map(|(_, a)| *a)
        .expect("a template");
    let name = system.atom_name(i);
    atoms
        .split_whitespace()
        .find_map(|w| {
            let (a, t) = w.split_once('=').expect("atom=type");
            (a == name).then_some(t)
        })
        .unwrap_or_else(|| panic!("{}:{name}: not in {template}", r.label()))
}

/// The line type `t` takes its Lennard-Jones from: its own, or its equivalent's. **Refused by name**
/// for a type with neither.
pub fn lennard_jones_line(t: &str) -> &'static str {
    let t = EQUIVALENT
        .iter()
        .find(|(a, _)| *a == t)
        .map_or(t, |(_, b)| *b);
    LENNARD_JONES
        .iter()
        .find(|(a, _, _)| *a == t)
        .map(|(a, _, _)| *a)
        .unwrap_or_else(|| panic!("no Lennard-Jones line for {t}"))
}

/// Type `t`'s `R*/2` (ångström) and `ε` (kcal/mol).
pub fn lennard_jones(t: &str) -> (f64, f64) {
    let line = lennard_jones_line(t);
    let (_, r, e) = LENNARD_JONES
        .iter()
        .find(|(a, _, _)| *a == line)
        .expect("a line");
    (*r, *e)
}
