//! **FreeSolv's GAFF/AM1-BCC benzene**, for W4's charge perturbation: the `.mol2` its calculated
//! hydration free energy was run from, kept here byte for byte, and the library's reader and map
//! under the names this directory's tests use. Data, cited, not code.
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
//! - **The reader and the map are the library's** ([`pantometry_forcefield::mol2`]): [`parse`] and
//!   [`charges_onto`] are [`Mol2::parse`] and [`Mol2::onto_graph`] with their refusals as text.
//!   They were written here first and moved there whole.
//! - **And octan-1-ol** ([`MOBLEY_36119`]), for step 4c's solvent: the same archive at the same
//!   commit, kept the same way, its own SHA-256 in [`OCTANOL_SHA256`]. Its AM1-BCC charges on
//!   UFF's octanol are one of the models `tests/liquid_octanol_in_a_box.rs` measures.
//! - **And n-octane, n-nonane and n-decane** ([`MOBLEY_3167746`], [`MOBLEY_4375719`],
//!   [`MOBLEY_2197088`]), for step C-1's oil: the same archive and commit, each with its SHA-256.
//!   Decane's charges, by position along the chain, are n-hexadecane's in `tests/hexadecane/`, and
//!   octane's and nonane's are what that construction is held to. Decane's sum to −0.0005 e and
//!   octane's to +0.0001 as written; nonane's to zero.
//! - **And p-xylene** ([`MOBLEY_3234716`]), the second solute partitioned into octanol: the same
//!   archive and commit ([`P_XYLENE_SOURCE`]), its SHA-256 in [`P_XYLENE_SHA256`]. FreeSolv's row
//!   for it is −0.80 ± 0.60 kcal/mol by experiment (the ± a default, FreeSolv says, not a measured
//!   uncertainty) and −0.66 ± 0.03 calculated. Its charges sum to +0.0001 e as written.

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

pub use pantometry_forcefield::mol2::{Mapped, Mol2};

/// An atom name's element ([`pantometry_forcefield::mol2::element_of`]), its refusal as text.
pub fn element_of(name: &str) -> Result<Element, String> {
    pantometry_forcefield::mol2::element_of(name).map_err(|e| e.to_string())
}

/// **One molecule of a `.mol2`, read strictly** by the library's [`Mol2::parse`], its refusal as
/// text so that a test can ask what it names.
pub fn parse(text: &str) -> Result<Mol2, String> {
    Mol2::parse(text).map_err(|e| e.to_string())
}

/// **`mol2`'s charges on the molecule whose atoms are `elements` and bonds `bonds`**, by the
/// library's [`Mol2::onto_graph`]; refused, as text, when they are not the same molecule.
pub fn charges_onto(
    mol2: &Mol2,
    elements: &[Element],
    bonds: &[[usize; 2]],
) -> Result<Mapped, String> {
    mol2.onto_graph(elements, bonds).map_err(|e| e.to_string())
}

/// A component's elements and bonds, as [`charges_onto`] takes them.
pub fn component_graph(c: &Component) -> (Vec<Element>, Vec<[usize; 2]>) {
    pantometry_forcefield::mol2::graph(c)
}

/// FreeSolv's benzene, read ([`MOBLEY_3053621`]).
pub fn benzene() -> Mol2 {
    parse(MOBLEY_3053621).expect("FreeSolv's mobley_3053621.mol2")
}

/// Where [`MOBLEY_36119`] came from: the same archive, at the same commit, as [`SOURCE`].
pub const OCTANOL_SOURCE: &str = "github.com/MobleyLab/FreeSolv v0.52, commit \
                                  9d8a8097d587a6c084befbb2abde6dd2a0f31ba3, mol2files_gaff.tar.gz \
                                  (SHA-256 15ece6114442a8eb5e720c0ad25bb48e4a6a35392da36f72b80a6f2f415fd133) \
                                  → mol2files_gaff/mobley_36119.mol2";

/// The SHA-256 of [`MOBLEY_36119`]'s bytes.
pub const OCTANOL_SHA256: &str = "ecfa179de0b94f48c38072f66bae1932fd0ac0e297b11969dafed89468929831";

/// **FreeSolv v0.52's `mobley_36119.mol2`, octan-1-ol, byte for byte** ([`OCTANOL_SOURCE`]),
/// CC BY 4.0 as [`MOBLEY_3053621`] is: GAFF's types (`c3`, `oh`, `hc`, `h1`, `ho`) and AM1-BCC's
/// charges to four places, as FreeSolv ran its hydration free energy. Its database row,
/// `mobley_36119; CCCCCCCCO; octan-1-ol`, is [`OCTANOL_EXPERIMENT`] and [`OCTANOL_CALCULATED`].
/// The charges sum to zero exactly, as written: O1 −0.6023, its H18 +0.3984, C8 +0.1314, the
/// other seven carbons −0.0918 to −0.0782 and the seventeen hydrogens on carbon +0.0194 to +0.0552.
pub const MOBLEY_36119: &str = r#"@<TRIPOS>MOLECULE
octan-1-ol
   27    26     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          2.5800    3.1400   -2.3930 c3        1 MOL     -0.091800
      2 C2          2.1180    1.7400   -2.7630 c3        1 MOL     -0.080400
      3 C3          0.9800    1.7770   -3.7840 c3        1 MOL     -0.078800
      4 C4          0.5210    0.3650   -4.1520 c3        1 MOL     -0.078600
      5 C5         -0.6200    0.3420   -5.1750 c3        1 MOL     -0.078300
      6 C6         -1.0790   -1.0710   -5.5420 c3        1 MOL     -0.078200
      7 C7         -2.2160   -1.0340   -6.5630 c3        1 MOL     -0.083300
      8 C8         -2.6810   -2.4350   -6.9360 c3        1 MOL      0.131400
      9 O1         -3.7340   -2.3320   -7.8840 oh        1 MOL     -0.602300
     10 H1          3.3940    3.1040   -1.6610 hc        1 MOL      0.032900
     11 H2          1.7580    3.7240   -1.9640 hc        1 MOL      0.032900
     12 H3          2.9440    3.6670   -3.2800 hc        1 MOL      0.032900
     13 H4          1.7880    1.2190   -1.8570 hc        1 MOL      0.037500
     14 H5          2.9670    1.1780   -3.1700 hc        1 MOL      0.037500
     15 H6          0.1440    2.3480   -3.3670 hc        1 MOL      0.039700
     16 H7          1.3150    2.3020   -4.6870 hc        1 MOL      0.039700
     17 H8          0.2130   -0.1730   -3.2470 hc        1 MOL      0.038900
     18 H9          1.3690   -0.1890   -4.5750 hc        1 MOL      0.038900
     19 H10        -1.4690    0.9200   -4.7910 hc        1 MOL      0.041600
     20 H11        -0.2740    0.8530   -6.0820 hc        1 MOL      0.041600
     21 H12        -1.4210   -1.6000   -4.6460 hc        1 MOL      0.039300
     22 H13        -0.2350   -1.6370   -5.9540 hc        1 MOL      0.039300
     23 H14        -1.8950   -0.5060   -7.4690 hc        1 MOL      0.055200
     24 H15        -3.0640   -0.4630   -6.1660 hc        1 MOL      0.055200
     25 H16        -1.8640   -3.0110   -7.3840 h1        1 MOL      0.019400
     26 H17        -3.0540   -2.9700   -6.0570 h1        1 MOL      0.019400
     27 H18        -4.0060   -3.2390   -8.1000 ho        1 MOL      0.398400
@<TRIPOS>BOND
     1    1    2 1   
     2    2    3 1   
     3    3    4 1   
     4    4    5 1   
     5    5    6 1   
     6    6    7 1   
     7    7    8 1   
     8    8    9 1   
     9    1   10 1   
    10    1   11 1   
    11    1   12 1   
    12    2   13 1   
    13    2   14 1   
    14    3   15 1   
    15    3   16 1   
    16    4   17 1   
    17    4   18 1   
    18    5   19 1   
    19    5   20 1   
    20    6   21 1   
    21    6   22 1   
    22    7   23 1   
    23    7   24 1   
    24    8   25 1   
    25    8   26 1   
    26    9   27 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's experimental hydration free energy of octan-1-ol and its uncertainty, kcal/mol: the
/// database row of [`OCTANOL_SOURCE`]'s commit, whose uncertainty is the default FreeSolv assigns
/// when none was available.
pub const OCTANOL_EXPERIMENT: (f64, f64) = (-4.09, 0.60);

/// FreeSolv's calculated GAFF/AM1-BCC hydration free energy of octan-1-ol in TIP3P and its
/// uncertainty, kcal/mol: the same row.
pub const OCTANOL_CALCULATED: (f64, f64) = (-2.69, 0.03);

/// FreeSolv's octan-1-ol, read ([`MOBLEY_36119`]).
pub fn octanol() -> Mol2 {
    parse(MOBLEY_36119).expect("FreeSolv's mobley_36119.mol2")
}

// ---------------------------------------------------------------------------------------------
// Three n-alkanes, for step C-1's hexadecane

/// Where the three n-alkanes' files came from: the same archive, at the same commit, as
/// [`SOURCE`], each `mol2files_gaff/<id>.mol2`, CC BY 4.0 as the rest. Used by
/// `tests/hexadecane/`, which builds n-hexadecane's charges from decane's by position and holds the
/// construction to octane's, nonane's and decane's own.
pub const ALKANE_SOURCE: &str = "github.com/MobleyLab/FreeSolv v0.52, commit \
    9d8a8097d587a6c084befbb2abde6dd2a0f31ba3, mol2files_gaff.tar.gz \
    (SHA-256 15ece6114442a8eb5e720c0ad25bb48e4a6a35392da36f72b80a6f2f415fd133) \
    → mol2files_gaff/mobley_3167746.mol2 (octane), mobley_4375719.mol2 \
    (nonane), mobley_2197088.mol2 (decane)";

/// The SHA-256 of [`MOBLEY_3167746`]'s bytes, octane.
pub const OCTANE_SHA256: &str = "f1b36a3b6b2e588725c43764e00421531c0633f3608f6ef57500fe2dd9846f9a";

/// **FreeSolv v0.52's `mobley_3167746.mol2`, n-octane (C8), byte for byte** ([`ALKANE_SOURCE`]):
/// GAFF's `c3` and `hc` and AM1-BCC's charges to four places, as FreeSolv ran them.
pub const MOBLEY_3167746: &str = r#"@<TRIPOS>MOLECULE
octane
   26    25     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          0.5810   -1.0130    0.8660 c3        1 MOL     -0.092200
      2 C2          0.2770   -0.4860    2.2600 c3        1 MOL     -0.080000
      3 C3          0.8900    0.8980    2.4780 c3        1 MOL     -0.079300
      4 C4          0.6120    1.4740    3.8700 c3        1 MOL     -0.078400
      5 C5          1.2250    2.8590    4.0890 c3        1 MOL     -0.078400
      6 C6          0.9170    3.3850    5.4930 c3        1 MOL     -0.079300
      7 C7          1.5280    4.7680    5.7150 c3        1 MOL     -0.080000
      8 C8          1.2240    5.2950    7.1080 c3        1 MOL     -0.092100
      9 H1          0.1400   -2.0040    0.7160 hc        1 MOL      0.032500
     10 H2          0.1860   -0.3410    0.0970 hc        1 MOL      0.032500
     11 H3          1.6620   -1.0970    0.7200 hc        1 MOL      0.032500
     12 H4         -0.8100   -0.4470    2.3890 hc        1 MOL      0.037900
     13 H5          0.6690   -1.1900    3.0030 hc        1 MOL      0.037900
     14 H6          1.9770    0.8320    2.3470 hc        1 MOL      0.038900
     15 H7          0.5210    1.5910    1.7130 hc        1 MOL      0.038900
     16 H8         -0.4710    1.5160    4.0410 hc        1 MOL      0.039400
     17 H9          1.0220    0.7800    4.6140 hc        1 MOL      0.039400
     18 H10         0.8330    3.5650    3.3480 hc        1 MOL      0.039400
     19 H11         2.3110    2.8070    3.9470 hc        1 MOL      0.039400
     20 H12         1.3090    2.6860    6.2420 hc        1 MOL      0.038900
     21 H13        -0.1700    3.4350    5.6360 hc        1 MOL      0.038900
     22 H14         1.1360    5.4720    4.9710 hc        1 MOL      0.037900
     23 H15         2.6140    4.7220    5.5770 hc        1 MOL      0.037900
     24 H16         1.6700    6.2850    7.2450 hc        1 MOL      0.032500
     25 H17         1.6310    4.6290    7.8760 hc        1 MOL      0.032500
     26 H18         0.1440    5.3830    7.2670 hc        1 MOL      0.032500
@<TRIPOS>BOND
     1    1    2 1   
     2    2    3 1   
     3    3    4 1   
     4    4    5 1   
     5    5    6 1   
     6    6    7 1   
     7    7    8 1   
     8    1    9 1   
     9    1   10 1   
    10    1   11 1   
    11    2   12 1   
    12    2   13 1   
    13    3   14 1   
    14    3   15 1   
    15    4   16 1   
    16    4   17 1   
    17    5   18 1   
    18    5   19 1   
    19    6   20 1   
    20    6   21 1   
    21    7   22 1   
    22    7   23 1   
    23    8   24 1   
    24    8   25 1   
    25    8   26 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's n-octane, read ([`MOBLEY_3167746`]).
pub fn octane() -> Mol2 {
    parse(MOBLEY_3167746).expect("FreeSolv's mobley_3167746.mol2")
}

/// The SHA-256 of [`MOBLEY_4375719`]'s bytes, nonane.
pub const NONANE_SHA256: &str = "04776cff464af14f61cb1b1dd2b07fc0b8c86904e2891a91de36f66e3b4993af";

/// **FreeSolv v0.52's `mobley_4375719.mol2`, n-nonane (C9), byte for byte** ([`ALKANE_SOURCE`]):
/// GAFF's `c3` and `hc` and AM1-BCC's charges to four places, as FreeSolv ran them.
pub const MOBLEY_4375719: &str = r#"@<TRIPOS>MOLECULE
nonane
   29    28     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          0.5000   -2.0150   -1.2890 c3        1 MOL     -0.092400
      2 C2          0.8190   -1.9280    0.1950 c3        1 MOL     -0.079800
      3 C3          2.0160   -2.7840    0.6140 c3        1 MOL     -0.079500
      4 C4          2.3340   -2.6940    2.1080 c3        1 MOL     -0.078300
      5 C5          3.5370   -3.5670    2.4730 c3        1 MOL     -0.078600
      6 C6          3.8560   -3.4790    3.9660 c3        1 MOL     -0.078300
      7 C7          5.0580   -4.3510    4.3310 c3        1 MOL     -0.079500
      8 C8          5.3770   -4.2630    5.8230 c3        1 MOL     -0.079800
      9 C9          6.5710   -5.1280    6.1910 c3        1 MOL     -0.092400
     10 H1         -0.3610   -1.3890   -1.5440 hc        1 MOL      0.032500
     11 H2          1.3560   -1.6690   -1.8780 hc        1 MOL      0.032500
     12 H3          0.2840   -3.0430   -1.5940 hc        1 MOL      0.032500
     13 H4          0.9930   -0.8820    0.4700 hc        1 MOL      0.038300
     14 H5         -0.0600   -2.2620    0.7600 hc        1 MOL      0.038300
     15 H6          1.7980   -3.8280    0.3560 hc        1 MOL      0.038900
     16 H7          2.8980   -2.4930    0.0300 hc        1 MOL      0.038900
     17 H8          1.4600   -3.0120    2.6900 hc        1 MOL      0.039200
     18 H9          2.5480   -1.6560    2.3850 hc        1 MOL      0.039200
     19 H10         3.3300   -4.6100    2.2030 hc        1 MOL      0.039000
     20 H11         4.4100   -3.2490    1.8910 hc        1 MOL      0.039000
     21 H12         4.0630   -2.4360    4.2370 hc        1 MOL      0.039200
     22 H13         2.9820   -3.7960    4.5480 hc        1 MOL      0.039200
     23 H14         5.9330   -4.0340    3.7490 hc        1 MOL      0.038900
     24 H15         4.8510   -5.3940    4.0600 hc        1 MOL      0.038900
     25 H16         5.5890   -3.2230    6.0970 hc        1 MOL      0.038300
     26 H17         4.5070   -4.5820    6.4090 hc        1 MOL      0.038300
     27 H18         6.7810   -5.0490    7.2620 hc        1 MOL      0.032500
     28 H19         7.4670   -4.8130    5.6450 hc        1 MOL      0.032500
     29 H20         6.3800   -6.1810    5.9580 hc        1 MOL      0.032500
@<TRIPOS>BOND
     1    1    2 1   
     2    2    3 1   
     3    3    4 1   
     4    4    5 1   
     5    5    6 1   
     6    6    7 1   
     7    7    8 1   
     8    8    9 1   
     9    1   10 1   
    10    1   11 1   
    11    1   12 1   
    12    2   13 1   
    13    2   14 1   
    14    3   15 1   
    15    3   16 1   
    16    4   17 1   
    17    4   18 1   
    18    5   19 1   
    19    5   20 1   
    20    6   21 1   
    21    6   22 1   
    22    7   23 1   
    23    7   24 1   
    24    8   25 1   
    25    8   26 1   
    26    9   27 1   
    27    9   28 1   
    28    9   29 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's n-nonane, read ([`MOBLEY_4375719`]).
pub fn nonane() -> Mol2 {
    parse(MOBLEY_4375719).expect("FreeSolv's mobley_4375719.mol2")
}

/// The SHA-256 of [`MOBLEY_2197088`]'s bytes, decane.
pub const DECANE_SHA256: &str = "3bb431594417968152b3a088822c2dac5f699990b1fb5029e5cc1884450a301a";

/// **FreeSolv v0.52's `mobley_2197088.mol2`, n-decane (C10), byte for byte** ([`ALKANE_SOURCE`]):
/// GAFF's `c3` and `hc` and AM1-BCC's charges to four places, as FreeSolv ran them.
pub const MOBLEY_2197088: &str = r#"@<TRIPOS>MOLECULE
decane
   32    31     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          0.4670   -0.7510    0.2710 c3        1 MOL     -0.092000
      2 C2         -0.2980   -2.0320    0.5650 c3        1 MOL     -0.079100
      3 C3         -1.2680   -2.3760   -0.5660 c3        1 MOL     -0.078300
      4 C4         -2.0370   -3.6650   -0.2670 c3        1 MOL     -0.077900
      5 C5         -3.0070   -4.0070   -1.4000 c3        1 MOL     -0.077700
      6 C6         -3.8040   -5.2920   -1.1510 c3        1 MOL     -0.077700
      7 C7         -4.7730   -5.6310   -2.2860 c3        1 MOL     -0.077900
      8 C8         -5.5700   -6.9160   -2.0370 c3        1 MOL     -0.078300
      9 C9         -6.5390   -7.2580   -3.1690 c3        1 MOL     -0.079100
     10 C10        -7.3050   -8.5390   -2.8750 c3        1 MOL     -0.091900
     11 H1          1.1550   -0.5240    1.0920 hc        1 MOL      0.032700
     12 H2         -0.2160    0.0970    0.1560 hc        1 MOL      0.032700
     13 H3          1.0550   -0.8470   -0.6470 hc        1 MOL      0.032700
     14 H4         -0.8490   -1.9140    1.5050 hc        1 MOL      0.038000
     15 H5          0.4150   -2.8530    0.7060 hc        1 MOL      0.038000
     16 H6         -0.7120   -2.4880   -1.5050 hc        1 MOL      0.038400
     17 H7         -1.9760   -1.5500   -0.7060 hc        1 MOL      0.038400
     18 H8         -2.5850   -3.5450    0.6740 hc        1 MOL      0.038400
     19 H9         -1.3280   -4.4900   -0.1280 hc        1 MOL      0.038400
     20 H10        -2.4410   -4.1310   -2.3310 hc        1 MOL      0.038500
     21 H11        -3.6980   -3.1700   -1.5610 hc        1 MOL      0.038500
     22 H12        -4.3470   -5.1970   -0.2040 hc        1 MOL      0.038500
     23 H13        -3.0910   -6.1170   -1.0280 hc        1 MOL      0.038500
     24 H14        -5.4720   -4.8030   -2.4470 hc        1 MOL      0.038400
     25 H15        -4.2070   -5.7560   -3.2170 hc        1 MOL      0.038400
     26 H16        -4.8570   -7.7410   -1.9140 hc        1 MOL      0.038400
     27 H17        -6.1190   -6.8310   -1.0910 hc        1 MOL      0.038400
     28 H18        -7.2590   -6.4450   -3.3120 hc        1 MOL      0.038000
     29 H19        -5.9880   -7.3760   -4.1090 hc        1 MOL      0.038000
     30 H20        -7.9930   -8.7660   -3.6960 hc        1 MOL      0.032700
     31 H21        -6.6220   -9.3870   -2.7610 hc        1 MOL      0.032700
     32 H22        -7.8930   -8.4430   -1.9570 hc        1 MOL      0.032700
@<TRIPOS>BOND
     1    1    2 1   
     2    2    3 1   
     3    3    4 1   
     4    4    5 1   
     5    5    6 1   
     6    6    7 1   
     7    7    8 1   
     8    8    9 1   
     9    9   10 1   
    10    1   11 1   
    11    1   12 1   
    12    1   13 1   
    13    2   14 1   
    14    2   15 1   
    15    3   16 1   
    16    3   17 1   
    17    4   18 1   
    18    4   19 1   
    19    5   20 1   
    20    5   21 1   
    21    6   22 1   
    22    6   23 1   
    23    7   24 1   
    24    7   25 1   
    25    8   26 1   
    26    8   27 1   
    27    9   28 1   
    28    9   29 1   
    29   10   30 1   
    30   10   31 1   
    31   10   32 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's n-decane, read ([`MOBLEY_2197088`]).
pub fn decane() -> Mol2 {
    parse(MOBLEY_2197088).expect("FreeSolv's mobley_2197088.mol2")
}

/// Where [`MOBLEY_3234716`] came from: the archive [`SOURCE`] names, at the same commit.
pub const P_XYLENE_SOURCE: &str = "github.com/MobleyLab/FreeSolv v0.52, commit \
    9d8a8097d587a6c084befbb2abde6dd2a0f31ba3, mol2files_gaff.tar.gz \
    (SHA-256 15ece6114442a8eb5e720c0ad25bb48e4a6a35392da36f72b80a6f2f415fd133) \
    → mol2files_gaff/mobley_3234716.mol2";

/// The SHA-256 of [`MOBLEY_3234716`]'s bytes, p-xylene, reproduced on 2026-11-02.
pub const P_XYLENE_SHA256: &str =
    "f87c244e2259defbaccefdfeabaa059dcb70d20487ad08aad39f7089508aad02";

/// **FreeSolv v0.52's `mobley_3234716.mol2`, p-xylene, byte for byte** ([`P_XYLENE_SOURCE`]):
/// GAFF's `c3`, `ca`, `hc` and `ha`, and AM1-BCC's charges to four places, as FreeSolv ran them.
pub const MOBLEY_3234716: &str = r#"@<TRIPOS>MOLECULE
p-xylene
   18    18     1     0     0
SMALL
No Charge or Current Charge


@<TRIPOS>ATOM
      1 C1          1.9310   -1.0130   -1.6030 c3        1 MOL     -0.052900
      2 C2          1.8850   -1.0360   -0.1120 ca        1 MOL     -0.081700
      3 C3          2.9210   -1.6310    0.6080 ca        1 MOL     -0.127300
      4 C4          2.8780   -1.6520    2.0020 ca        1 MOL     -0.127300
      5 C5          1.8000   -1.0790    2.6760 ca        1 MOL     -0.081700
      6 C6          0.7640   -0.4840    1.9550 ca        1 MOL     -0.127300
      7 C7          0.8070   -0.4630    0.5610 ca        1 MOL     -0.127300
      8 C8          1.7530   -1.1040    4.1660 c3        1 MOL     -0.052800
      9 H1          2.4380   -0.1080   -1.9480 hc        1 MOL      0.043200
     10 H2          2.4630   -1.8880   -1.9910 hc        1 MOL      0.043200
     11 H3          0.9200   -1.0360   -2.0250 hc        1 MOL      0.043200
     12 H4          3.7650   -2.0800    0.0920 ha        1 MOL      0.129800
     13 H5          3.6890   -2.1180    2.5550 ha        1 MOL      0.129800
     14 H6         -0.0810   -0.0360    2.4720 ha        1 MOL      0.129800
     15 H7         -0.0050    0.0020    0.0080 ha        1 MOL      0.129800
     16 H8          2.2380   -0.2110    4.5790 hc        1 MOL      0.043200
     17 H9          0.7190   -1.1380    4.5320 hc        1 MOL      0.043200
     18 H10         2.2600   -1.9910    4.5660 hc        1 MOL      0.043200
@<TRIPOS>BOND
     1    1    2 1   
     2    2    7 ar  
     3    2    3 ar  
     4    3    4 ar  
     5    4    5 ar  
     6    5    6 ar  
     7    6    7 ar  
     8    5    8 1   
     9    1    9 1   
    10    1   10 1   
    11    1   11 1   
    12    3   12 1   
    13    4   13 1   
    14    6   14 1   
    15    7   15 1   
    16    8   16 1   
    17    8   17 1   
    18    8   18 1   
@<TRIPOS>SUBSTRUCTURE
     1 MOL         1 TEMP              0 ****  ****    0 ROOT
"#;

/// FreeSolv's p-xylene, read ([`MOBLEY_3234716`]).
pub fn p_xylene() -> Mol2 {
    parse(MOBLEY_3234716).expect("FreeSolv's mobley_3234716.mol2")
}
