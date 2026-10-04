# Chemical components, unmodified

Each file is one entry of the wwPDB Chemical Component Dictionary as RCSB serves it, byte for byte,
from `https://files.rcsb.org/ligands/download/<ID>.cif`. Nothing is trimmed, reordered or
regenerated: a fixture that has been through a script can no longer be checked against the thing
it came from. The tests that need a malformed file make one from this text in memory, one edit at a
time, and assert the edit's anchor occurs exactly once before making it.

| file | component (`_chem_comp.name`) | formula | modified | fetched | SHA-256 as fetched |
| --- | --- | --- | --- | --- | --- |
| `AIN.cif` | 2-(acetyloxy)benzoic acid — aspirin | C9 H8 O4 | 2020 | 2026-10-01 | `e24ecc21a59aa8e26e1143234f92f3e1e309a589508f38d4d91e078a456c57ed` |
| `NME.cif` | METHYLAMINE | C H5 N | 2024-09-27 | 2026-10-02 | `fb95b4cd28ade57c1c72e41fbdbdc243bb5770d4f8d45586f6a1c52d21767447` |
| `MOH.cif` | METHANOL | C H4 O | 2024-09-27 | 2026-10-02 | `2217d45623bd86341883043ae4e58a42bb0972376362a70b58ef4deb4812377b` |
| `MEE.cif` | METHANETHIOL | C H4 S | 2024-09-27 | 2026-10-02 | `6a38f431c53910d37bb89764b2d8b6611c66a2a67e77286dcdefbbc9d295d848` |
| `PEO.cif` | HYDROGEN PEROXIDE | H2 O2 | 2024-09-27 | 2026-10-02 | `12a6699248a0106e51862e1e1902eaa1803aa4b13e69416b14de9c7ddca7fd8a` |
| `S2H.cif` | Hydrogen disulfide | H2 S2 | 2024-09-27 | 2026-10-02 | `b4d89cf21e65479b2652ef73594fed24b31a1f0207af0fcadca2bd58acc4f0d7` |
| `A1JFW.cif` | methoxybenzene — anisole | C7 H8 O | 2026-03-27 | 2026-10-02 | `523fc55503e31e737542eee690c138315aadc702068c618db89d5d3171f8814b` |
| `16R.cif` | (methylsulfanyl)benzene — thioanisole | C7 H8 S | 2021-03-01 | 2026-10-02 | `dca705dcac41e6d9ef45ea4f15dda75f39a32526e31ee31fecde7fa46a1c8c34` |
| `ACE.cif` | ACETYL GROUP — used as acetaldehyde, see below | C2 H4 O | 2025-12-15 | 2026-10-02 | `463cfab0e3ed37e0c6308546855d30c14408b2031fe9c0996acd67185c314074` |
| `61G.cif` | 2-methylbuta-1,3-diene — isoprene | C5 H8 | 2016-02-26 | 2026-10-02 | `d6a1bb918716111d72e257b3a8582e187d9703ce8ba6f00fc73e6666f5952189` |
| `PYJ.cif` | PHENYLETHANE — ethylbenzene | C8 H10 | 2024-09-27 | 2026-10-02 | `f480ba3775822cb4fbf88d08009e9a43216807b4ec3b0145d8b31e07ad922c11` |
| `2F2.cif` | dimethyl ether | C2 H6 O | 2014-07-11 | 2026-10-02 | `7e77c3dd3780fba0b4e24699d6559cb594f32a299e9f7d0ea38dbed0daf27c4c` |
| `2ME.cif` | METHOXYETHANE — methyl ethyl ether | C3 H8 O | 2011-06-04 | 2026-10-02 | `8b83702be5724865af6d66199c5c89aebdb7a4c4d18a84c6bb846a26d8a1819f` |
| `ACN.cif` | ACETONE | C3 H6 O | 2024-09-27 | 2026-10-02 | `82f12f2c170ba00442dc0582080d232b3160b7d2fb03dae1cecdd5bf2c2c958d` |
| `DMN.cif` | DIMETHYLAMINE | C2 H7 N | 2026-03-27 | 2026-10-02 | `80bf432f4b74e2351b1f6567f2aa9d2f82b47918ca927bf4211f954dc3ee085c` |
| `KEN.cif` | N,N-dimethylmethanamine — trimethylamine | C3 H9 N | 2011-06-04 | 2026-10-02 | `b128f87740c25ee8c9a25af8c19b997af26df787982f9cf14509ba7bf9551e1d` |
| `CCN.cif` | ACETONITRILE | C2 H3 N | 2011-06-04 | 2026-10-02 | `ba528d6b625f29c57ebc805afec6fbab7f9338d4866d4b88f43fcd5e63f48639` |
| `ACM.cif` | ACETAMIDE | C2 H5 N O | 2024-09-27 | 2026-10-02 | `c4f6f007bbd5ab7840a27238002a6f149312da5791758c422322fdb100e9d2c3` |
| `TME.cif` | PROPANE | C3 H8 | 2024-09-27 | 2026-10-02 | `456caaa8dcced177e773d58f571980790d6e19b7e47d6b4a4dd96dcce107ad50` |
| `6AC.cif` | prop-2-enenitrile — acrylonitrile | C3 H3 N | 2021-03-01 | 2026-10-02 | `88e9b633186164c4a93712a67e9813dfd1bd21d0ddb89ba7f08e5cd3a3217622` |
| `DMF.cif` | DIMETHYLFORMAMIDE — N,N-dimethylformamide | C3 H7 N O | 2011-06-04 | 2026-10-02 | `e2ccace4fe9e14d8285e71cd804498ff5fa763911abcd55343bc47b232e6eb49` |
| `ALA.cif` | ALANINE | C3 H7 N O2 | 2024-09-27 | 2026-10-04 | `6d32b34d4f7b3ddf0cd3dff3f98ddaf7649bc5303ff9a8bd95ba62283f47a1ca` |
| `ARG.cif` | ARGININE | C6 H15 N4 O2 | 2024-09-27 | 2026-10-04 | `191c437cfa501e0a87d275ca3d39974a1a9c8dc704df61bfc442a4c314e196ac` |
| `ASN.cif` | ASPARAGINE | C4 H8 N2 O3 | 2024-09-27 | 2026-10-04 | `81491bbdb49f96e2945192c07510324d206298d3479029b76defba2ed6ab2dc7` |
| `ASP.cif` | ASPARTIC ACID | C4 H7 N O4 | 2024-09-27 | 2026-10-04 | `82c4a2d0b5eda0a32c36f106f0bb00edbc5bc51172db414fdc494025b7b8c858` |
| `CYS.cif` | CYSTEINE | C3 H7 N O2 S | 2024-09-27 | 2026-10-04 | `b28c79381f9f41fe57a9e97468fc8101d5fd1064272a22da1eed6d05d27a24cc` |
| `GLN.cif` | GLUTAMINE | C5 H10 N2 O3 | 2024-09-27 | 2026-10-04 | `b542ea5c013fbe0c7acee8e3f42d5dc4e392bd16c0aa545edbf6cb3c335ca30f` |
| `GLU.cif` | GLUTAMIC ACID | C5 H9 N O4 | 2024-09-27 | 2026-10-04 | `99fa19054d31b07e7101b33ad9cc847cc6a209da60de279aed5e6a9236bb6da4` |
| `GLY.cif` | GLYCINE | C2 H5 N O2 | 2024-09-27 | 2026-10-04 | `c49458946b0ebc057db6ad0a4e1557a1caaed4c80a203accd458efddccbf92ff` |
| `HIS.cif` | HISTIDINE | C6 H10 N3 O2 | 2023-11-03 | 2026-10-04 | `89355e09d9c622dd514738e559d7a75ab0a1659df9a179978e88e1a44704a46f` |
| `ILE.cif` | ISOLEUCINE | C6 H13 N O2 | 2023-11-03 | 2026-10-04 | `8166adcc6acb4661f4a167b0a4295879612dfad3f1a6b6f1b7bd0589c5bded93` |
| `LEU.cif` | LEUCINE | C6 H13 N O2 | 2024-09-27 | 2026-10-04 | `fc56e022bc8baa25dc6e543b2c46664e365779aa1dda9831576fdf10b30ab84b` |
| `LYS.cif` | LYSINE | C6 H15 N2 O2 | 2024-09-27 | 2026-10-04 | `318a0f6fa8f98682b35e1793517d5b683b9045c1e6174bc564b7d43e07e4c07c` |
| `MET.cif` | METHIONINE | C5 H11 N O2 S | 2024-09-27 | 2026-10-04 | `a89adea8738d069b0a7bb0c3efc26b362dcc624ef9abd3b3f4bfaedd785f3afd` |
| `PHE.cif` | PHENYLALANINE | C9 H11 N O2 | 2024-09-27 | 2026-10-04 | `002c626d24aec78fd622459d4f4470c682c7c3f9f97d32414fec4313d780496e` |
| `PRO.cif` | PROLINE | C5 H9 N O2 | 2024-09-27 | 2026-10-04 | `08674e4d60bc40ad5702b8de6e5e4abb66c2bbb603832c021c31194f4b4e593c` |
| `SER.cif` | SERINE | C3 H7 N O3 | 2024-09-27 | 2026-10-04 | `46db0296a07732854d2c3ca67ed0f74f4c82ddc5698d0e59f433039a5d7c88a0` |
| `THR.cif` | THREONINE | C4 H9 N O3 | 2023-11-03 | 2026-10-04 | `feb75dfac5bfcc7045de7a75a31575f3bf3e5b411716a06e92712924759dd554` |
| `TRP.cif` | TRYPTOPHAN | C11 H12 N2 O2 | 2023-11-03 | 2026-10-04 | `4cea7e281ee404f3e72fc20d3ba98064373cdd9309a7f5c75479a04c4b6e36d6` |
| `TYR.cif` | TYROSINE | C9 H11 N O3 | 2024-09-27 | 2026-10-04 | `fa0ef32f9d623e087a527b0a146c0494a5a656037c6448d956072dd2406fa554` |
| `VAL.cif` | VALINE | C5 H11 N O2 | 2023-11-03 | 2026-10-04 | `9c3c807970c817b1e1a4b43164db91a5bf1e7bf438b3269c5649f6ce84d4c0ed` |
| `BNZ.cif` | BENZENE | C6 H6 | 2011-06-04 | 2026-10-04 | `01bcf7c3ce9befdb4078e9832252eb5fe99e2598f320f87358ea1a9a247f7c61` |

**How each entry was found and checked.** By exact name or by formula through RCSB's search API,
then read: its name, its formula, and its bond table, which must be the molecule's — every heavy
atom's bonds and every hydrogen's one. The test `every_component_is_the_molecule_it_says` holds
each to its own stated formula and one bond per hydrogen. **`ACE` is named "ACETYL GROUP"**: it is
the dictionary's N-terminal acetyl cap, and as an entry it is CH₃–CH=O with all four hydrogens — C=O,
C–CH₃, the aldehyde C–H and three methyl C–H — which is acetaldehyde atom for atom; it is used as
acetaldehyde for that reason, and the name is the dictionary's, not this crate's.

Ideal coordinates are what the crate uses. Aspirin's model coordinates come from PDB entry `1OXR`;
the fallback test uses them, and so does the comparison of aspirin's relaxed structure with them.

**Not in the dictionary, and so hand-built** in `../tests/hand_built/` — each file says so in its
first line: ethane, methylphosphine, N-methylformamide, methyl formate, methyl vinyl ether,
*trans*-dimethyldiazene, propene, *s-trans*-butadiene and propyne. Their coordinates are a
z-matrix's, and every use minimises them before reading any geometry, so only their bonds and their
conformer carry information.

**The twenty amino acids and benzene** are the residue templates of `crate::pdb`: their bond
orders, aromatic flags, formal charges and hydrogens are what a crystal structure lacks. Each is
the free amino acid — NH₂ and COOH, Asp and Glu neutral, Lys, Arg and His charged (`_chem_comp.
pdbx_formal_charge` 1) — and the module's protonation table says what is removed and recharged.
`every_component_is_the_molecule_it_says` holds each to its own formula, and
`the_dictionarys_templates_are_free_amino_acids` to that protonation.

# Structures, unmodified

| file | entry | fetched from | revision | fetched | SHA-256 as fetched |
| --- | --- | --- | --- | --- | --- |
| `181L.pdb` | T4 lysozyme L99A (pseudo-wild-type C54T, C97A) with benzene, X-ray 1.80 Å | `https://files.rcsb.org/download/181L.pdb` | REVDAT 5, 07-FEB-24 | 2026-10-04 | `77018feaaa65bb22dea47c784e8c059b0ccc09cd6dc7442b79cce83f3170985f` |

**PDB format rather than mmCIF**: fixed columns are the simpler of the two to read strictly. Fetched
twice and identical both times. Its primary citation is A. Morton and B. W. Matthews,
*Biochemistry* **34**, 8576 (1995), doi:10.1021/bi00027a007; its `REMARK 1` cites A. Morton, W. A.
Baase and B. W. Matthews, *Biochemistry* **34**, 8564 (1995), the binding thermodynamics.
**What the file holds besides the protein**: 136 waters, two chloride ions (`CL` 173 and 178), one
2-hydroxyethyl disulfide (`HED` 170, a crystallisation additive) and the benzene (`BNZ` 400); no
alternate locations, no hydrogens, and residues 163–164 (Asn, Leu) not located (`REMARK 465`).

**The congeners, recorded and not yet used.** RCSB's search for entries whose citation is
doi:10.1021/bi00027a007 (Morton and Matthews 1995) returns nine, each the same L99A pseudo-wild
type at 1.80 Å with `CL` and `HED`:

| entry | ligand | code |
| --- | --- | --- |
| 181L | benzene | `BNZ` |
| 182L | benzofuran | `BZF` |
| 183L | indene | `DEN` |
| 184L | isobutylbenzene | `I4B` |
| 185L | indole | `IND` |
| 186L | n-butylbenzene | `N4B` |
| 187L | p-xylene | `PXY` |
| 188L | o-xylene | `OXE` |
| 1NHB | ethylbenzene | `PYJ` |

The same search for doi:10.1021/bi00027a006 (Morton, Baase and Matthews 1995) returns 181L only,
as a secondary citation. Queried 2026-10-04 on `rcsb_primary_citation.pdbx_database_id_DOI` and
`citation.pdbx_database_id_DOI`, which agreed; the ligands are each entry's `HETNAM`. Ethylbenzene's
`PYJ.cif` is already above.
