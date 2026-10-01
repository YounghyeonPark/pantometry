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
