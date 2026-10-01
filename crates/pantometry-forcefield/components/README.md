# One chemical component, unmodified

`AIN.cif` is aspirin (acetylsalicylic acid, C₉H₈O₄) as the wwPDB Chemical Component Dictionary
serves it, byte for byte:

| file | component | source | fetched |
| --- | --- | --- | --- |
| `AIN.cif` | `AIN`, 2-(acetyloxy)benzoic acid | <https://files.rcsb.org/ligands/download/AIN.cif> | 2026-10-01 |

SHA-256 `e24ecc21a59aa8e26e1143234f92f3e1e309a589508f38d4d91e078a456c57ed` as fetched. Nothing is
trimmed, reordered or regenerated: a fixture that has been through a script can no longer be checked
against the thing it came from. The tests that need a malformed file make one from this text in
memory, one edit at a time, and assert the edit's anchor occurs exactly once before making it.

The entry was created by RCSB in 2004 and last modified in 2020 (`_chem_comp.pdbx_modified_date`).
Its ideal coordinates are the ones the crate uses; its model coordinates come from PDB entry `1OXR`
and are used only by the test that checks the fallback.
