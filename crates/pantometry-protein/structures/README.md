# Four structures, unmodified

These are Protein Data Bank entries as `files.rcsb.org` serves them, byte for byte. Nothing here
is trimmed, reordered or regenerated, which is the point: a fixture that has been through a script
can no longer be checked against the thing it came from, and `a_protein_a_crystallographer_measured.rs`
asserts each file's own `HEADER` line names the entry it is supposed to be.

| file | | residues | why it is here |
| --- | --- | --- | --- |
| `1CRN.pdb` | crambin | 46 | The small one. A whole normal-mode analysis in under a second, so the check runs in a debug build too |
| `1UBQ.pdb` | ubiquitin, 1.8 Å | 76 | Small, well refined, and not a special case the way crambin is |
| `193L.pdb` | hen lysozyme, 1.33 Å | 129 | A modern refinement of a protein this model should get right |
| `4LYZ.pdb` | hen lysozyme, 1975 | 129 | **The control.** The same protein, refined before `REFINEMENT. PROGRAM` had a value |

## The pair is the reason there are four

Predicted fluctuations against the deposited B-factors, at the anisotropic network model's
published 15 Å cutoff:

| | correlation |
| --- | --- |
| 193L | **0.659** |
| 1UBQ | **0.571** |
| 1CRN | **0.395** |
| 4LYZ | **−0.404** |

The last two lines are the interesting ones. 4LYZ and 193L are the same protein with the same
fold and the same 129 residues, and the model says opposite things about them — so what the
comparison is sensitive to is the **measurement**, not the model, and a check that only ever
reported the agreeing cases would not have shown that. 4LYZ's B-factors run from 0.46 to 12.35 Å²
against 193L's 9.1 to 32.9, and its `REMARK 3` refinement program is `NULL`.

Crambin's 0.395 is low and is not a defect either: it is 46 residues held by three disulfides,
about as rigid as a protein gets, and what little its B-factors vary by is mostly crystal contact
rather than internal motion. A cutoff chosen to make that number larger — 7 Å gives 0.685 — would
be a fit, and the cutoff here is the published one.

## Terms

Protein Data Bank coordinate files carry no usage restriction; the wwPDB releases them into the
public domain (CC0). Cite the entries themselves rather than this directory:

- **1CRN** — Hendrickson & Teeter (1981), *Nature* 290, 107. <https://doi.org/10.2210/pdb1CRN/pdb>
- **1UBQ** — Vijay-Kumar, Bugg & Cook (1987), *J. Mol. Biol.* 194, 531. <https://doi.org/10.2210/pdb1UBQ/pdb>
- **193L** — Vaney et al. (1996), *Acta Cryst.* D52, 505. <https://doi.org/10.2210/pdb193L/pdb>
- **4LYZ** — Diamond (1974), *J. Mol. Biol.* 82, 371. <https://doi.org/10.2210/pdb4LYZ/pdb>

## What the B-factor comparison is worth, measured

Two predictors that are one line of geometry each, against the same deposited B-factors at the
same 15 Å cutoff:

| | model | distance from centroid | neighbours within the cutoff | residue index |
| --- | --- | --- | --- | --- |
| 1CRN | +0.395 | +0.402 | +0.423 | +0.375 |
| 1UBQ | +0.571 | **+0.804** | **+0.768** | +0.346 |
| 193L | +0.659 | **+0.699** | **+0.697** | +0.262 |
| 4LYZ | −0.404 | −0.301 | −0.267 | −0.234 |

**The model does not beat either geometric predictor on any of them.** That is not a defect in the
implementation and it is not news in the field — a residue's temperature factor is mostly a
statement about how buried it is — but it is the kind of thing a crate reports about itself only if
somebody writes the check, and `the_trivial_predictors_do_as_well` is that check.

What follows is that these four structures are a **necessary** test and not a sufficient one. The
sufficient one needs the two files below.

## And the pair that a scalar cannot fake

| file | | residues | |
| --- | --- | --- | --- |
| `4AKE.pdb` | adenylate kinase, **open** | 214 (×2, a dimer) | what the model is given |
| `1AKE.pdb` | the same enzyme **closed** on an inhibitor | 214 (×2) | what it is asked to have predicted |

Superposed, the two conformations are **7.1 Å** apart: the enzyme folds two lids over its
substrates. The model sees only the open one. Its **single softest mode** overlaps the observed
motion by **0.799**, and ten modes reach **0.966**, in a space of `3 × 214 = 642` dimensions where
a direction chosen without looking scores `1/√642 = 0.039` — measured over two hundred draws at
`0.0312`, against the closed form `√(2/πn) = 0.0315`.

That is a statement about a **direction**, and it is the claim an elastic network model is actually
making. A per-residue scalar cannot produce it, which is why both kinds of check are here.

- **4AKE** — Müller et al. (1996), *Structure* 4, 147. <https://doi.org/10.2210/pdb4AKE/pdb>
- **1AKE** — Müller & Schulz (1992), *J. Mol. Biol.* 224, 159. <https://doi.org/10.2210/pdb1AKE/pdb>
