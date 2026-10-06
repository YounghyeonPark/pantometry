# Changelog

Notable changes, in the format of [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This workspace follows [semantic versioning](https://semver.org/). It is `0.x`, so the API is
explicitly not stable and a minor bump may break you. The first consumer exists now, and it
has already found forty-nine places it is awkward, forty-three of which have been changed — see
`app/pantometry-world/FRICTION.md`.

**Entries below 0.16.0 name crates as `pantometry-*` and they were published as `dualis-*`.** The
tree was renamed in one pass rather than annotated version by version, because the alternative is a
changelog in two vocabularies; what a reader needs is the one sentence above and the fact that
crates.io holds both.

Entries record what was *found* as well as what was added, because several of the more useful
changes here were corrections to a mistaken assumption rather than new features. The commit
messages carry the full account.

`CONTRIBUTING.md` holds the gate and indexes the sixteen times it has reported a pass it had not
earned; `CLAUDE.md` carries the kinds they were, and why a `set -euo pipefail` pasted into a shell
protects nothing.

## [Unreleased]

### Added

- **`pantometry-forcefield`, the fourteenth domain: a small molecule, read and typed — and not yet
  an energy.** The first step towards molecular mechanics on drug-sized molecules, and deliberately
  only the part that can be checked before any energy exists. `Component::from_ccd` reads one wwPDB
  Chemical Component Dictionary entry (mmCIF) with a strict loop parser and no dependency: a missing
  column, an element outside H C N O F P S Cl Br I, a bond order other than `SING`/`DOUB`/`TRIP`, a
  bond naming an atom that is not there, a repeated atom or bond, a row one value short, or several
  components in one file are each a `CcdError` that names what it refused. Positions are the
  entry's ideal coordinates, or its model coordinates **for every atom** when any ideal one is
  missing, because the two are in different frames — twenty ångström apart for aspirin.
  `uff::assign` gives every atom its Universal Force Field type from element, bond orders and the
  dictionary's own aromatic flags, with no ring perception; `UffType::parameters` is Table I of
  Rappé et al., *J. Am. Chem. Soc.* 114, 10024 (1992), for the twenty-two types those ten elements
  need, in SI, with the torsion barriers beside it. `Molecule` is the `Domain`: atoms as `Bodies`
  with the dictionary's names and bonds, so `capture` draws it ball-and-stick without naming it.
  **It does not move**: `step` leaves every atom in place, the ledger is empty, and both say so.
  Checked against aspirin's (`AIN`) own stated formula and net charge, one bond per hydrogen,
  every heavy atom's valence with aromatic bonds at 1.5, the types a chemist reads off its
  structure, a textbook molecule for each rule aspirin does not reach (acetamide against
  methylamine for the amide nitrogen, HCN, CO, CO₂, pyridine, thiophene, furan, phosphoric acid,
  phosphine), the tetrahedral, trigonal and linear angles `θ₀` has to be, and UFF's own rule that
  the non-bonded parameters belong to the element. Left out and said so where it would be missed:
  the GMP electronegativity χ, whose source paper is not yet verified, so nothing here computes a
  UFF natural bond length; and the hypervalent sulfur types, so a sulfone types as `S_2` until the
  energy step refuses or extends it.
- **`pantometry-forcefield` computes an energy: UFF's bond stretch, van der Waals and
  electrostatics, with analytic forces — and still does not move.** `ForceField` gives the energy
  by term and the force on every atom: harmonic stretch about the natural length of eq 2 with the
  bond-order (eq 3, aromatic 1.5, amide C–N 1.41) and electronegativity (eq 4) corrections and the
  force constant of eq 6; Lennard-Jones 12-6 with geometric combination (eqs 20–22); Coulomb at
  332.0637 kcal mol⁻¹ Å e⁻² (eq 43); 1-2 and 1-3 pairs excluded. **`r_EN` is subtracted, not added
  as eq 2 prints it**, as Open Babel and RDKit do and as the paper's minimised dimethyl ether and
  trimethylamine favour — and the paper's own amide sentence favours the printed plus: with the
  amide carbon now typed `C_R`, as that sentence's "C_R and N_R single bond radii" require
  (`uff::assign` typed it `C_2` before), n = 1.41 gives 1.35684 Å with the minus and 1.36845 Å
  with the plus, against the paper's 1.366. Recorded beside the decision for the minimisation step
  to settle. χ is Open Babel's copy, not read from its source paper; the paper's worked Si–O
  correction (0.0533 Å) is the one partial check, and one natural length per element against
  `C_3` holds the transcription. Charges are an input and default to zero until the charge model
  arrives. A sulfur whose bond orders sum past two — sulfoxide, sulfone, sulfonium — is now refused
  by name (`Unsupported`) instead of getting a divalent radius, and the message says which. An
  N-acyl aromatic nitrogen is not treated as an amide (order 1, carbon `C_2`), a recorded choice.
  Checked against closed forms — 1.514 Å for `C_3–C_3`, `k` and `r_IJ` by hand for aromatic and
  amide `C_R–N_R`, the well's depth and zero crossing, the Coulomb constant against CODATA's ε₀,
  the exclusions counted on chains and rings — and the forces against central differences of the
  energy on aspirin, to a tolerance built from per-term rounding bounds, with every atom required
  to have most of its terms well above it; each of thirty-one sabotages was caught. **Found on the way:** the dictionary's
  ideal coordinates for aspirin put the acetyl oxygen O4 1.645 Å from ring hydrogen H1, and that
  one contact is 134 of aspirin's 173 kcal/mol of van der Waals energy there.
- **`pantometry-forcefield` has all six of UFF's terms: angle bend, torsion and inversion join
  the three, with analytic forces — and it still does not move.** New module `angular`, every
  formula cited to its page and equation of Rappé et al. (1992) and read off the scanned pages.
  Angle bend: the general three-term Fourier form of eqs 11–12, the special forms of eq 10 for
  linear and trigonal centres, and the force constant of eq 13 from the two bonds' natural
  lengths; which form a centre takes is decided by its Table I θ₀ (180° linear, exactly 120°
  trigonal), so `N_2` and `O_R`, whose angles p. 10028 says were fitted, keep the general form.
  **Eq 10 with n = 1 is a misprint**: `K (1 − cos θ)` has its minimum at 0°, and a linear centre's
  is at 180°, so the linear term is `K (1 + cos θ)`. Torsion: eq 15 shared among the torsions
  about each bond, with eq 16 and Table III for sp³–sp³, eq 17 for sp²–sp², and the paper's
  sp²–sp³, group-6 and propene cases in a stated order, nothing about an sp centre, and resonant
  types counted as sp². The paper's "set to zero" near a straight central angle is done as a
  smoothstep from 170° to 180° so that neither the energy nor the force jumps — a number the
  paper does not give. Inversion: eq 18 about each of the three bonds, divided by three; `C_2` and
  `C_R` at 6 kcal/mol or 50 bonded to `O_2`; **no nitrogen at all**, as the paper says; and
  phosphine's 22 kcal/mol barrier with ω₀ = 84.4339° derived from PH₃'s 93.8°. Checked against
  closed forms: every angle form zero, flat and of curvature K at θ₀ for all twenty-two types; K
  for methane's H–C–H (75.4988 kcal mol⁻¹ rad⁻²) and for H–C–Cl by hand; every torsion row by
  hand; ethane's nine torsions summing to exactly 2.119 kcal/mol eclipsed; eq 17 at 38.974 for
  ethylene's twist and 26.948 for benzene; the six-fold row on methyl isocyanate; the trigonal
  form away from 120°; K read off the force field at bond orders 1.5 and 2 (222.595 and 170.396);
  the switch against its documented formula on both central angles, across its edge, and its
  slope at both ends; three-ring paths; inversion zero when flat and by hand when not, including
  an asymmetric pyramid, and phosphine at 0 and 22. Forces against central differences on
  aspirin, on a probe molecule built to put torsions inside the switch, a six-fold torsion and
  two linear centres off their minimum and a pyramidal phosphorus, and on phosphine near its
  minimum where `cos ω` ≈ 0.1, with the tolerance model extended to measured third derivatives
  for the angular terms; 79 tests in the crate. Each of thirty-four sabotages was caught —
  eight of them, from a review, passed the first version of these tests and were each closed
  with the check that now fails them. **Found on the
  way:** the paper's one explicit rule for an sp³ oxygen on an sp² atom — eq 17, n = 2, φ₀ = 90° —
  puts the minimum of an ester's, an acid's and an anisole's C–O torsion at perpendicular, 10
  kcal/mol below planar; implemented as written, it is 29.7 of aspirin's 30.7 kcal/mol of torsion
  energy at the dictionary's planar geometry. The paper also says its minima are DREIDING's, which
  is what step 1d's anisole barrier against Table II will test. Ethane's rigid-rotation barrier is
  3.082 kcal/mol (torsion 2.119, van der Waals 0.963); Table II's 2.90 is relaxed, and that
  comparison waits for the minimiser.
- **`pantometry-forcefield` minimises, and a run of `Molecule` is the molecule relaxing frame by
  frame.** New module `minimise`: L-BFGS (memory 8) with an Armijo backtracking line search,
  chosen over FIRE because a step is taken only if it lowers the energy, so **the energy never
  rises** — tested step by step — and over steepest descent because a molecule's Hessian spans
  bonds at ~700 kcal mol⁻¹ Å⁻² to rotors at ~1. Converged means the largest per-atom force is at
  most the tolerance; when no step lowers the energy in floating point it says `Stalled` and does
  not move. `ForceField::minimise`, `Molecule::minimise`, and `Domain::step` as one iteration, with
  `max force`, `rms force`, `converged` and `minimiser steps` readings and a checkpoint that
  restores the minimiser's history too. `DihedralRestraint` holds a dihedral for relaxed scans.
  Only `+ − × ÷ sqrt` in the unrestrained minimiser (a restraint adds `atan2`): aspirin's minimum
  is bit-identical in debug and release on one machine (FNV-1a `dac2edd5b0a72aa3` after 443
  steps, both). Not across platforms: the energy calls the platform's `cos` per torsion and `ln`,
  `sin`, `cos`, `asin` at construction, so the digest is printed, not pinned. Checked against minima known exactly, each
  within the displacement its final force allows (`|F|/k`, or `|F|/λ_min` from a measured Hessian):
  a C–C bond at 1.514000000000000 Å, water at its natural lengths and θ₀ = 104.51°, methane the
  regular tetrahedron; methanol's minimum the same from a moved and turned start; the restraint's
  force against central differences. **Found while writing it:** the step compared its input
  with the last point through a metres → Å round trip, which is not bit-exact, so every step
  discarded its history and the minimiser was steepest descent — still monotone, still
  converging: aspirin had not converged after 6000 steps (18.806 kcal/mol, force 0.12); compared
  in metres, 443 steps to 18.573. Nothing caught it; `the_curvature_history_survives_between_steps`
  does now. **And a review found the L-BFGS itself was not checked**: a reversed pairing in the
  two-loop recursion converged aspirin in 1562 steps to a different minimum (18.8037 kcal/mol), and
  a fixed `H₀` and an ignored history were caught only by accident. `two_loop` is now compared with
  the dense BFGS update (agreement 2.2e-16), and L-BFGS with the same minimiser at memory 0 —
  steepest descent with the same search, 13904 steps to 18.8037 against L-BFGS's 443 to 18.5734;
  an ignored history would make the two runs the same arithmetic. The closed-form minima now check
  the force independently of the minimiser's report — water and methane passed a minimiser
  stopping at 1e4 times the tolerance — and accept a stall only below the force at which rounding
  makes a correct minimiser stall, `2 √(λ_max δE)`, instead of failing it. **The paper's Table II,
  relaxed** (scans at 30°, golden-section refinement, restraint energy ≤ 1e-14 kcal/mol at every
  asserted row's extrema): ethane 2.8976 against 2.90, CH₃NH₂ 1.979 against
  2.0, CH₃PH₂ 2.027, CH₃OH 1.036 against 1.0, CH₃SH 1.300 against 1.3, trans HO–OH 1.678 against
  1.7, trans and cis HS–SH 6.822 and 7.216 against 6.8 and 7.2 — eight rows within half the last
  printed figure, asserted. Not asserted, with reasons recorded: cis HO–OH 6.52 against 6.6;
  acetaldehyde 0.17 against 0.83 (the propene rule's six terms cancel to a constant); isoprene and
  ethylbenzene, where the table does not say which rotation (neither the central bond nor the
  methyl agrees); and anisole and thioanisole — **19.9 and 14.5 under the group-6 rule as printed,
  0.26 and 1.50 with it made planar, and 3.628 and 1.666 with the heteroatom typed `O_R` and
  `S_R`, against the paper's 3.6 and 1.7.** **What the eight asserted rows can see is the paper's
  fit, not its method**: p. 10029 says the `C_3`–sp³ barriers were fitted and H₂O₂ and H₂S₂ are
  compromise values, and the rows it names as tests of eq 17 are the ones that do not agree. They
  are blind to the combination rule (every 1-4 pair in them is H···H; propane's C–C in the structure
  comparison is the asserted check that sees a C···H pair, and arithmetic combination moves it to
  1.5271 against 1.526), to `V_O` anywhere from about 0.006 to 0.021 and to `V_S` = 0.448 (held
  instead by the transcription check of Table III), and to eq 17's constants, held by hand values.
  The paper gives no calculated number for eq 17 at bond order above one — it fitted the 5 and 4.18
  to vibrational modes and to N,N-dimethylformamide's ΔH‡ = 19.7 kcal/mol, with no residual — so
  none is asserted; DMF's relaxed barrier, 19.807 kcal/mol, is printed for information. The paper's methyl vinyl ether points the same way: its
  `O_R` angle was fitted to that molecule's 118.3°, typed `O_R` it relaxes to 118.05°, typed `O_3`
  to 107.8°. **The paper's Figures 4–7 under both signs of `r_EN`, and without it**: of 24
  heteronuclear bonds the subtracted sign is nearer the paper than the added one for all 24 —
  which alone would also hold if the paper had used no correction, so the correction dropped is
  the third reading: subtracted is nearer than none for all 24, and none nearer than added for all
  24. 36 measures agree to the printed figure —
  dimethyl ether C–O 1.4096 against 1.410 (1.4497 added), trimethylamine C–N 1.4708 against 1.471,
  N-methylformamide C–N 1.3660 against 1.365; added, a C=O is 0.04 Å long. **Both questions are
  measured and neither is changed in the code**: `Variant` switches them for measuring, its default
  is unchanged, and the maintainer decides. Aspirin relaxed from the dictionary's ideal geometry:
  223.10 → 18.57 kcal/mol (van der Waals 173.28 → 13.68, torsion 30.67 → 1.63), O4–H1 1.645 →
  4.242 Å, heavy-atom RMSD to the crystal (`1OXR`) 0.905 Å at the start and 1.113 Å relaxed, with
  the ester and the acid's hydrogen turned perpendicular by the group-6 rule; typed `O_R`, 0.894
  Å. Twenty more Chemical Component Dictionary entries in `components/`, each with its SHA-256
  and held to its own formula; nine molecules the dictionary lacks hand-built under `tests/` and
  labelled so; N,N-dimethylformamide among the twenty, for the eq 17 measurement. Fig 6's dimethyl ether
  methyl, printed 1.109 and 1.113, was compared as one mean against 1.111; it is now the shortest
  and longest C–H, 1.1094 and 1.1130. 97 tests in the crate; each of seventeen sabotages was caught,
  and the review's eleven were rerun before and after the changes it asked for.
- **A scene for the fourteenth domain: `33-aspirin-relaxing-out-of-a-clash`, and an example,
  `aspirin_relaxes`.** `molecule` is the scene format's twenty-second kind: one Chemical Component
  Dictionary entry named by `ccd`, a path beside the scene served the way `pdb` and `stl` are, so a
  preset carries it (`WithFiles`) — the generator found the file by the rule it already had. The
  build refuses an entry the force field cannot describe, naming the domain, the file and the atom,
  because `Molecule::new` never fails and a scene asking for a relaxation would otherwise run with
  nothing moving and NaN in the readings. The scene photographs every one of 500 minimiser steps,
  so frame to frame is step to step: 185.672 → 29.583 kcal/mol, converged at step 454 with the
  largest force 7.65e-5 kcal mol⁻¹ Å⁻¹, drawn ball-and-stick with the dictionary's 21 bonds.
  **Checked against what minimisation guarantees and a published threshold, not a pinned energy**:
  the energy never rises and a frame that moved is strictly lower; it ends converged; and the
  O4–H1 contact, 1.645 Å at the start, opens past Bondi's radii less MolProbity's 0.4 Å
  serious-clash overlap, 2.32 Å — to 2.666 Å, which is **still inside the Bondi contact of
  2.72 Å** and inside UFF's own zero crossing for the pair, 2.831 Å, and is said rather than
  asserted. The entry beside the scenes is a byte-for-byte copy of the crate's, held by a test.
  The example prints the six terms before and after, the steps and the contact, and checks the
  atoms against the entry's own `_chem_comp.formula`, the energy at every step, electrostatics
  exactly zero with no charges given, the largest force recomputed from the forces, and the
  contact against the same threshold from inside it. Each of twenty sabotages was caught. All
  fourteen domains have a scene; 33 scenes, 22 kinds, 29 tiles, 16 examples run by CI and 17
  example files. The scene is filed under the chooser's protein area, renamed *Molecules*: a
  fourteenth area put the last one's tile below the fold of a 950-point window, which
  `every_area_of_the_start_screen_lays_out_at_every_width` caught.
- **`pantometry-forcefield` computes charge-equilibration (QEq) charges, the ones UFF prescribes —
  as an option; the default stays zero.** New module `qeq`, from Rappé and Goddard, *J. Phys.
  Chem.* 95, 3358 (1991), read off the scanned pages: Table I's χ, J, R and ζ for the ten
  elements; the shielded Coulomb `J_AB` as the exact Coulomb integral between normalised `ns`
  Slater densities (n = 1 for H, 2 for C–F, 3 for P–Cl, 4 for Br, 5 for I), a finite sum of
  closed forms whose `∫ v^q e^(−δv)` part is a one-signed power series in `(ζ_A − ζ_B)R`, so equal
  and nearly equal exponents lose nothing; Table I's idempotential on the diagonal, not the Slater
  self-integral (carbon's 2s self-repulsion is 8.615 eV against J⁰ 10.126); the dense linear solve
  by Gaussian elimination with partial pivoting, no dependency; the eq 5/5′ ranges with the
  paper's fix-at-the-boundary procedure and eq 13; and hydrogen's charge-dependent `ζ_H = ζ⁰ + Q_H`
  (eq 20), in **every** J it enters, with `J_HH(Q)` (eq 21), iterated from zero to 10⁻¹⁰ e.
  `qeq::charges(&component, &at)` takes the dictionary's formal charges as the total;
  `ForceField::with_qeq_charges(&at)` puts them in the electrostatic term, **fixed at that
  geometry**: the force treats them as constants, which is the paper's usage. **Three things the
  paper does not say plainly, each settled by reproducing its numbers.** (1) **Eq 12 prints
  `C_1j = Q_j`; it is `C_1j = 1`**, the row that is eq 9. (2) **Table I's printed ζ is not the ζ
  the paper's charges were computed with.** It is eq 17 with the fitted λ = 0.4913 — the λ each
  printed value implies is 0.49123–0.49129 for C, N, F, P, S, Cl, Br and I — while hydrogen's
  1.0698 is λ = ½ and oxygen's 0.9745 is neither (λ = 0.49280; it is 0.4913 at R = 0.667 Å, not the
  printed 0.669). The paper rounds to λ = ½ (eq 17′), and only that reproduces Table III: HF's
  hydrogen 0.4623 against 0.462, 0.4568 with the printed ζ; H₂O 0.3531 against 0.353 (0.3453); NH₃
  0.2412 against 0.243 (0.2295); CH₄ 0.1497 against 0.149 (0.1342). (3) **The paper's charges are
  its iteration, not the stationary point of its eq 23**: eq 23's derivative puts `1.5 Q/ζ⁰` where
  the iteration has `Q/ζ⁰`, and gives HF's hydrogen 0.405. **Provenance closed:** the crate's GMP
  χ, transcribed from Open Babel, equal Table I exactly for all ten elements, which prints them
  citing the GMP paper (its ref 9) — `uff`'s note that they were not from their primary source
  is replaced. Checked against closed forms: every `ns` self-repulsion for n = 1 to 5 against the
  exact fractions 5/8, 93/256, 793/3072, 26333/131072 and 43191/262144 of ζ (from
  `J = 2 ∫ p q / r`, a different sum from the code's); two unequal 1s exponents at one centre;
  Roothaan's equal-exponent 1s–1s `J(R)` at eight distances; `J·R = 1` to 8 ε at 120 bohr for
  every n pair (`14.39964` eV Å, the paper's 14.4); six points of very unequal exponents at a
  distance against two references that share no code with the crate — a SymPy closed form for
  1s–1s and a 40-digit Fourier quadrature — to 10⁻¹² relative; symmetry under exchanging the atoms
  (the arithmetic is asymmetric); positivity and `J < 1/R`; continuity at R = 0; nearly equal
  exponents to first order; eq 18 for C–O, C–F and C–Cl; both ends of a heavy atom's range (F
  pushed below −1, C above +4) and of hydrogen's (±1, through the iteration), with eq 13 for the
  rest; total charge to `n ε Σ|Q|` at 0, ±1; aspirin's charges unchanged by a rigid motion to
  6e-16 e; eq 8's equal chemical potentials to 2e-12 eV; the matrix positive definite; and every
  error reached by an input that should reach it. **Against the paper**: Table II's eighteen
  alkali-halide charges by eq 18 (Na, K, Rb with Cl, Br, I, both λ columns, Huber–Herzberg r_e)
  to 0.00048 e against 0.0005 — the only physical check of the 4s and 5s integrals; all four
  Table III rows and 36 Table IV rows over seventeen molecules at stated experimental geometries
  (NIST CCCBDB and WebBook, each citing its source), every one within half its last printed
  figure plus `Σ |∂Q/∂g| δg`. For HF and HCl δ is 10⁻⁴ Å, Huber–Herzberg's precision, so they are
  held to the printed figure alone: HF 0.4623 against 0.462. For polyatomics δ is 0.01 Å per
  length and 1° per angle — **an assumed allowance, stated as one**: it was set after a Python
  reproduction had shown the misses, not derived from these molecules' r_e/r_0/r_s spreads, and
  the smallest common scale of it at which every row passes is 0.962 (formamide's carbon, 0.0125
  against 0.0128); `∂Q/∂g` is measured with the code under test. None missed. Every molecule's
  charges are checked against eq 8 at default settings. **The paper disagrees with itself**
  twice on p. 3361: its text gives H₂CO's carbon 0.21 and H₂O's hydrogen 0.36, where Table IV
  prints 0.19 and 0.35 and Table III 0.353; the tables are compared (0.1968 here). An unlabelled
  Table IV row under PH₃ is left out. The paper's "six to ten iterations" is aspirin's 6 at 10⁻²
  e and 8 at 10⁻³; 26 at 10⁻¹⁰. **Aspirin at its UFF minimum**: O1 −0.703, O2 −0.445, O3 −0.542,
  O4 −0.496, C7 +0.603, C8 +0.589, HO1 +0.354; UFF's electrostatic energy with them at that fixed
  geometry −83.88 kcal/mol, against a total of 29.58 without charges. No asserted number moved.
  **A review found the first version's tests could not see five defects**, each now a test that
  fails without it and was shown passing before: n capped at 3 (no check reached 4s or 5s), a
  default tolerance of 10⁻² (HF's hydrogen then 0.4601, inside the generous diatomic allowance),
  the upper clamp fixing at the lower bound, a branch of the integral returning zero, and a NaN
  total charge coming back `Ok` with NaN charges (`max` discards a NaN and `<=` is false of one;
  now `QeqError::NonFinite`). **And the new checks found a defect in the integral**: expanding
  the shielding polynomial about the origin of the integration rectangle loses `R^(2m−1)` to
  cancellation when one density is compact and the other diffuse, and 4s(ζ 8) against 2s(ζ 0.7)
  at 8 bohr came out 5.4 × 10⁻¹⁰ relative wrong (3.7 × 10⁻⁷ of the shielding term). Above
  `|ζ_A − ζ_B| R` = 2 it is now expanded about the corner where the weight is, and the
  antiderivative branch is gone; the switch is where the two expansions' f64 errors, measured
  against the same formula at 60 digits over nineteen cases, are comparable, and the worst either
  then makes is 4.0 × 10⁻¹⁴ relative. **Found on the way:** two neutral hydrogens are singular at
  0.7915 Å, where `J_HH` equals J⁰_HH — longer than H₂'s 0.741 Å bond, so at its bond length eq
  6′'s form for moving charge between them is negative at fixed orbitals and QEq's zero for H₂ is
  a stationary point, not a minimum. Each of sixteen sabotages of the first version was caught —
  eleven of them by the paper's own tables — and each of the eight for the review's items since; chlorine's n,
  which Table IV cannot see (n = 2 gives HCl 0.3203, nearer the printed 0.32), Table II does
  (NaCl, KCl and RbCl off by 0.0011–0.0013). Sabotage runs use `--no-fail-fast`: the first did
  not, and `cargo test` stopped at the first failing binary, so "caught by the paper" was not a
  measurement until it was rerun. 122 tests in the crate.
- **`pantometry-forcefield` has implicit water: generalized Born solvation, OBC II, with its
  analytic force — opt-in; the default stays vacuum.** New module `solvation`, from Onufriev,
  Bashford and Case, *Proteins* 55, 383 (2004), read off its pages: eq 2 over every ordered pair,
  `i = j` included, with Still's `f` of eq 3; the solute's dielectric 1 (p. 384), the same vacuum
  UFF's Coulomb term is in; ε_w = 80, the paper's "80 for water at 300 K" (p. 384); optional
  Debye–Hückel κ = 0.316 √[salt] Å⁻¹ (p. 384), default none; `ρ̃ = ρ − 0.09 Å` (p. 385); eq 6 with
  OBC II's α, β, γ = 1.0, 0.8, 4.85 (eq 8), and OBC I and HCT's eq 4 beside it.
  `ForceField::with_generalized_born()` adds `Energy::solvation` to the total and its force to the
  forces, computed from the force field's charges, so it is zero until QEq's or a caller's are
  given; `Molecule` and its readings are unchanged. **The force is analytic at fixed charges, by
  the chain rule through every Born radius**, and against central differences of the energy on
  aspirin (OBC II, OBC II in 0.15 M salt, HCT) and on a cluster built to reach every branch of the
  pair integral, to a tolerance from a measured third derivative and a rounding bound: the worst
  miss is 0.018 of it, and every coordinate's force is over a hundred tolerances. **Secondary, and
  said so**: eq 5's closed pair form and the scale factors (H 0.85, C 0.72, N 0.79, O 0.85, F 0.88,
  P 0.86, S 0.96, else 0.8) are OpenMM's `customgbforces.py`, and the radii are AMBER's "Bondi"
  set the paper fitted with (Table I, p. 387; its OBC II trajectory is "parm99, GB^OBC (II), Bondi",
  Table II) — **not the mbondi2 set AMBER and OpenMM use for igb=5**, which puts a hydrogen on
  nitrogen at 1.3 Å, and **not Bondi's own table** (O 1.5 against
  1.52, F 1.5/1.47, P 1.85/1.80, Cl 1.7/1.75) and has **nothing for Br or I**, where AMBER's 1.5 Å
  default would make them smaller than carbon; Bondi's 1.85 and 1.98 Å are used, not read in
  primary and outside anything OBC was fitted to. **The pair integral was re-derived, not taken on
  trust, and OpenMM's form is missing a term**: when atom i's sphere lies wholly inside j's scaled
  sphere the shells from ρ̃_i to `s − r` add `1/ρ̃_i − 1/(s − r)`, which its `L = max(ρ̃, |r − s|)`
  leaves out; the derivative needs only the explicit `∂/∂r`, because the moving limits' terms
  vanish or cancel. Checked against closed forms: eq 5 against a quadrature in the test in every
  branch to 10⁻¹²; a lone ion's radius `ρ − 0.09 Å` exactly and Born's energy, with salt too;
  two ions 10, 20 and 40 Å apart against two Born terms and the screened Coulomb, off by 1.1e-2,
  7.0e-4 and 4.4e-5 kcal/mol within a bound from the descreening of 1.8e-2, 8.9e-4 and 4.9e-5;
  two coincident charges as one of their sum; eq 6 typed from the paper for all three
  rescalings; the paper's 30 Å bound for carbon (30.41; a carbon deep in a 50 Å sphere 30.30);
  HCT's radius NaN, not negative, when overlap pushes `I` past `1/ρ̃`; rigid motions; zero net
  force. **A charge inside a dielectric sphere**, as one dense atom with a = 3 Å and the charge
  1.2 Å from its centre: HCT gives the Coulomb-field radius `½[a/(a² − d²) + ln((a+d)/(a−d))/(2d)]⁻¹`
  = 2.667143 Å to 10⁻¹² for every probe radius — only with the inside term — and **Kirkwood's
  exact series at ε = 80 gives 2.522651 Å: the Coulomb-field approximation is 5.73% long in R and
  5.42% short in energy**, the model's own error, pinned as a measurement, not tolerated away (the
  brief said 5.8%; it is 5.73). OBC II puts the same charge at 3.16–5.32 Å depending on the probe's
  radius. The nonpolar `0.005 kcal mol⁻¹ Å⁻² × SASA` (p. 392) is computed by Shrake–Rupley on a
  golden-spiral point set, Bondi radii plus a 1.4 Å probe (a choice; the paper names none), and is
  **not in the force field**: a counted surface has no gradient. Its worst error over eight
  overlapping pairs against the exact two-sphere area is 1.86%, 0.49%, 0.29%, 0.096% and 0.025% at
  100 to 25 600 points — no clean rate, and none is claimed. **Real molecules**, relaxed by UFF in vacuum, QEq charges at that minimum,
  against FreeSolv v0.52 (Mobley and Guthrie 2014; each value's own reference beside it): **aspirin
  ΔG_GB −18.22 kcal/mol**, plus 1.86 nonpolar, −16.36 against −9.94. Over eighteen molecules the
  error has mean −2.14 and RMS 3.36 kcal/mol, Pearson r 0.71; methylamine −4.75 against −4.55 and
  DMF −8.17 against −7.81, but methanol −9.15 against −5.10, dimethyl ether −7.81 against −1.91,
  anisole −6.51 against −2.45, propane −1.15 against +2.00, and acetamide −6.99 against −9.71.
  **What that means: with QEq charges, OBC II over-solvates** — mean −2.1, RMS 3.4 kcal/mol, the
  ethers and anisole by about −4 to −6 — **so absolute binding energies built on it are not
  quantitative.** Reported, not asserted: α, β, γ and the radii were fitted to Poisson–Boltzmann
  energies of proteins with AMBER charges, not QEq's. Each of fifteen sabotages was caught — ρ̃ without the offset, the inside term,
  OBC I for OBC II, the `i = j` term, Still's exponential, the chain rule through R, the `sech²`
  in dR/dI, the sign of a term in dI/dr, κ in g′, κ doubled, κ halved in the energy, the cross
  term counted once, HCT's NaN guard, the surface on the bare radius, and solvation left out of
  the total — each restored byte for byte and its SHA-256 checked. **A review then found five
  kinds of error those tests passed, each shown passing before and caught after**: Still's
  exponent 4 → 8 in energy and gradient together (the coincident ions sit at r = 0, the far ones
  inside the allowance, and the force check sees only consistency); κ dropped from the pair terms
  only; `with_solvent_dielectric` ignoring its argument; the surface tension 0.005 → 0.05, which
  passed the whole crate; and seven radii and scale factors moved, which passed the whole crate
  too, because the lone-ion test read ρ from the code it checked. Now: two bare ions at 1, 2.3
  and 4 Å, ε 80 and 4, no salt and 0.15 M, against eq 2 and eq 3 typed in the test to 10⁻¹³, with
  `still_distance` and the energy sharing one eq 3; `nonpolar_energy(100 Å²)` = 0.5 kcal/mol;
  every element's ρ and S against values typed from `customgbforces.py` — opened, with line
  numbers — or Bondi; the lone ion from typed ρ; and the spiral's points at the midpoints of n
  equal bands, which an off-by-half spiral (z = 1 − 2k/n) had passed. The Kirkwood gap is
  relabelled a documented measurement — it calls no crate code — and the pair integral gains the
  case `s > r` with `ρ̃ > s − r`. An unverified sentence that AMBER's integral has the inside term
  is removed. 140 tests in the crate.
- **`pantometry-forcefield` reads a protein and its ligand into one typed system: T4 lysozyme
  L99A with benzene, PDB 181L — step 2c-1 of the binding energy.** New module `pdb`:
  `System::from_pdb(text, templates, &Selection)` reads a PDB-format entry (fixed columns, the
  simpler of the two formats to read strictly) and gives back one `Component` that `uff::assign`,
  `ForceField::new`, `qeq` and `solvation` take unchanged, with every atom marked `Part::Protein`
  or `Part::Ligand`. Bond orders, aromatic flags and hydrogens come from the dictionary entries of
  the twenty amino acids and the ligand, passed in; `Coordinates::Structure` is new and says so.
  It does not depend on `pantometry-protein`: no domain does on another. **Hydrogens** are carried
  from each template by Horn's closed-form quaternion superposition (deterministic Jacobi), **not
  of the whole residue but of the hydrogen's parent and its heavy neighbours** — a fragment that
  never spans a torsion — because the whole-residue fit misses by up to 1.94 Å (Arg8; mean
  0.89 Å) where a side chain's χ angles are not the dictionary conformer's; the fragments fit to
  0.029 Å on average and 0.106 Å at worst (Trp126's CA). The backbone amide H is rebuilt on the
  external bisector of `C(i−1)–N–CA` at the template's N–H length, and a rotor (methyl, hydroxyl,
  NH₃⁺) is turned in 10° steps to the step whose nearest heavy atom is farthest: **the
  dictionary's own torsion put Lys162's HZ3 1.495 Å from Asp159's OD1**, which the rule clears;
  186 rotor hydrogens moved. **Protonation at pH 7 by a stated table on the templates' own state,
  which is mixed**: every template is a free amino acid (NH₂, COOH), Asp and Glu are drawn neutral
  and Lys, Arg and His charged (+1), so Asp/Glu lose HD2/HE2 (−1), Lys and Arg are kept (+1), the
  termini are NH₃⁺ (an `H3` added) and COO⁻, and **a histidine with either ring N within 3.2 Å
  of an Asp or Glu carboxylate O is a salt bridge, doubly protonated (+1); any other is neutral on
  Nε2**. **T4 lysozyme's one histidine, His31, has Nδ1 2.66 Å from Asp70's Oδ2, so it is +1** —
  the pair of Anderson, Becktel and Dahlquist, "pH-Induced Denaturation of Proteins: A Single Salt
  Bridge Contributes 3–5 kcal/mol to the Free Energy of Folding of T4 Lysozyme", *Biochemistry*
  29, 2403 (1990), title and citation verified, body not read. **The rule was first written the
  other way round** — neutral, with the H on the ring N facing the carboxylate — which made His31
  δ-protonated and the total +8; a carboxylate at 2.66 Å is an ion pair, and the rule now says so.
  `Selection::histidine` sets any of the three states for one residue. **Strict**: more than one model, an
  insertion code, a hydrogen in the file, an atom the template lacks or of another element, a
  residue missing any heavy atom (every name listed), a chain break by number or by a C–N past
  2.0 Å, a residue that is not `SEQRES`'s, an unnamed `HETATM` residue, a ligand other than once,
  an N-terminal proline, a template without the atoms the rules name, and a sulfur with two
  disulfide partners are each a `PdbError` naming it. Disulfides are found by SG–SG ≤ 2.5 Å and are
  divalent `S_3+2`, not the hypervalent sulfur `ForceField::new` refuses. **181L**: 162 of 164
  residues modelled (Asn163 and Leu164 are `REMARK 465`'s, so Lys162 is the C-terminus and its
  `OXT`, the one heavy atom not read, is placed by the same superposition); 136 waters, two
  chlorides and one 2-hydroxyethyl disulfide dropped and counted; no alternate locations;
  2616 atoms (2604 protein, 12 benzene), 161 peptide bonds, no disulfide (C54T, C97A); **formal
  charge +9**: Arg 13 + Lys 13 − Asp 10 − Glu 8 = +8 counted from `SEQRES`, and His31 +1 by the
  rule applied in the test to the file's coordinates. Typed with no
  refusal — `C_3` 509, `C_R` 276, `C_2` 32, `N_R` 182, `N_3` 40, `N_2` 13, `O_2` 196, `O_3` 37,
  `O_R` 6, `S_3+2` 5, `H_` 1320 — and `ForceField::new` builds 3 413 030 pairs. **Two typings the
  rules give, known and left for 2c-2**: a carboxylate is one `O_2` and one `O_3`, not two
  equivalent oxygens, and arginine's guanidinium is `N_3`, `N_3`, `N_2` around a `C_2`, not
  resonant; both follow from the Kekulé structure the dictionary writes. **QEq on the pocket**
  (benzene and the 18 residues within 6 Å: 322 atoms, total +1): converged in 31 solves, none at
  a bound, 8.5 s unoptimised, and unchanged by His31's charge, which is outside it; its 322 atoms,
  +1, empty bound list and His31's absence are asserted. Within 8 Å — 43 residues and benzene,
  719 atoms, +3 — benzene's charges move by at most 0.0023 e (C −0.10 to −0.12 e, H +0.095 to
  +0.146 e), asserted ≤ 0.004 e by a test ignored by default because it takes 66 s unoptimised
  (9.5 s with `--release -- --ignored`). The 0.004 e first reported was measured before rotors
  were turned. **Open**:
  the hydrogen iteration takes 31 solves on both pockets, against the paper's six to ten; it
  converges, and why it is slower is not yet known. The whole protein is
  not solved here: the N³ scaling of the pocket puts it past an hour unoptimised. **Checked
  against the entry and closed forms, not pinned outputs**: every template heavy atom present at
  the file's coordinates bit for bit; the sequence `SEQRES`'s and the unmodelled residues
  `REMARK 465`'s; one bond per hydrogen; every heavy atom's Kekulé bond orders equal to its
  valence with its formal charge; each residue's formula the textbook residue's adjusted by the
  rule; the total the sequence's; residues − 1 peptide bonds; one N-terminus at residue 1 and one
  C-terminus at 162; the placement rules run on the atoms they name — 158 backbone H (161 residues
  after the first, less three prolines, counted from `SEQRES`) and one N-terminal H3, which makes
  more than 100° with each of the other three bonds (measured 111.7–112.4°); every superposed
  hydrogen at its template's bond length — asserted below 1e-12 relative, measured 1.0e-14, the
  rounding of coordinates near 50 Å — and within `asin(2 √k rmsd / |X − P|)` of its template's angle
  to each heavy neighbour — the bound its own fit earns by the triangle inequality, met at 0.585
  of it at worst (8.65°, Trp126 HA–CA–N); the backbone H coplanar, equidistant in angle and
  external; every fit no worse than a three-atom frame superposition built in the test **and no
  better than distances allow** — `rmsd ≥ max |Δd_ij| / √(2k)`, from `|Δd_ij| ≤ |e_i| + |e_j|`
  and `Σ |e|² = k rmsd²`, held for every whole-residue and fragment fit (the whole-residue fits
  sit at 1.32 or more times it), so a residual reported as zero is caught; every rotor's chosen
  step at least as clear as the dictionary's torsion, put back in the test by its reported turn,
  every turn whole 10° steps (72 rotors, 186 hydrogens), and Lys162's NH₃⁺ named: 1.495 Å from
  Asp159's OD1 at the dictionary's torsion and 1.923 Å turned, asserted > 1.8; and no hydrogen
  within 1.5 Å of a heavy atom. **Measured, not hidden**: one contact below 1.6 Å
  (Thr54 H to Asp47 OD1, 1.540 Å, from an N···O of 2.53 Å in the crystal), 14 below 1.8, 66
  below 2.0, 131 below 2.2. A disulfide is checked on two cysteines built for it. Every refusal is
  fed its input, the two template guards — a histidine without its CE1–NE2 bond, a methionine
  with an H3 already — included. Each of forty sabotages was caught, each restored by copying the
  original back, touched, and its SHA-256 checked. **Eight passed, or were caught only by
  accident, until a test was added for them**: a dropped type's atoms counted once per residue
  (181L's repeated types are one atom each, so a second 8-atom HED was added); the δ histidine's
  Kekulé structure left unmoved (once the rule stopped reaching Nδ1 on 181L, only the overrides do,
  so their valence is now checked); the backbone rule switched off and a rotor keeping its worst
  step, each caught only by the clash test; the N-terminal H3 pointed inward; the whole-residue
  RMSD reported as zero; hydrogen's QEq range narrowed so charges clamp; and the entry's code
  ignored. A review (`unearned-pass-hunter`) found the last six; each now fails a test of its own
  rule, run before and after the fix. Fixtures: the twenty amino acids and `BNZ` from the dictionary
  and 181L as fetched, each with its SHA-256 in `components/README.md`, beside the eight other
  entries of the 1995 papers (182L–188L, 1NHB) recorded and not used. 185 tests in the
  crate, and one ignored.
- **`pantometry-forcefield` computes a binding energy: benzene in T4 lysozyme L99A's rigid
  pocket, and benzene minimised there — step 2c-2.** New module `binding`. `Binding::new(&System,
  cutoff)` cuts a pocket of **whole residues** (any heavy atom within the cutoff of any ligand
  atom, 2c-1's rule, so no residue is split and each fragment's total is an integer) and gives the
  pocket and the ligand **QEq charges separately**, each at its own formal charge: fixed-charge
  practice, which leaves out polarisation and charge transfer between them, and keeps the charges
  exact under any rigid motion of the ligand. Every atom keeps the whole system's UFF type.
  `interaction()` returns `ΔE_bind = E(complex) − E(pocket) − E(ligand)` both as the direct
  protein–ligand cross sums and as the term-by-term difference of three evaluations;
  `desolvation(points)` returns OBC II's complex − pocket − ligand and the buried area;
  `minimise_ligand` relaxes the ligand with every protein atom frozen, on the complex's force field
  less its protein-only terms (a constant), so the force on a ligand atom is the complex's bit for
  bit; `ligand_at(&RigidMotion)`, `with_ligand_positions` and `moved` share the force fields behind
  an `Arc`, so a path for 2c-3 copies positions and nothing else. **`Minimiser::with_frozen`** is new.
  A frozen atom's gradient is zeroed, and its position in every trial is copied rather than rebuilt
  from the minimiser's ångström vector. **Measured: `(x / 1e-10) · 1e-10` is not `x` for 6–11% of
  values**, only where `x / Å` sits near the bottom of its binade, and it is `x` for every
  three-decimal PDB coordinate, so on 181L a sabotage removing the copy passed. The test that
  catches it now places aspirin's frozen coordinates where the round trip fails. A guard zeroing a
  frozen atom's search direction was removed: no input reaches it, because the two-loop recursion
  only scales and adds components that are already exactly zero. **Checked against identities and
  derived bounds, not against experiment**:
  - the vacuum ΔE_bind is exactly the cross terms: the four valence fields of the difference are
    zero, and van der Waals, Coulomb and the total agree with the cross sums. The allowance is pure
    summation rounding, because every term is the same bits in the complex as in its fragment;
    measured worst 0.005 of it, 4.5e-13 kcal/mol. The cross sum also matches one built in the test
    over all 310 × 12 pairs with no exclusion list. **This holds the bookkeeping and the pair
    coverage, not the pair formula**: the test's own sum uses the crate's `Pair::energy` and
    `coulomb`, which `the_energy_terms_against_closed_forms.rs` holds;
  - far away the energy vanishes: at 10³ and 10⁴ Å, van der Waals under `3 D (x/(R−a))⁶` and
    Coulomb under a first-order multipole bound with the Hessian-of-1/r remainder. The measured
    values are 4.5e-7 and 2.7e-9 kcal/mol against bounds of 2.1e-4 and 4.0e-7. **The Coulomb
    bound is about 470× loose**: it bounds the pocket's side by Σ|q| over 310 charges, with no
    cancellation, while the real tail is the pocket's +1 against benzene's 8.5e-4 e Å crystal
    dipole. Expanding the pocket's side as well brings a remainder in the pocket's ~12 Å extent,
    which is looser at these distances, so no tighter bound is claimed. A Coulomb energy wrong by
    1% passes this test and is caught only by the identity above, when the error is in one place
    and not the other. The polar desolvation is under a bound from the screened cross
    terms plus the mean value theorem over the box of Born radii (measured −4.4e-7 against
    2.2e-4), and the buried area is exactly zero;
  - a change of frame (2.1 rad and 13.6 Å) moves ΔE_bind by 5.4e-14 kcal/mol, inside a
    coordinate-rounding bound;
  - frozen protein atoms do not move, to the bit, while benzene moves and converges;
  - the minimised ligand's force from the whole complex's force field is ≤ 1e-4 kcal mol⁻¹ Å⁻¹;
  - the interaction force is −∇ΔE_bind against central differences, with the truncation and
    rounding tolerance of `forces_are_the_gradient.rs`. The worst error is 3.5e-5 of a tolerance,
    and every atom's force is at least 847 tolerances;
  - binding grows no atom's surface area and shrinks no Born radius, both exact (OBC II's `R(I)` is
    increasing, its cubic's discriminant negative). The radius check calls
    `GeneralizedBorn::born_radii` on the fragments directly, so it holds the solvation module on
    this system and says nothing about `Binding::desolvation`;
  - **added after review** (`numerics-reviewer`, fifteen sabotages, three classes passing):
    - each fragment's charges are bit for bit what `Qeq::equilibrate` gives on that fragment's
      own atoms read from the `System`. Before this, taking the ligand's charges from the
      pocket's first twelve atoms passed every test, even though it moved benzene's solvation from
      −2.27 to −18.89 kcal/mol;
    - each charge sum is checked against its formal charge, with an allowance from Higham's
      backward-error bound for the elimination (Theorem 9.4, growth factor taken as 8) in place of
      an unsourced 1e-9. Measured 8.7e-15 e against 1.3e-10;
    - `RigidMotion` is checked against closed forms: a quarter turn about z centred at (1, 2, 3)
      takes (2, 2, 3) to (1, 3, 3), oblique turns match Rodrigues' vector formula written in the
      test, and `ligand_at` leaves the centroid within `(2n + 8) ε X`. Before this, a transposed
      matrix and a rotation about the origin both passed the whole crate;
    - ΔG_GB of the complex, the pocket and the ligand is rebuilt in the test from `born_radii` and
      eq 2 and eq 3, at the crystal and the minimised pose: agreement 1.4e-12, 2.0e-12 and
      2.6e-15 kcal/mol, within `2 (n + 10) ε Σ|t|`. Every atom's area is rerun through
      `surface_area`, and ΔSASA is summed over all atoms, bit for bit. Before this, the polar term
      halved, ΔSASA summed over the ligand's atoms only, and the ligand-alone GB taken at the
      crystal pose (0.059 kcal/mol wrong at the minimised one) all passed;
    - the polar desolvation's change on a change of frame is asserted, no longer only printed: it
      must lie within `Σ |F_GB| √3 · 10 ε X` plus each evaluation's rounding, and measured 1.2e-8
      of that. A frame-dependent descreening distance is caught by this check alone.

  Fourteen sabotages, each restored by copying the original back, touching it and checking its
  SHA-256: thirteen are caught — the pocket rule counting protein hydrogens, a dropped ligand
  atom's pairs, the pocket's charges scaled 0.999, stale positions in the cross sum, a
  non-orthogonal rotation, the frozen copy, the frozen gradient, the reduced force field dropping
  cross pairs, the Coulomb force at 0.99, surface radii enlarged in the complex, the ligand
  fragment solvated at ε = 4, the complex GB taken as its fragments' sum, and the descreening's
  log term at 1.05 — and the fourteenth targeted the unreachable direction guard, which was
  removed. **Three had passed the first version**: the frozen copy, that guard, and the
  ε = 4 fragment, which only the far-field GB bound sees. Nine more were run after review —
  the reviewer's six that had passed (each run again here before and after its fix: passing
  before, caught after), the ligand charged at the pocket's total, a complex GB at stale
  positions, and the frame-dependent descreening — and all nine are caught. **Benzene, reported and not asserted**
  (kcal/mol):

  | cutoff | residues | atoms | pocket charge | vdW | elec | ΔE vacuum | polar | ΔSASA Å² | nonpolar | ΔE + solv | minimised RMSD Å | ΔE minimised |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 6 Å | 18 | 322 | +1 | −11.295 | −1.133 | −12.428 | +7.322 | −294.6 | −1.473 | −6.579 | 0.809 | −22.346 |
  | 8 Å | 43 | 719 | +3 | −13.048 | −1.104 | −14.153 | +10.607 | −294.6 | −1.473 | −5.019 | 0.814 | −24.149 |
  | 10 Å | 60 | 987 | +6 | −13.327 | −1.104 | −14.431 | +12.211 | −294.6 | −1.473 | −3.693 | 0.815 | −24.440 |
  | whole | 162 | 2616 | +9 | −13.500 | −1.104 | −14.604 | +14.962 | −294.6 | −1.473 | −1.115 | 0.816 | −24.623 |

  **The vacuum terms converge** (van der Waals −13.33 at 10 Å against −13.50 whole; Coulomb −1.10
  from 8 Å) and **the polar term does not**: it rises by 7.6 kcal/mol from 6 Å to the whole
  protein. It is also not the small term an apolar ligand in an apolar cavity was expected to give.
  At 6 Å it splits exactly, since GB is quadratic in the charges: +2.27 is the pocket desolvated by
  benzene's volume, +1.94 is benzene's own (of the 2.27 it has alone), and +3.11 is the screened
  cross terms. GB solvates the empty apo cavity as though it were water. **QEq on the whole complex
  instead** moves 0.010–0.012 e onto benzene and changes no benzene charge by more than 0.048 e.
  The cross Coulomb becomes −2.22 instead of −1.13 at 6 Å, but −1.18 instead of −1.10 for the
  whole protein, so most of the 6 Å change comes from the cut. **Minimised in the rigid pocket**,
  benzene slides 0.49 Å and turns (0.81 Å RMSD), and van der Waals goes from −11.30 to −21.28.
  Benzene's own energy goes from 15.89 to 12.14, against 11.88 at its vacuum minimum, so 0.26 of
  strain. The crystal pose has a 1.89 Å contact between a benzene H and Val111's placed HG13,
  which 2c-1's rotor rule, scoring heavy atoms only, cannot see. Freeing every protein hydrogen
  as well moves benzene 0.83 Å, so the move is UFF's own and not the placement's. **Against
  experiment**: ΔG = −5.19 ± 0.16 kcal/mol for benzene, as Mobley et al. tabulate it (*J. Mol.
  Biol.* 371, 1118 (2007), Table 1, measured at 302 K; read at PMC2104542), from Morton, Baase and
  Matthews, *Biochemistry* 34, 8564 (1995), whose own page was not read. The 8 Å row's −5.02 is
  near that by coincidence: the column runs from −6.58 to −1.12 with the cutoff, and is an energy
  at one rigid geometry with no entropy, no protein flexibility, no ligand strain against
  solution, QEq charges with which OBC II over-solvates small molecules, and a polar term that has
  not converged. It cannot be compared with a free energy, and no test does. 2c-1's carboxylate
  typing (`O_2` with `O_3`) cannot reach ΔE_bind: UFF gives both 3.500 Å and 0.060 kcal/mol. QEq
  on the whole protein takes 195 s with `--release`. Ignored by default: the cutoff table, the
  whole protein and the hydrogens' share. The default tests take 13 s unoptimised. 197 tests in the crate, and four ignored.
- **`pantometry-forcefield` on the congener series: nine ligands in T4 lysozyme L99A, set beside
  experiment, and the rigid-pocket energy does not rank them — step 2c-3a.** No library code
  changed. A new test file, `the_congener_series.rs`, builds 182L benzofuran (`BZF`), 183L indene
  (`DEN`), 184L isobutylbenzene (`I4B`), 185L indole (`IND`), 186L n-butylbenzene (`N4B`), 187L
  p-xylene (`PXY`), 188L o-xylene (`OXE`) and 1NHB ethylbenzene (`PYJ`) with 2c-1's pipeline and
  2c-2's `Binding`, beside 181L benzene. The eight entries and seven new dictionary entries are
  fetched from files.rcsb.org, each fetched twice and identical, with SHA-256 in
  `components/README.md`, and embedded with `include_str!`. **The entries, checked by string
  operations on the files**: each `SEQRES` is identical to 181L's, and each has exactly 181L's
  three `SEQADV` conflicts against P00720 (C54T, C97A, L99A). Each has `REMARK 465` Asn163 and
  Leu164, no alternate location, and `HETATM` residues that are waters (116–134), `CL` ×2, one
  `HED` and the ligand once with every heavy atom. The drop list is the same for all nine. None was
  refused, none has a missing atom, and each builds to 181L's 162 residues and 2604 protein atoms,
  with His31 +1 by the rule (Nδ1/Nε2 to Asp70 2.64–2.71 Å) and a total of +9. All of this is
  asserted. Reported and not asserted: helix F moves 2.27–2.62 Å at Ala112 (Cα, no superposition)
  for indene, isobutylbenzene and o-xylene, and no Cα moves more than 0.72 Å in the other six.
  **Each ligand is the molecule its entry names**: a formula cannot tell o- from p-xylene or n- from
  isobutylbenzene, so each is held to its bond graph, written by hand, plus the xylenes'
  methyl–methyl bond count. The bond graph is every heavy atom's element, aromatic flag, heavy
  neighbours and hydrogens. **Experiment**: ΔG° from Mobley et al., *J. Mol. Biol.* 371, 1118
  (2007), Table 1, read at PMC2104542. Its caption gives the values as from Morton, Baase and
  Matthews, *Biochemistry* 34, 8564 (1995), by ITC at 302 K, and that paper itself was not read.
  Benzofuran −5.46 ± 0.03, benzene −5.19 ± 0.16, ethylbenzene −5.76 ± 0.07, indene −5.13 ± 0.01,
  indole −4.89 ± 0.06, isobutylbenzene −6.51 ± 0.06, n-butylbenzene −6.70 ± 0.02, o-xylene
  −4.60 ± 0.06, p-xylene −4.67 ± 0.06 kcal/mol. The table's toluene and n-propylbenzene have no
  structure among the 1995 entries and are not used. **The crystal pose's vacuum ΔE measures the
  hydrogen placement for three of the nine.** Both partners' hydrogens are placed: the ligand's by
  its template, the protein's by 2c-1's rotor rule, which clears heavy atoms only. Neither
  placement sees the other partner's hydrogens. Indene's H12 sits 1.41 Å from Val111's HG13, so
  van der Waals is **+224.4 kcal/mol**. N-butylbenzene (H3′2–Leu118 HD22 1.65 Å) is +80.8, and
  ethylbenzene (HCD1–Val111 HG13 1.64 Å) is +12.8. So van der Waals < 0 is not asserted at the
  crystal pose. It is asserted at two other poses: with **every hydrogen relaxed** on the crystal's
  heavy atoms (complex force field, heavy atoms frozen, closest H–H then 1.94–2.13 Å), and after
  minimisation. Also asserted: ΔSASA < 0 at both poses, every protein atom unmoved to the bit, every
  minimisation converged, and **181L's 6 and 8 Å numbers unchanged from 2c-2's table**, to half a
  unit in its last printed digit. 6 Å, kcal/mol:

  | entry | ligand | ΔG°exp | vdW | elec | ΔE vac | polar | nonpolar | total | ΔE vac, H relaxed | RMSD Å | vdW min | ΔE vac min | polar min | total min |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 181L | benzene | −5.19 | −11.29 | −1.13 | −12.43 | +7.32 | −1.47 | −6.58 | −20.42 | 0.809 | −21.28 | −22.35 | +7.33 | −16.50 |
  | 182L | benzofuran | −5.46 | −14.04 | −2.65 | −16.69 | +10.98 | −1.79 | −7.50 | −27.77 | 0.566 | −24.81 | −27.65 | +11.46 | −17.99 |
  | 183L | indene | −5.13 | +224.44 | −1.49 | +222.95 | +10.34 | −1.88 | +231.41 | −18.62 | 0.749 | −22.81 | −24.18 | +10.58 | −15.48 |
  | 184L | isobutylbenzene | −6.51 | −6.50 | −0.66 | −7.16 | +12.88 | −2.20 | +3.52 | −35.80 | 0.667 | −31.67 | −33.16 | +19.27 | −16.11 |
  | 185L | indole | −4.89 | −20.31 | −1.90 | −22.21 | +10.26 | −1.78 | −13.72 | −25.43 | 0.148 | −24.08 | −26.54 | +11.05 | −17.26 |
  | 186L | n-butylbenzene | −6.70 | +80.80 | −1.01 | +79.79 | +16.22 | −2.17 | +93.83 | −14.99 | 0.551 | −29.16 | −30.02 | +18.89 | −13.36 |
  | 187L | p-xylene | −4.67 | −14.00 | −0.02 | −14.01 | +9.61 | −1.91 | −6.31 | −25.46 | 0.483 | −29.13 | −29.64 | +10.43 | −21.14 |
  | 188L | o-xylene | −4.60 | −18.39 | −0.73 | −19.12 | +9.59 | −1.91 | −11.44 | −24.81 | 0.446 | −28.08 | −29.14 | +10.15 | −20.91 |
  | 1NHB | ethylbenzene | −5.76 | +12.77 | −0.83 | +11.95 | +11.81 | −1.86 | +21.89 | −24.84 | 0.370 | −26.37 | −27.51 | +12.65 | −16.74 |

  The 6 Å pockets are 18–26 residues and 322–439 atoms, with pocket charges 0 to +2. At 8 Å they
  are 40–49 residues, 671–804 atoms and +3 each; the full 8 Å table is printed by the test.
  **Against ΔG°exp, n = 9, 95% Fisher intervals** (Spearman's with the 1.06 factor usually
  attributed to Fieller, Hartley and Pearson 1957, quoted from memory, not read):

  | column | 6 Å Pearson r | Spearman ρ | slope | RMS after offset | 8 Å Pearson r | Spearman ρ | slope | RMS after offset |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | ΔE vac, hydrogens relaxed | +0.04 [−0.64, +0.69] | −0.03 [−0.69, +0.66] | +0.35 | 5.61 | +0.08 [−0.62, +0.71] | −0.02 [−0.69, +0.67] | +0.64 | 5.60 |
  | ΔE vac, minimised | +0.46 [−0.29, +0.86] | +0.32 [−0.46, +0.82] | +1.98 | 2.79 | +0.49 [−0.25, +0.87] | +0.33 [−0.44, +0.82] | +2.20 | 2.91 |
  | total, minimised | −0.77 [−0.95, −0.22] | −0.72 [−0.94, −0.08] | −2.53 | 2.93 | −0.90 [−0.98, −0.58] | −0.88 [−0.98, −0.51] | −3.97 | 3.81 |

  (The crystal-pose columns are printed too, and they are the clashes: r = −0.14, spread 245.)
  **What this can and cannot mean.** The nine ligands are of one size class, and ΔG°exp spans 2.10
  kcal/mol, while every computed column spans 8–21. With n = 9, a Pearson r is different from zero
  at the two-sided 5% level only past 0.666. The vacuum columns' intervals all contain zero, so the
  minimised vacuum ΔE shows no correlation that nine points can distinguish from none. The solvated
  total's interval excludes zero **on the wrong side**: it ranks the series roughly backwards. The
  cause is in the table. The polar term is largest for the two largest ligands (+19.3 and +18.9
  minimised, against +7.3 for benzene), which also bind best. It grows from 6 to 8 Å for every
  ligand (isobutylbenzene +19.3 to +26.7), as 2c-2 measured for benzene, and with QEq charges OBC II
  over-solvates. So the total's anticorrelation is a measurement of the polar term's size
  dependence, which has not converged, and it is not a measurement of binding. **None of these
  statistics is asserted.** The statistics functions are held to closed forms instead: `y = x²` on
  1…5 gives r = 60/√3740, slope 6, ρ = 1 and offset-removed RMS √52.8, plus a reversal, tied ranks
  and Fisher's interval at r = ½, n = 12. Ignored by default: the series, 100 s with `--release
  -- --ignored`. The three default tests take 1.1 s unoptimised. Each new assertion was sabotaged:
  26 sabotages, each restored by copying the original back, touching it and checking its
  SHA-256, and each caught by the assertion it targeted. The 19 on the default tests: o-xylene's
  methyl moved to meta (caught by the methyl count alone), isobutylbenzene rewired to n-butyl with
  formula and valence intact (caught by the graph alone), a `SEQRES` residue, a `SEQADV` line, a
  water renamed, the His rule at 2.6 Å, dropped residues double-counted, a ligand H removed,
  lysine neutral, serine's HG removed, residue 1 deleted, the `HEADER` code, a ligand atom deleted,
  an alternate location added, ties not averaged, Fisher at 2σ, Pearson ×0.999, the RMS without
  its offset, and the slope's denominator. The 7 on the series: one pocket atom unfrozen, ΔSASA's
  sign, the cross van der Waals' sign, hydrogens not relaxed, the hydrogen relaxation cut to 5
  steps, the minimiser cut to 20 steps, and the surface tension at 0.0051. The last is caught only by 181L's committed numbers. 200 tests in
  the crate, and five ignored. *Superseded by step A-1 below: the "hydrogens relaxed" pose was a
  test helper outside `Binding`, with no solvation and the fragments taken at the complex's
  positions. `Binding::relaxing_hydrogens` replaces it, and the two sabotages of that helper went
  with it.*
- **`pantometry-forcefield` relaxes each system's own hydrogens as part of the binding
  calculation, and the three crystal-pose hydrogen clashes are gone — step A-1.** **What was
  wrong**: every hydrogen in a `Binding` is placed, the protein's by 2c-1's rotor rule (which
  scores heavy atoms only) and the ligand's by its template, and neither placement sees the other
  partner's. At the crystal pose three of the nine T4 lysozyme L99A complexes clash hydrogen to
  hydrogen: indene (183L) H12–Val111 HG13 at 1.41 Å, van der Waals **+224.4 kcal/mol**;
  n-butylbenzene (186L) 1.65 Å, +80.8; ethylbenzene (1NHB) 1.64 Å, +12.8. 2c-3a worked around
  this in a test helper outside `Binding`. That helper relaxed every hydrogen on the complex's
  force field and took the fragments at the complex's positions. A `Binding` could not move a
  pocket atom, so that pose had no solvation numbers. Its fragments also kept whatever their
  hydrogens had bent to make room for the partner, strain the apo pocket would not have, and
  counted it as binding. **New**: `Binding::relaxing_hydrogens(max_steps, tolerance)` and
  `relaxing_hydrogens_of(mask, …)`. Each relaxes, in each of the three systems, that system's own
  hydrogens with every heavy atom frozen (`Minimiser::with_frozen`): first the complex, then the
  pocket alone and the ligand alone, each from where the complex left it. Then `ΔE_bind =
  E(complex, H relaxed) − E(pocket, its own H relaxed) − E(ligand, its own H relaxed)`, and the OBC
  II solvation and ΔSASA are taken at each system's own positions. `Binding` holds the fragments'
  positions beside the complex's (`pocket_alone_positions`, `ligand_alone_positions`).
  `Interaction::reorganisation` is new: `[E(pocket) at the complex's positions − at its own] +
  [the same for the ligand]`, term by term. `Interaction::total()` adds it. Also new:
  `HydrogenRelaxation` (where each of the three stopped), `free_hydrogens()` and
  `Binding::HYDROGEN_TOLERANCE`. **An unrelaxed binding is unchanged to the bit**: its
  reorganisation is exactly zero and not computed, its fragments' positions are the complex's
  slices, and 2c-2's committed 181L numbers are still asserted by the series. `minimise_ligand`
  on a relaxed binding moves the pocket's free hydrogens with the ligand (heavy atoms still
  frozen), relaxes both fragments again, and stays relaxed. `ligand_at` returns an unrelaxed binding, since its
  hydrogens were relaxed against the old pose. `moved` carries the fragments with it. **The
  choices**, each stated in the `binding` module documentation:
  - **Each fragment at its own minimum** is the consistent definition, because ΔE_bind is a
    reaction energy and each side is taken at its own minimum.
  - **Each fragment starts from the complex's relaxed hydrogens, not from the placement.**
    Started apart, two minimisers on the same flat rotor landscape can stop in different minima
    for reasons unrelated to binding. Started from the complex, a fragment moves only where the
    partner's removal pushes it. The reorganisation is then never negative, and it vanishes
    exactly when the partner exerts no force. This is not a search: a fragment whose best
    arrangement lies in another basin is not found there.
  - **Pocket hydrogens beyond the cutoff are held — a truncation, with a measured cost.** A pocket
    hydrogen is free only when its heavy parent is within the cutoff of a ligand atom, which is
    the pocket's own rule applied to atoms; on benzene that frees 106 of the pocket's 167. The
    held ones do feel the ligand: at the relaxed complex the partner's force on them is 0.001–0.14
    kcal mol⁻¹ Å⁻¹ on four entries measured, nearly all of it above the tolerance. Freeing all of
    them moves benzene's relaxed ΔE by 0.25 kcal/mol at 6 Å and its reorganisation by 1.4 (+2.62 to
    +4.03). On the other eight the relaxed ΔE moves by at most 0.10 at 6 Å and 0.03 at 8 Å, and
    benzene's relaxed-and-minimised ΔE by 0.016. The rule is kept because it is the pocket's own
    and because the outer hydrogens relax against the cut as well as against the ligand, not
    because they carry no information.
  - **The tolerance is 2e-3 kcal mol⁻¹ Å⁻¹, not 2c-2's 1e-4, and the reason is empirical.** At
    1e-4 the minimiser stalled four times on the series, at 1.57, 1.61, 1.74 and 4.27e-4
    kcal mol⁻¹ Å⁻¹, and 2e-3 converges on all 18 entry-cutoff pairs: every relaxation and every
    minimisation from a relaxed complex, every-hydrogen-free comparison included. **Against 1e-4
    on benzene at 6 Å, where 1e-4 converges, a test asserts that the looser tolerance moves no
    printed digit**: ΔE_bind and the reorganisation within 5e-3 kcal/mol (measured 6e-7 relaxed,
    2.8e-5 minimised), the minimised ligand within 5e-4 Å RMS (measured 1.7e-4). Each system's
    two energies also lie within a convexity bound computed from their own forces and
    displacements. The series prints the same comparison for all nine at 6 Å: the relaxed ΔE
    within 6e-7 on eight and 0.040 on p-xylene, whose two runs end in neighbouring minima; the
    minimised ΔE within 3.2e-4; the minimised ligands within 4.2e-4 Å RMS. **Why a pocket stalls,
    as a scale and not a bound**: these energies are sums of 0.5–3.8·10⁵ terms whose totals reach
    −39,600 kcal/mol, and a 1e-9 Å nudge of one free atom moves the sum by a rounding noise δE that
    varies from nudge to nudge. Over 300 nudges of benzene's 6 Å pocket the median is 3.6e-12
    kcal/mol and the largest 2.1e-11, about six times the first draw; single draws on the 8 Å
    pockets reach 1.8e-10. The last steps along an N–H or O–H stretch (k ≈ 1190 kcal mol⁻¹ Å⁻²,
    measured) lower the energy by only `g²/2k`, so `√(4kδE)`, between 1.3e-4 and about 2e-3 over
    that spread, is the order of the force below which the decrease can drown in the noise. That
    is where the stalls were. It does not predict where a given run stops: 1e-4 converges on
    benzene at 6 Å, where the formula gives 3.2e-4. **A stall is now a failure everywhere in the
    series**, since accepting one accepted it at any force. The first version allowed the
    every-hydrogen-free comparison to stall, and with its tolerance sabotaged to 1e-6 every one of
    those relaxations stalled and the series passed.
  - The charges stay QEq's at the crystal's positions, as everywhere in `binding`.

  **Checked by identities and exact facts**, in `benzene_in_its_pocket.rs` (eight new default
  tests and one ignored), and on all nine entries in the series:
  - the mask is the rule, recomputed from the system's bonds. Every atom it holds — every heavy
    atom, and every pocket hydrogen beyond the cutoff — is at its crystal position bit for bit in
    the complex, the pocket alone and the ligand alone. The ligand alone's heavy atoms are the
    minimised complex's bit for bit. In each system some hydrogen moved from the crystal (112, 106
    and 6 on benzene), and in each fragment some hydrogen moved from where the complex left it
    (106 and 6). Each fragment's relaxation took steps at the crystal pose;
  - each system's free atoms are converged on that system's **whole** force field (not the reduced
    one the minimiser used), heavy atoms' forces excluded: the largest is 1.8e-3 against 2e-3. Each
    relaxation had work to do, with a free atom at 113, 21 and 13 kcal mol⁻¹ Å⁻¹ at its start;
  - each relaxation lowers its own system's energy, to twice the two evaluations' summation
    rounding: the complex by 74.6 kcal/mol, the pocket by 1.06 and benzene by 1.56;
  - **ΔE_relaxed = the cross terms at the complex's positions + the pocket's and the ligand's
    reorganisation**, field by field. Every part is evaluated in the test (the cross sums over all
    310 × 12 pairs, the five energies by the fragments' force fields), to the rounding of the
    sums. The difference agrees to 2.7e-13 kcal/mol relaxed and 3.4e-12 relaxed and minimised. The
    crate's `reorganisation` is the two brackets. The identity holds for any fragment positions,
    so it catches an energy taken at the wrong positions or a missing bracket, not fragments left
    unrelaxed;
  - **the desolvation is rebuilt at each system's own positions** on the relaxed binding and the
    one minimised from it, as 2c-2's test does at the crystal and minimised poses: ΔG_GB of each
    system from `born_radii` and eq 2–3, every atom's area, and ΔSASA bit for bit. The fragments'
    positions are asserted to differ from the complex's;
  - **far away the fragments' relaxations are the complex's exactly.** Each fragment starts where
    the force on it is the complex's less the partner's. The test asserts the bound `|F_complex| +
    |F_partner| ≤ tolerance` on every free atom, from which the fragments take no step: their
    positions are the complex's bit for bit and the reorganisation is exactly zero. ΔE_bind, the
    polar term and ΔSASA then lie inside 2c-2's tail bounds. This uses 10⁴ and 10⁵ Å, not 10³: at
    10³ Å the pocket's +1 still pulls a benzene hydrogen with 3.3e-5 kcal mol⁻¹ Å⁻¹, against
    3.3e-7 at 10⁴;
  - a relaxed binding moved whole takes its fragments with it, bit for bit;
  - **relaxed, no ligand hydrogen is within 1.8 Å of a pocket hydrogen**. At that distance UFF's
    H···H pair is +11.2 kcal/mol, which the test computes and asserts repulsive. Measured minima
    are 1.94 Å relaxed (indene) and 2.14 Å relaxed and minimised, against 1.41 at the crystal pose.

  Seventeen sabotages, each restored by copying the original back, touching it and checking its
  SHA-256, and each caught by the check it targeted. Four came from `numerics-reviewer`, and each
  passed every test before its fix and is caught after it: `HYDROGEN_TOLERANCE` at 1e-1 (by the
  tolerance test), the desolvation's fragments put back at the complex's slices (by the rebuilt
  desolvation), both fragment relaxations skipped (caught before only by the convergence test,
  now also by the moved-from-the-complex check and the rebuilt desolvation), and the
  every-hydrogen-free relaxation at 1e-6 (by the series, now that a stall fails). The other
  thirteen: a heavy atom unfrozen; the mask taken from
  the hydrogen's own distance and not its parent's; the pocket's relaxation cut to 5 steps; the
  pocket relaxed into a copy and reported converged, which only the whole-force-field force check
  sees; a ligand hydrogen displaced 0.1 Å after relaxing; the ligand alone's energy at the
  complex's positions; the reorganisation without the ligand's bracket; the pocket relaxed from
  the placement and not the complex, which only the far-field test sees (106 steps); `ligand_at`
  keeping the relaxation; `moved` leaving the pocket behind; `minimise_ligand` not relaxing the
  fragments again; an unrelaxed ligand alone at the crystal pose; and the complex left unrelaxed,
  caught by the H–H threshold on 183L at 1.41 Å. **Before and after, 6 Å, kcal/mol**:

  | entry | closest H–H Å, crystal → relaxed → relaxed and minimised | vdW crystal → relaxed | ΔE vac crystal | relaxed: ΔE vac (reorganisation, ligand's share) | total | RMSD Å from crystal H → from relaxed | ΔE vac minimised from crystal H → from relaxed | total minimised from relaxed |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 181L benzene | 1.89 → 2.07 → 2.33 | −11.29 → −17.17 | −12.43 | −16.14 (+2.62, +1.56) | −9.25 | 0.809 → 0.902 | −22.35 → −22.55 | −15.86 |
  | 182L benzofuran | 1.85 → 2.03 → 2.29 | −14.04 → −24.23 | −16.69 | −25.18 (+2.68, +1.34) | −14.59 | 0.566 → 0.511 | −27.65 → −29.38 | −18.53 |
  | 183L indene | **1.41** → 1.94 → 2.14 | **+224.44** → −16.66 | +222.95 | −7.26 (+11.38, +4.22) | +2.68 | 0.749 → 0.629 | −24.18 → −24.25 | −14.11 |
  | 184L isobutylbenzene | 1.68 → 2.12 → 2.18 | −6.50 → −34.40 | −7.16 | −33.20 (+2.56, +0.79) | −17.13 | 0.667 → 0.643 | −33.16 → −32.06 | −13.05 |
  | 185L indole | 1.97 → 2.05 → 2.16 | −20.31 → −22.89 | −22.21 | −23.30 (+2.15, +0.73) | −13.23 | 0.148 → 0.209 | −26.54 → −25.67 | −15.30 |
  | 186L n-butylbenzene | **1.65** → 1.95 → 2.21 | **+80.80** → −13.50 | +79.79 | +8.80 (+23.80, +12.50) | +26.12 | 0.551 → 0.509 | −30.02 → −29.16 | −10.23 |
  | 187L p-xylene | 1.88 → 2.03 → 2.27 | −14.00 → −25.06 | −14.01 | −23.12 (+2.33, +1.08) | −13.56 | 0.483 → 0.443 | −29.64 → −28.31 | −18.56 |
  | 188L o-xylene | 1.87 → 2.06 → 2.19 | −18.39 → −24.02 | −19.12 | −21.91 (+2.92, +1.19) | −12.75 | 0.446 → 0.398 | −29.14 → −28.63 | −19.28 |
  | 1NHB ethylbenzene | **1.64** → 2.12 → 2.19 | **+12.77** → −23.55 | +11.95 | −16.00 (+8.85, +3.20) | −4.06 | 0.370 → 0.403 | −27.51 → −26.76 | −14.70 |

  **The clash is gone and van der Waals is attractive in every entry, but on the crystal's heavy
  atoms two complexes still pay for it.** Indene's relaxed hydrogens cost +11.4 kcal/mol of
  reorganisation and n-butylbenzene's +23.8, of which +12.5 is the ligand's own. The crystal
  heavy-atom pose is tighter than UFF's 2.886 Å hydrogen allows, and minimising the ligand
  releases it (+1.5 and +5.4). The cross terms alone, 2c-3a's "hydrogens relaxed" column, were
  −18.62 and −14.99. That column is reproduced exactly: benzene's cross terms with every hydrogen
  free are −20.42 here as there. What it left out was this reorganisation. **The pose moves
  little**: −0.12 to +0.09 Å RMSD against the pose minimised from the placed hydrogens. The
  minimised vacuum ΔE moves by −1.73 to +1.33 kcal/mol, partly because it now includes the
  reorganisation that remains (+0.57 to +5.36). **Against ΔG°exp, n = 9, 95% Fisher intervals**
  (computed as 2c-3a did; nothing asserted):

  | column | 6 Å Pearson r | Spearman ρ | slope | RMS after offset | 8 Å Pearson r | Spearman ρ | slope | RMS after offset |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | ΔE vac, relaxed | −0.33 [−0.82, +0.43] | −0.15 [−0.75, +0.59] | −5.30 | 11.76 | −0.32 [−0.81, +0.44] | −0.15 [−0.75, +0.59] | −5.03 | 11.65 |
  | total, relaxed | −0.52 [−0.88, +0.22] | −0.17 [−0.76, +0.58] | −9.27 | 13.18 | −0.56 [−0.89, +0.17] | −0.18 [−0.77, +0.56] | −10.50 | 13.84 |
  | ΔE vac, relaxed and minimised | +0.48 [−0.27, +0.87] | +0.40 [−0.38, +0.85] | +1.84 | 2.49 | +0.51 [−0.23, +0.88] | +0.45 [−0.33, +0.86] | +2.03 | 2.55 |
  | total, relaxed and minimised | −0.81 [−0.96, −0.33] | −0.80 [−0.96, −0.27] | −3.16 | 3.38 | −0.87 [−0.97, −0.50] | −0.70 [−0.93, −0.04] | −4.90 | 4.65 |

  Beside 2c-3a's minimised columns (+0.46 and −0.77 at 6 Å, unchanged here), **relaxing the
  hydrogens does not make the series rank.** The vacuum intervals still contain zero. With n = 9,
  r must pass 0.666 to differ from zero at the two-sided 5% level. The solvated total still
  correlates on the wrong side, for 2c-3a's reason: the polar term grows with the ligand and with
  the cutoff. At the crystal's heavy atoms the relaxed columns anticorrelate weakly. The driver is
  n-butylbenzene: it binds best, and its +23.8 reorganisation makes its relaxed ΔE the worst,
  +8.80. The series test is updated to the relaxed `Binding`, with the every-hydrogen-free
  comparison and the ligand's share of the reorganisation as columns. It takes 112 s with
  `--release -- --ignored`. The benzene file's default tests take 13.4–13.8 s unoptimised, against
  13 s before. 208 tests in the crate, and six ignored.
- **`pantometry-forcefield`: why the polar desolvation does not converge with the pocket, an
  empty-cavity reference for it, and what is left — step A-2.** **What was not known**: 2c-2
  measured OBC II's polar desolvation of benzene in T4 lysozyme L99A growing from +7.32 kcal/mol at
  6 Å to +14.96 for the whole protein, and 2c-3a found it largest for the largest ligands, which
  turned the solvated total against experiment (r ≈ −0.8). 2c-2 split benzene's +7.32 exactly into
  the pocket desolvated by benzene's volume (+2.27), benzene's own (+1.94) and the screened cross
  terms (+3.11), and said GB solvates the empty apo cavity as water. **That is true, and it is not
  what grows most.** Split the same way at the crystal pose, kcal/mol:

  | | 6 Å | 8 Å | 10 Å | 15 Å | 20 Å | whole |
  | --- | --- | --- | --- | --- | --- | --- |
  | benzene: polar | +7.32 | +10.61 | +12.21 | +14.17 | +14.86 | +14.96 |
  | the pocket by benzene's volume | +2.27 | +3.62 | +4.31 | +5.43 | +5.89 | +5.93 |
  | benzene's own | +1.94 | +2.23 | +2.30 | +2.35 | +2.36 | +2.37 |
  | cross | +3.11 | +4.77 | +5.61 | +6.39 | +6.61 | +6.66 |
  | n-butylbenzene: polar | +16.22 | +22.27 | +25.09 | | | +31.01 |
  | the pocket by its volume | +5.11 | +7.23 | +8.30 | | | +11.17 |
  | its own | +3.66 | +4.02 | +4.10 | | | +4.20 |
  | cross | +7.44 | +11.02 | +12.69 | | | +15.65 |

  All three parts grow, the cross terms most, and **the term does converge, slowly: within 0.10
  kcal/mol of the whole protein at 20 Å (1979 atoms), 0.80 at 15 Å.** The diagnosis is a new
  ignored file, `the_polar_term_diagnosed.rs` (751 s with `--release`). It attributes each part to the protein's atoms (an
  atom's Born term and half of each pair term it is in), bins them by distance from the ligand,
  and asserts only that the attribution sums to the part within its computed rounding (an
  a-priori bound, 3.2e-5 kcal/mol at 6 Å against 4.4e-13 measured):
  - **Nothing in the descreening is wrong.** Past 10 Å from the ligand its descreening of each
    protein atom is within 1.2% of the sphere's asymptote `Σ s³/3r⁴` (1.1% for n-butylbenzene).
  - **The far tail is small and falls as it should.** In the whole protein, the atoms more than
    15 Å from benzene carry +0.18 of the pocket part and cross terms' +12.59 (+0.29 of 26.82 for
    n-butylbenzene), and the 789 atoms past 20 Å carry +0.017: a protein atom 20 Å away is not
    desolvated by the ligand. Per atom the net falls 7 and 9.5 times from the 10–15 Å shell to the
    15–20 Å one, about as `r⁻⁶`, a dipole's field squared, because GB's pair terms cancel most of
    the incoherent `Σ q_i² ΔI_i`: Born terms +42.78 against pairs −36.85 for benzene.
  - **OBC II amplifies a buried atom's response**: `d(1/R)/dI` is 1.06–1.88 times HCT's, median
    1.78, and that linearisation gives the actual `Σ q² Δ(1/R)` past 10 Å to 5e-5 of itself
    (−2.142e-2 e² Å⁻¹, against HCT's −1.261e-2). That is OBC II's fit, not an error.
  - **The growth is the cut.** A cut pocket solvates its own cut surface: benzene's 6 Å pocket's
    charges inside the whole protein's dielectric (every other atom present, uncharged) give +9.47
    against +7.32 in the fragment, and the remaining +5.49 is the charges the cut leaves out, nearly
    all within 15 Å (a 15 Å pocket's charges in the whole dielectric give +14.79). QEq's charges
    moving with the cutoff matter little: the whole protein's on the 6 Å pocket give +7.80
    (n-butylbenzene +16.08 against +16.22).

  No defect was found, so nothing was fixed and there is no test that fails on an old behaviour.

  **Is the cavity empty? Read, not assumed.** **PDB 1L90**, apo L99A (Eriksson, Baase and Matthews,
  *J. Mol. Biol.* 229, 747 (1993), whose `REMARK 1` cites the cavity's report, *Science* 255, 178
  (1992)), was fetched twice, identical, and is committed with its SHA-256. A new default test
  reads it by string operations: the same `SEQRES` and three `SEQADV` conflicts as 181L, an
  isomorphous cell (c 96.8 Å against 97.0), Cα 0.27 Å RMS from 181L's without superposition, **no
  water within 5 Å of any of 181L's benzene carbons** (the nearest is 7.79 Å), and no protein atom
  within 3 Å of their centroid (Ala99 CB at 3.26). **Collins, Hummer, Quillin, Matthews and
  Gruner, *PNAS* 102, 16668 (2005)**, read at PMC1283839: the L99A cavity "is believed to be
  entirely empty under ambient conditions"; by high-pressure crystallography about 0.5 waters
  enter it at 100 MPa and two to four at 200 MPa, and their simulation finds it predominantly
  empty below about 100 MPa. The 1992 and 1993 papers themselves were not read.

  **New: `ApoCavity` and `Binding::desolvation_with_cavity(charges, points, cavity)`.**
  `ApoCavity::Solvent` is GB's own reference and the default: `desolvation` and `desolvation_with`
  are unchanged to the bit, and the series still asserts 181L's committed numbers. `ApoCavity::Empty`
  solvates the pocket alone with the ligand's atoms present at the complex's ligand positions and
  uncharged, so the ligand's volume descreens the apo pocket as it does the complex; the polar term
  is then the ligand's own desolvation plus the cross terms. **An option, not the default**:
  whether a cavity is wet is a fact about the pocket that the force field cannot know. The areas
  and the nonpolar term are the same under both (the cavity's walls still count as
  solvent-accessible, a choice). `Desolvation` gains a `cavity` field saying which reference it
  holds. **Checked by identities**, in `benzene_in_its_pocket.rs`:
  - under both references, with the binding's charges and with either partner's set to zero, the
    complex's and the ligand's ΔG_GB and every area are the default's bit for bit, and the ghost
    pocket term is rebuilt in the test from `born_radii` and eq 2–3 to `2 (n + 10) ε Σ|t|`, at the
    crystal pose and relaxed (where the pocket alone's hydrogens are its own);
  - **with the ligand's charges zero**, the solvent reference's polar term less the empty one's is
    the ghost pocket term less the ghost-free one — bit for bit at the crystal pose, where the
    empty-cavity polar term is exactly 0.0 and its pocket term is the complex's bits;
  - **the decomposition sums**: `polar = pocket part + ligand's part + cross` under each
    reference, the cross terms summed directly in the test, to 1.6e-13 kcal/mol (5.8e-9 of the
    allowance);
  - far away (10³ and 10⁴ Å; relaxed, 10⁴ and 10⁵) the empty-cavity polar term lies inside 2c-2's
    tail bound: −4.4e-7 kcal/mol at 10³ Å.

  Seven sabotages, each restored by copying the original back, touching it and checking its
  SHA-256, and each caught: the ghost at the ligand alone's positions (by the rebuilt ghost term,
  relaxed); the ghost keeping the ligand's charges (also by both far-field tests); no ghost at all;
  the ghost at the crystal pose (both far-field tests and the rebuild); the default switched to
  `Empty` (three tests, the default's own check among them); an uncharged atom descreening nothing
  in `GeneralizedBorn`, which the rebuilt-desolvation tests cannot see because QEq's charges are
  never exactly zero (caught by the new test and by
  `burying_benzene_buries_area_and_costs_polar_solvation`, whose split zeroes them); and a water moved into 1L90's cavity (by the 1L90 test). An eighth,
  the test's own direct cross sum with `R_i R_j` 1% large, confirms that the decomposition
  assertion can fail.

  **Added after review** (`unearned-pass-hunter`; each sabotage was run before the fix and passed,
  and is caught after it, restored by copy, touch and SHA-256):
  - **nothing pinned `#[default]` on `ApoCavity`**, because `desolvation_with` names `Solvent`
    explicitly. Moving `#[default]` to `Empty` passed every test; now
    `ApoCavity::default() == Solvent` is asserted;
  - **the bit-exact checks were conditional on a flag nothing asserted.** They run only where
    the three systems share positions, so a "relaxed" binding that was not relaxed passed them
    vacuously: `relaxed()` returning the crystal binding passed. Now the flag is asserted true
    for the crystal binding and false for the relaxed one;
  - **the 20 Å test passed with every entry skipped**, printing NaN statistics on n = 0, because
    it skipped any error from `Binding::new`; a cutoff of zero passed. Now only a QEq failure is
    a skip and anything else panics (the cutoff sabotage is caught there). The skipped codes
    must equal the expected set, empty today, and at least three entries must be measured. QEq
    capped at one solve, so that every entry fails as a QEq failure, is caught by the
    skipped-set assertion. That sabotage was run only after the fix; by construction the old
    test skipped it;
  - the whole-protein numbers once given here for 182L and 183L came from no committed code and
    are removed; 181L's and 186L's are regenerated by `the_polar_term_diagnosed.rs`;
  - the `binding` documentation said the 1L90 test asserts no water within 7.79 Å; it asserts
    none within 5 Å, and 7.79 is the measured nearest.

  **The series, relaxed and minimised, kcal/mol** (the series test also prints the crystal pose's
  split; nothing is asserted):

  | entry | ΔG°exp | 6 Å: polar | pocket by volume | ligand's own | cross | empty cavity | total | total, empty | 8 Å: polar | empty cavity | total | total, empty |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 181L benzene | −5.19 | +8.17 | +2.77 | +1.96 | +3.45 | +5.41 | −15.86 | −18.62 | +12.04 | +7.52 | −13.59 | −18.11 |
  | 182L 2,3-benzofuran | −5.46 | +12.64 | +4.58 | +4.15 | +3.91 | +8.10 | −18.53 | −23.07 | +15.39 | +8.81 | −17.03 | −23.61 |
  | 183L indene | −5.13 | +12.01 | +5.04 | +2.35 | +4.62 | +6.99 | −14.11 | −19.13 | +16.02 | +8.87 | −11.91 | −19.05 |
  | 184L isobutylbenzene | −6.51 | +21.20 | +6.45 | +5.47 | +9.29 | +14.77 | −13.05 | −19.49 | +29.71 | +20.55 | −6.86 | −16.01 |
  | 185L indole | −4.89 | +12.14 | +4.60 | +3.25 | +4.28 | +7.56 | −15.30 | −19.87 | +16.37 | +9.54 | −12.81 | −19.64 |
  | 186L n-butylbenzene | −6.70 | +21.11 | +6.52 | +4.68 | +9.91 | +14.88 | −10.23 | −16.46 | +28.62 | +19.66 | −4.59 | −13.55 |
  | 187L p-xylene | −4.67 | +11.66 | +4.96 | +2.45 | +4.25 | +6.72 | −18.56 | −23.50 | +16.59 | +9.21 | −15.41 | −22.78 |
  | 188L o-xylene | −4.60 | +11.26 | +5.04 | +2.22 | +4.00 | +6.18 | −19.28 | −24.36 | +15.41 | +8.30 | −16.67 | −23.79 |
  | 1NHB ethylbenzene | −5.76 | +13.92 | +5.11 | +3.12 | +5.69 | +8.84 | −14.70 | −19.78 | +19.63 | +12.34 | −10.68 | −17.98 |

  **Where it has converged: a 20 Å pocket** (119–122 residues, 1945–2005 atoms), at the crystal
  pose and minimised from it, in a new ignored series test. Not the whole protein, because **QEq on
  184L's whole protein does not settle**: its hydrogen iteration still changes a charge by
  1.8e-9 e after the default 100 solves, against a tolerance of 1e-10. That is recorded, not fixed
  here. The test names an entry whose QEq fails and leaves it out of the statistics rather than
  stopping; it fails on any other error, on a skipped set other than the expected one (empty), and
  on fewer than three entries measured. At 20 Å all nine converged. The whole protein, measured for the two entries
  `the_polar_term_diagnosed.rs` builds it for, is within 0.10 (benzene) and 0.25 (n-butylbenzene)
  kcal/mol of 20 Å's polar term (in brackets). Kcal/mol:

  | entry | ΔG°exp | crystal: polar (whole protein) | pocket by volume | ligand's own | cross | empty cavity | minimised: polar | empty cavity | total | total, empty |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 181L benzene | −5.19 | +14.86 (+14.96) | +5.89 | +2.36 | +6.61 | +8.97 | +15.19 | +9.11 | −10.91 | −16.99 |
  | 182L 2,3-benzofuran | −5.46 | +17.67 | +8.30 | +4.78 | +4.59 | +9.37 | +19.11 | +10.08 | −12.89 | −21.92 |
  | 183L indene | −5.13 | +19.56 | +9.41 | +2.73 | +7.42 | +10.15 | +19.95 | +10.92 | −8.76 | −17.78 |
  | 184L isobutylbenzene | −6.51 | +24.28 | +10.89 | +5.21 | +8.18 | +13.39 | +36.43 | +24.18 | −2.16 | −14.42 |
  | 185L indole | −4.89 | +19.24 | +8.80 | +3.68 | +6.76 | +10.44 | +20.54 | +11.45 | −10.54 | −19.63 |
  | 186L n-butylbenzene | −6.70 | +30.76 (+31.01) | +11.08 | +4.19 | +15.50 | +19.68 | +34.81 | +22.55 | −0.29 | −12.55 |
  | 187L p-xylene | −4.67 | +19.31 | +9.41 | +2.70 | +7.19 | +9.90 | +20.63 | +10.92 | −13.40 | −23.11 |
  | 188L o-xylene | −4.60 | +18.49 | +9.39 | +2.53 | +6.57 | +9.10 | +19.10 | +9.98 | −14.32 | −23.44 |
  | 1NHB ethylbenzene | −5.76 | +22.63 | +9.18 | +3.43 | +10.02 | +13.45 | +23.95 | +14.55 | −8.02 | −17.42 |

  **Against ΔG°exp, n = 9, 95% Fisher intervals** (nothing asserted):

  | column | 6 Å Pearson r | Spearman ρ | 8 Å Pearson r | Spearman ρ | 20 Å Pearson r | Spearman ρ |
  | --- | --- | --- | --- | --- | --- | --- |
  | ΔE vac, minimised from crystal | +0.46 [−0.29, +0.86] | +0.32 [−0.46, +0.82] | +0.49 [−0.25, +0.87] | +0.33 [−0.44, +0.82] | +0.50 [−0.25, +0.87] | +0.32 [−0.46, +0.82] |
  | total, minimised from crystal | −0.77 [−0.95, −0.22] | −0.72 [−0.94, −0.08] | −0.90 [−0.98, −0.58] | −0.88 [−0.98, −0.51] | −0.93 [−0.99, −0.69] | −0.85 [−0.97, −0.41] |
  | the same, cavity empty | −0.50 [−0.87, +0.25] | −0.48 [−0.87, +0.29] | −0.75 [−0.94, −0.18] | −0.77 [−0.95, −0.19] | −0.87 [−0.97, −0.47] | −0.87 [−0.97, −0.46] |
  | total, relaxed and minimised | −0.81 [−0.96, −0.33] | −0.80 [−0.96, −0.27] | −0.87 [−0.97, −0.50] | −0.70 [−0.93, −0.04] | | |
  | the same, cavity empty | −0.68 [−0.92, −0.02] | −0.70 [−0.93, −0.04] | −0.82 [−0.96, −0.34] | −0.83 [−0.97, −0.36] | | |
  | polar alone, relaxed and minimised (20 Å: crystal) | −0.88 [−0.98, −0.53] | −0.80 [−0.96, −0.27] | −0.88 [−0.97, −0.50] | −0.52 [−0.88, +0.25] | −0.81 [−0.96, −0.30] | −0.58 [−0.90, +0.16] |
  | the same, cavity empty | −0.93 [−0.98, −0.69] | −0.82 [−0.96, −0.31] | −0.91 [−0.98, −0.63] | −0.62 [−0.91, +0.10] | −0.85 [−0.97, −0.41] | −0.65 [−0.92, +0.05] |

  **The honest conclusion: generalized Born with QEq charges in a rigid pocket cannot rank this
  series, and neither converging it nor emptying the cavity changes that.** The polar term alone
  anticorrelates with experiment (r −0.81 to −0.93), because it is largest for the ligands that
  bind best, the two largest (+34.8 and +36.4 minimised at 20 Å, against +15.2 for benzene), and
  its spread across the series (9.5–17.7 kcal/mol) is four to eight times experiment's 2.10. Converging it
  makes the solvated total rank the series more firmly backwards (−0.77, −0.90, −0.93 at 6, 8 and
  20 Å). The empty cavity takes out the pocket's part, the one 2c-2 named, and the total still
  anticorrelates past the 0.666 at which nine points differ from zero, at 8 and 20 Å: what is left
  is the ligand's own desolvation and **the cross terms, the largest part for the two largest
  ligands** (+18.1 and +17.6 minimised at 20 Å). The cross terms are GB's reaction field between
  the pocket's QEq charges and the ligand's, which nothing here has checked against a
  Poisson–Boltzmann or explicit-solvent reference, and with QEq's charges OBC II already
  over-solvates small molecules by a mean −2.1 kcal/mol. The vacuum ΔE's correlation, +0.46 to
  +0.50, has an interval containing zero at every cutoff. So the series needs either a solvation
  model checked on buried pockets or charges OBC II was fitted with, and neither is in this step.
  The series test takes 124 s with `--release -- --ignored`, the 20 Å test 883 s; the benzene
  file's default tests 13.2–13.8 s unoptimised, against 13.4–13.8 before. 210 tests in the crate,
  and nine ignored.
- **`pantometry-forcefield` moves: molecular dynamics, BAOAB Langevin and velocity Verlet, and
  `Molecule` as a domain whose step can be a time step.** This is step 3a of sampling and complex
  stability. New module `dynamics`: `MolecularDynamics` integrates any `Potential`. `ForceField`
  is one, and a test's harmonic well another.
  - **The integrator.** BAOAB is the splitting of Leimkuhler and Matthews, *Appl. Math. Res.
    Express* 2013, 34, read in arXiv:1203.5428, whose Appendix gives the five lines. Its
    Ornstein–Uhlenbeck `O` is solved exactly. `Bath::Isolated` gives velocity Verlet instead. One
    force evaluation per step.
  - **Masses** are CIAAW's abridged standard atomic weights (2024 table) over `N_A`, as
    `Element::atomic_weight` and `Element::mass`.
  - **Randomness.** Every kick is `Rng::for_index(seed ^ KICK_STREAM, step·N + atom)`, with `N`
    counting frozen atoms. So at a fixed step a run is the same bits however it is cut into calls,
    and freezing one atom does not change another's noise. A `Simulation` whose frames are not
    whole multiples of the step takes other substeps, and is another run.
  - **Frozen atoms** have zero velocity and are never written.
  - **Initial velocities.** `thermalised` draws Maxwell–Boltzmann velocities. With nothing frozen
    it removes the centre-of-mass **and** the angular momentum: velocity Verlet conserves both for
    an isolated molecule, so a tumble drawn at the start would never reach the internal modes. The
    count is `3N − 6` in NVE (`3N − 5` if linear), and `3 N_free` under a bath, which
    re-thermalises the removed motion.
  - **The books.** The bath's work is booked per `O`. `ledger()` is `energy` = kinetic +
    (potential − its value at the start) − work, in three entries.
  - **The domain.** `Molecule::thermalised` or `with_dynamics` puts a molecule in motion. Its
    `step` then advances `dt` in steps of at most `Molecule::time_step` (0.5 fs), and
    `max_stable_dt` says so, so the kernel cuts each frame.
  - **New readings:** `temperature` (°C, after `O` under a bath), `kinetic energy`, `thermostat
    work`, `conserved energy` (the ledger's books) and `dynamics steps`. In motion the diagnostics
    are `max force`, `rms force`, `conserved energy` and `dynamics steps`. `domain_readings` now
    also holds a molecule in motion to them.
  - **The minimiser is still the default**, so no scene, reading or pinned label changes.
  - **The enum is `Bath`, not `Thermostat`.** The prelude already exports `pantometry-molecular`'s
    `Thermostat`, which has variants of the same names and different fields.

  **The time step, from measurement.** NVE aspirin from its relaxed geometry at 300 K, for 1 ps at
  each step:

  | dt (fs) | RMS ΔE (kcal/mol) | max \|ΔE\| | drift (kcal/mol/ps) |
  | --- | --- | --- | --- |
  | 0.125 | 0.0040 | 0.0070 | −0.0003 |
  | 0.25 | 0.0162 | 0.0281 | −0.0012 |
  | 0.5 | 0.0654 | 0.1147 | −0.0059 |
  | 1 | 0.2738 | 0.4468 | +0.0091 |
  | 1.5 | 0.6637 | 1.0748 | −0.0143 |
  | 2 | 1.3667 | 2.1629 | −0.0007 |
  | 2.5 | 2.8027 | 4.5741 | −0.0448 |
  | 3 | blew up at 27 fs | | |

  The error quarters per halving up to 1 fs and departs from `h²` above 2 fs. The stability edge is
  between 2.5 and 3 fs, below the 3.7 fs that `ωh = 2` gives a lone C–H (UFF's `k` = 662 kcal
  mol⁻¹ Å⁻², period 11.5 fs).
  - **0.5 fs**: its worst departure is a fifth of `k_BT`, where 1 fs's is three quarters. BAOAB's
    whole-step kinetic bias on a C–H, `(ωh)²/4`, is 1.9% at 0.5 fs and 7.5% at 1 fs.
  - **No constraints.** They would buy about a fourfold step, and nothing here needs that yet.

  **Checked against closed forms, not pinned outputs** (`tests/the_dynamics.rs`,
  `tests/a_molecule_in_motion.rs`):
  - **Velocity Verlet on a harmonic well.** Its largest energy departure is the shadow
    Hamiltonian's `E₀(ωh)²/4`, met to 3e-9 at ωh = 0.4, 0.2 and 0.1. The measured order is 2.000.
  - **One C–H bond**, released along its axis, follows Verlet's exact discrete cosine at
    `√(k/μ)`. Measured: 4.8e-13 of the amplitude over 2000 steps. Bound: 1e-11, the rounding summed
    linearly: 5.2e-12 from the positions and 1.6e-12 from the closed form's `acos`.
  - **BAOAB on a harmonic well**, 150 000 steps at each of ωh = 0.5, 1 and 1.9:
    - `k⟨x²⟩/k_BT` = 1.0012, 1.0016 and 1.0006 against 1 exactly;
    - the whole-step `⟨mv²⟩/k_BT` = 0.9412, 0.7510 and 0.0977 against `1 − (ωh)²/4`;
    - after `O` it is 1.0038, 1.0014 and 1.0008 against 1;
    - every |z| is under 1.3, against a bound of 4 standard errors from the series' own
      autocorrelation time.
  - **Two unequal masses (C and H) under one bath** each sample `k⟨x²⟩ = k_BT`, and their kicks
    are uncorrelated between atoms and between components.
  - **Each BAOAB step moves the books by exactly `(a²/2)(P₁² − P₂² + Q₂² − Q²)`** on a harmonic
    well, with `P = √m v`, `Q = √k x` and `a = ωh/2`, worked by hand from the five lines. Measured
    to 9.7e-16 of the scale over 1000 steps, against a rounding bound of 1e-13. Its leading part is
    the `h²(P₁² − P₂²)/8` that makes a bath's books a random walk.
  - **A 0 K bath at `γh = ln 2`** is the Appendix worked by hand.
  - **Thermalised velocities are χ² on `3N − 6`** over 8000 seeds: mean 1.0025 (σ 0.0021) and
    variance 1.037 of `2(3N − 6)` (σ 0.0166). Both bounds are 4σ, with no autocorrelation, since
    the draws are independent.
  - **Equipartition on aspirin** in a 20 ps⁻¹ bath, 30 000 steps, in total and **per element**:
    - total: 1.0023 ± 0.0094 with nothing frozen, 1.0078 ± 0.0114 with the six ring carbons
      frozen;
    - hydrogens 1.008 ± 0.014, carbons 1.004 ± 0.012 and oxygens 0.986 ± 0.016 with nothing
      frozen;
    - bound: 4σ from each series' autocorrelation time.
  - **Energy-error orders.** NVE aspirin's energy error falls as `h²` and its per-step change as
    `h³`. The bath's books also fall as `h²`, but that bound (0.8) is weak; the closed form above
    is the real check of the bath's books.
  - **Masses.** Every atomic weight equals the CIAAW table typed into the test. The masses sum to
    all 49 entries' `formula_weight` within CIAAW's stated uncertainties. Each atom gets its own
    element's mass.
  - **Determinism and invariance.**
    - A run is the same bits in one call, three, or 300, and from a rebuilt start.
    - Frozen atoms keep their bits.
    - A translated and rotated start gives the moved trajectory to 7.8e-15 Å. That holds in NVE
      only: the kicks are lab-frame vectors.
  - **The domain.**
    - Kernel frames of 0.5, 1 or 5 fs give the same bits, and so does a direct 5 fs `step`.
    - A restored run repeats itself, noise included.
    - New charges replace the cached forces and the books' zero.

  **Found on the way.**
  - **The books' scale held UFF's arbitrary zero, and a half-counted bath passed the audit.**
    `ledger()` first booked the whole potential, and aspirin's minimum is 29.6 kcal/mol. The audit
    judges a change against the largest entry, so the zero was part of the scale every change was
    measured against. With the bath's work half counted, a one-step frame moved the books by
    5.0e-3 of that scale, inside the 1e-2 tolerance, and the audited run passed.
    - The potential is now booked from its value at the start (`potential_reference`, retaken by
      `forget_forces`).
    - `Molecule::ENERGY_TOLERANCE` is re-earned on the new scale at **7e-3**. Over 20 ps at 300 K
      and four seeds, the worst correct frame moved the books by 2.63e-3 (one step), 3.20e-3
      (50 fs) and 4.35e-3 (1 ps). Half-counted work moves them by at least 1.06e-2, 6.7e-2 and
      0.146.
    - The audited test now runs frames of one step and of 50 fs, and fails on half-counted work.
    - The tolerance scales as `(dt / 0.5 fs)²` for a longer step.
  - **A bath makes the books a random walk.** One start's departure measures its order with sd
    0.21 over 32 orders, against NVE's 0.029. Pooling four starts did not narrow it.
  - **Holes the first tests had**, each found by a sabotage that passed them:
    - one noise vector for every atom passed equipartition;
    - the noise width from the mean free mass passed, with hydrogens near 155 K and heavy atoms
      near 386 K and the total exact;
    - velocities drawn 5% cold passed the 400-seed χ² at +1.5σ;
    - iodine at 12.690, chlorine at 53.45 and fluorine at 19.998 passed, because no aspirin-sized
      entry weighs them closely;
    - the kinetic energy read before `O` passed, because it is the same in distribution for a
      harmonic well;
    - an unbounded `max_stable_dt` passed, because `step` also cuts a long `dt` itself.

    Each now fails a test written for it. **One sabotage is equivalent**: drifting a frozen atom
    moves it by zero.

  **Sabotage.** Thirty-eight sabotages were run, each restored by copying the original back,
  `touch` and SHA-256. The new ones were run before their fix and after it:

  | sabotage | before | after |
  | --- | --- | --- |
  | noise width from the mean free mass | passed | caught by per-element equipartition and the two-mass well |
  | velocities 5% cold | passed | caught by the 8000-seed χ² |
  | I 12.690, Cl 53.45, F 19.998 | passed | caught by the CIAAW table |
  | half-counted work, against `a_molecule_in_motion` | passed (with the whole potential booked) | caught by the audited run |
  | the kinetic energy after `O` booked | caught (by the aspirin books tests) | also caught by the closed-form books |
  | drift by `dt(1 + 1e-6)` | caught (by the shadow test) | also caught by the C–H bond at 1e-11 |
  | masses in reverse order | caught (by the books tests) | also caught by the per-atom mass assertion |

  The earlier sabotages are all still caught after the changes. They are:
  - B/A/O reordered (BOAAB);
  - the friction exponent's sign, and the noise factor's;
  - kicks reused across steps, and one kick for every atom;
  - the second half-kick dropped, and a full first kick;
  - H's mass set to C's, and a mass in grams;
  - the bath's work not counted, or half counted;
  - the kinetic energy read before `O`;
  - frozen atoms kicked, and the noise index counting free atoms only;
  - the rigid motion not removed, or only its linear part;
  - the bath's degrees of freedom less the removed six;
  - a `2k_BT` width, and an anisotropic kick;
  - stale forces after new charges, and a checkpoint without the dynamics;
  - an empty ledger in motion;
  - an unbounded step, no internal substeps, and the step count advanced per call;
  - the temperature taken at the whole step;
  - a declared diagnostic renamed (caught by `domain_readings`).

  **Counts.** Unoptimised, the new default tests take 2.1 s in `the_dynamics.rs` and 0.1 s in
  `a_molecule_in_motion.rs`. `cargo test -p pantometry-forcefield -- --list` counts **245 tests,
  twelve of them ignored**: 233 run by default. 26 are new, three of them ignored measurements
  that print tables and assert nothing. The previous entry's "210 tests and nine ignored" counted
  the tests that run, and 210 + 23 is the 233.

- **`pantometry-forcefield` puts a complex in motion: a mobile zone around the ligand, a frozen
  buffer around that, and benzene measured staying bound in T4 lysozyme L99A.** This is step 3b.
  New module `complex`.
  - **`Complex::new(binding, mobile_cutoff, solvent)`** frees the ligand and every residue with a
    heavy atom within the cutoff of a ligand atom at the crystal pose. That is the pocket's own
    rule, so a 6 Å zone in a 10 Å binding is exactly 2c-1's 6 Å pocket: 18 residues and 322 atoms.
    The rest of the binding is a frozen buffer that keeps its bits. `Complex::with_mobile` takes
    any mask.
  - **A mobile atom bonded across the cut is refused**, as `ComplexError::MobileAtTheCut`: its
    backbone would end free in vacuum. Whole residues reach far, so this decides the buffer:
    - the 6 Å zone needs a 10 Å binding for benzene and 12 Å for n-butylbenzene, whose Phe153
      is cut at 10 Å;
    - the 8 Å zone needs 14 Å, because Trp126 is cut at 10 and 12.

    `Binding` gained `residue_spans`, `bonded_outside`, `crystal_positions` and `cross_terms_at`
    to carry this.
  - **Charges are the binding's**: QEq on the pocket and on the ligand separately, fixed, for
    mobile and frozen atoms alike. They are the same model as the binding energies.
  - **The potential** is the complex's force field less the terms among frozen atoms. A test
    holds its force on every mobile atom to the whole complex's, bit for bit.
    `Solvent::GeneralizedBorn` adds OBC II over every atom, since the Born radii are many-body.
  - **`Complex::run`** writes a `Frame` into a `Record` every so many steps, keyed by the
    dynamics' own step count, so a run cut into calls records the same frames. Each frame holds:
    - the ligand's RMSD without superposition;
    - its **site RMSD**, the same blind to which atom sits on which crystal site;
    - the centroid displacement, and the distance to the cavity centre (the lining atoms where
      they are now);
    - the heavy-atom contacts below 4 Å;
    - the vacuum interaction;
    - the mobile atoms' temperature after `O`;
    - the books.

    The `Record` also keeps each atom's mean-square fluctuation. `Estimate::of` gives a mean and
    its standard error from Sokal's autocorrelation time. `complex::mean_square_displacement` is
    `3B/(8π²)`.

  **The cost of a step** (release, one core):

  | zone in binding | atoms (mobile) | ms/step | ns/day |
  | --- | --- | --- | --- |
  | 6 Å in 10 Å, vacuum | 987 (322) | 1.63 | 26.5 |
  | 6 Å in 14 Å, vacuum | 1486 (322) | 2.56 | 16.9 |
  | 8 Å in 14 Å, vacuum | 1486 (719) | 5.00 | 8.6 |
  | 6 Å in 12 Å, vacuum, n-butylbenzene | 1388 (439) | 3.09 | 14.0 |
  | 6 Å in 10 Å, OBC II | 987 (322) | 87.2 | 0.50 |

  **Generalized Born is 53 times the vacuum step.** 100 ps of it would take nearly five hours, so
  it was run for 8 ps, which is one correlation time of the slowest vacuum observable.

  **What the runs measured.** Each run starts from the hydrogen-relaxed binding (A-1) and
  minimises the mobile zone. It then runs at 300 K, 0.5 fs, in a 1 ps⁻¹ bath, discards 10 ps
  (2 ps under GB) and takes a frame every 10 fs. Means are given ± the standard error from the
  autocorrelation time:

  | | 6 in 10 Å, 100 ps | 6 in 14 Å, 50 ps | 8 in 14 Å, 50 ps | OBC II, 8 ps |
  | --- | --- | --- | --- | --- |
  | RMSD (Å) | 2.03 ± 0.33 | 2.14 ± 0.32 | 2.065 ± 0.007 | 0.79 ± 0.03 |
  | site RMSD (Å) | 0.579 ± 0.014 | 0.593 ± 0.015 | 0.656 ± 0.003 | 0.604 ± 0.014 |
  | centroid displacement (Å) | 0.500 ± 0.016 | 0.528 ± 0.013 | 0.557 ± 0.007 | 0.518 ± 0.026 |
  | to the cavity centre (Å) | 1.197 ± 0.019 | 1.234 ± 0.024 | 1.200 ± 0.009 | 1.289 ± 0.027 |
  | contacts < 4 Å | 22.76 ± 0.10 | 22.87 ± 0.11 | 22.22 ± 0.12 | 21.78 ± 0.20 |
  | interaction (kcal/mol) | −23.42 ± 0.03 | −23.58 ± 0.04 | −24.21 ± 0.05 | −22.85 ± 0.10 |
  | mobile temperature (K) | 299.96 ± 1.37 | 301.7 ± 1.7 | 300.9 ± 0.9 | 299.6 ± 4.2 |

  - **Benzene never leaves.** Its centroid's largest displacement in any run is 1.16 Å, against
    a stated threshold of 3 Å.
  - **The 2 Å RMSD is the ring turning in its plane, and it does not sit where the crystal puts
    it.** In the 6 Å runs the RMSD swings between 0.14 and 2.91 Å with τ = 8 ps, while the centroid
    stays near 0.5 Å. Each atom's fluctuation about its mean is 1.94 Å², about `a²` (1.93), the
    value for a ring that turns all the way round.
    - **The in-plane turn**, a new observable (`Frame::turn`, `complex::in_plane_turn`), measures
      it directly. Folded into the 60° sector of the crystal orientation it sits at, the turn
      spreads with an RMS of 22.6°, 23.1° and 25.4° in the three vacuum runs, against 17.3° for a
      ring that fills its sectors evenly. **It peaks 20–30° from the crystal orientation**, half a
      sector away, where the ring's atoms sit between the crystal's sites. Only 13%, 10% and 1.2% of
      frames are within ±10° of the crystal orientation, against a third if the turn were uniform.
    - **It turns as well as offsetting.** In the 6 Å zone in 10 Å it turned −300° net over
      100 ps, five sixths of a revolution, and in 6 in 14 Å −101° net. In the 8 Å zone it turned
      −30° net, staying in two neighbouring sectors.
    - 100 ps holds about six of the RMSD's correlation times, so its ±0.33 is itself uncertain by
      about a factor of two. The other observables decorrelate within a picosecond.
  - **The zone and the buffer move the means little.** A 14 Å buffer moves every mean by about its
    error. The 8 Å zone binds 0.8 kcal/mol more strongly.
  - **Under GB the interaction is 0.6 kcal/mol weaker, and the ring did not turn in its 8 ps.** GB
    is reported as a check of how far a solvent moves the pose, not as the better answer: with
    QEq charges OBC II over-solvates small molecules, and A-2 found its polar term does not rank
    the series.
  - **The books** moved between −0.85 and +1.06 kcal/mol over the runs: the random walk a bath
    makes of them (3a).
  - **n-butylbenzene (186L), 50 ps in vacuum, also stays bound.**
    - Centroid displacement 0.764 ± 0.007 Å, at most 1.32.
    - Interaction −40.74 ± 0.06 kcal/mol, of which van der Waals is −38.79.
    - Minimising moved it 0.99 Å from the crystal, and it stayed there: RMSD 1.190 ± 0.005 Å.

  **Against the B-factors, reported, not asserted.**
  - **Benzene: the simulated ring disagrees with the crystal's B-factors, and that is a finding
    about the model.** Its carbons have B = 20.05–26.80 Å² in 181L, a mean `3B/(8π²)` of
    0.894 Å², an RMS displacement of 0.95 Å. That is a ring held at one orientation with its atoms
    smeared by about an ångström. The simulation's ring sits half a sector from that orientation
    and turns, so its per-atom fluctuation is 1.94 Å² in the 6 Å zone in 10 Å, and 1.09 for the same zone in a 14 Å binding.
    Neither run's ring is the density the crystal measured. Vacuum, the frozen buffer and UFF's
    cavity are each a candidate cause, and none is separated here.
    - **No folded fluctuation is reported.** Folding 60° jumps back onto the crystal's sites is
      only meaningful if the ring jumps between them, and it does not. It sits between them.
    - **The site RMSD is not a fluctuation to set beside B.** It is a nearest-site distance, and it
      has a ceiling. A ring of the crystal's radius (1.367 Å) turning uniformly reads
      `2a²(1 − 3/π)` = 0.168 Å², however freely it turns. An isotropic Gaussian with the crystal's
      own `⟨|u|²⟩` of 0.894 Å² reads back only 0.828 of itself, by Monte Carlo over the crystal
      ring. The runs' 0.350, 0.365 and 0.433 Å² are above the free-turning ceiling because the
      ring is offset: at 30° from its sites each atom is `2a sin 15°` = 0.71 Å from one, which is
      0.50 Å². **The earlier reading, "0.39–0.48 of the crystal's, so no looser than the crystal",
      is withdrawn.**
  - **The mobile protein's heavy atoms fluctuate by a tenth of their B-factors**: 0.08–0.12 Å²
    against 0.97.
  - **The correlation with B** (6 Å in 10 Å):
    - r = +0.67 atom by atom, 95% interval [+0.57, +0.75];
    - r = +0.79 by residue, [+0.52, +0.92];
    - led by the surface residues Lys83, Lys85 and Asp89 in both.
  - **In n-butylbenzene's run** the correlation is weaker: +0.35 by atom, +0.29 by residue. Two of
    its butyl carbons carry B = 100.00 Å², which looks like a refinement's ceiling.
  - **What the protein's comparison can mean**: the simulation ranks the same atoms as mobile as
    the crystal does.
  - **What it cannot mean**: a B-factor includes static disorder and the lattice's motion. Here
    the buffer is frozen, there is no lattice and no water, and the run is 100 ps. So a smaller
    simulated fluctuation is expected for the protein, and its size is not a test of the force
    field. The ligand is the opposite case: its fluctuation is larger than the crystal allows.

  **Checked as exact facts and closed forms** (`tests/the_complex_in_motion.rs`; thirteen default
  tests, 2.8 s unoptimised, on a 3 Å complex of benzene and four side chains with every backbone
  frozen):
  - **Bits and masses.**
    - The buffer keeps its bits and has zero velocity.
    - Each atom has its own element's mass.
    - The degrees of freedom are `3 N_mobile`, with and without a bath.
  - **The potential's force** on each mobile atom is the whole complex's, bit for bit, in vacuum
    and under OBC II.
  - **The frozen-only terms are gone.** No term of the potential has all its atoms frozen. The
    whole complex's energy less the potential's is the 528 left-out terms' sum, computed term by
    term in the test. The allowance is traced: `n ε Σ|t|` for each of the three sums, plus one ε of
    the subtraction. The measured gap is 6e-34 J against an allowance of 1.5e-28 J.
  - **Determinism.** A run is the same bits in one call, in four, in 120, and rebuilt from the
    file, with a frame interval (7) that divides none of the chunks.
  - **NVE order.** The mobile zone's energy error falls as `h²`, within ±0.15 of 2. That bound
    is earned over sixteen starts: mean 2.041, sd 0.025, range [2.015, 2.070].
  - **Every frame is recomputed from its positions**, with the cross terms summed pair by pair in
    the test.
  - **Observables on hand-built inputs.**
    - A translation's RMSD is its length.
    - A ring turned by θ has RMSD `2a sin(θ/2)`, and site RMSD zero at 60°.
    - The site RMSD of a uniformly turning ring is `2a²(1 − 3/π)`, by a 6000-point midpoint
      average whose error is bounded by `h² · 2a²/24`. Two atoms on one site read zero: nearest,
      not one-to-one.
    - The in-plane turn of a turned, translated ring is the turn.
    - Contacts count strictly below the cutoff, and a 2 × 3 grid gives the 28 counted by hand.
    - At the crystal pose the contacts and the cavity distance are the file's, read with string
      operations.
  - **B-factor.** `3B/(8π²)` matches a typed B to 4 ε.
  - **Autocorrelation.**
    - On one AR(1) series, τ and the standard error are within four of their own spreads of
      `(1 + φ)/(2(1 − φ))`.
    - Over 32 seeds, the mean of `τ̂/τ` is 0.9989 ± 0.0051, within four of its standard errors
      (2.1%). That is what sees Sokal's window constant: stopping at `3τ` biases τ by −4.6%,
      which passes one series. The bias at `6τ`, 0.27%, is inside the bound.
  - **Refusals.** A zone with no buffer, or one across the cut, is refused by name.
    - `Binding::bonded_outside` is checked atom by atom against the test's own neighbour scan.
    - Every atom the scan marks, each residue's backbone `N` **and** `C`, is freed alone, and each
      is refused, naming its residue.
  - **A cutoff's zone is the pocket's rule, at debug speed** (`Complex::zone`, factored out of
    `new`).
    - In the 3 Å binding it is the pocket `Binding::new` cuts at 2.7 Å (Val87 alone) and at 2 Å
      (none).
    - The radii are chosen so that each part of the rule matters: a radius half as large again
      reaches all four residues, and counting the residues' hydrogens would take in Val111.
  - **In release and ignored:**
    - a cutoff's zone is the pocket `Binding::new` cuts at that radius (18 residues and 322 atoms
      at 6 Å, 4 and 82 at 3 Å);
    - equipartition holds per element over the 82 mobile atoms of the 3 Å zone in an 8 Å binding:
      H 0.998 ± 0.006, C 1.008 ± 0.006, N 1.003 ± 0.014, O 0.988 ± 0.017, and the total
      1.001 ± 0.005, each |z| below 1.3 against 4.
    - **No gate runs this.** The per-mass noise width it would catch is held in the default run by
      3a's `aspirin_shares_its_kinetic_energy_equally`, not by this test.

  **Found on the way.**
  - **The first cross-term check verified the code against itself.** It compared each frame's
    cross terms with `Binding::cross_terms_at`, the function that produced them. Sabotaged to read
    the binding's stored positions, it passed. The test now sums every pair itself.
  - **A benzene's RMSD from its crystal pose is the wrong number to compare with a B-factor**, as
    above. The site RMSD was added for this, and then found unearned for it as well, by review. A
    nearest-site distance has a ceiling of 0.168 Å² for a ring turning freely, and reads back 0.83
    of a Gaussian with the crystal's own width. The first version of this entry read "0.39–0.48 of
    the crystal's" as "no looser than the crystal". That reading is withdrawn. The in-plane turn
    was added to measure the orientation directly, and what it measured is the finding above: the
    ring sits half a sector from the crystal orientation.
  - **Five checks did not notice what they claimed to hold**, found by a review's sabotages:
    - the estimator test could not see Sokal's window constant;
    - the cut refusal was checked on a backbone `N` only, so checking only a residue's first atom
      passed;
    - the zone's rule was held only by an ignored test, so a radius half as large again, or no
      heavy-atom filter, passed the defaults;
    - nothing asserted that the frozen-only terms were dropped, so a potential that kept them
      passed;
    - the energy-gap allowance, `1e-12 (|E_all| + |E_zone|)`, was not traced to anything.

    Each now has a check that fails on its sabotage, listed below.
  - **The NVE test cannot see a wrong potential that conserves energy**: a zone that ignored the
    buffer's pull would conserve its own energy just as well. The bit-for-bit force test is the
    check for that.

  **Sabotage.** Twenty sabotages were run first, each restored by copying the original back,
  `touch` and SHA-256. All twenty are caught:
  - the dynamics without the frozen mask;
  - frames keyed by the call's own count;
  - GB dropped from the potential, and a potential of the ligand's terms only;
  - the RMSD taken after removing the centroids, and the site RMSD blind to the element;
  - contacts counted at or below the cutoff;
  - `B/(8π²)` with the 3 dropped;
  - the standard error without its factor of 2, and an autocorrelation window of one lag;
  - the fluctuation's mean squared over `n` instead of `n²`;
  - the cavity centre taken from the crystal;
  - the temperature taken at the whole step;
  - no cut check, and nothing bonded outside;
  - masses in reverse order;
  - the cross terms taken at the binding's positions, caught only after the fix above;
  - the second half kick dropped;
  - the zone rule from the ligand's heavy atoms only (release);
  - the noise width from the mean mass (release).

  Then the review's eight, run against the review's snapshot of the tests (before) and against the
  fixed ones (after), each restored the same way:

  | sabotage | before | after |
  | --- | --- | --- |
  | Sokal's window at `3τ` instead of `6τ` | passed | caught by the 32-seed mean |
  | the cut checked on a residue's first atom only | passed | caught: freeing a `C` is refused |
  | `bonded_outside` marking only bonds to an earlier atom | passed | caught by the atom-by-atom scan |
  | the zone's radius half as large again | passed | caught at debug speed (2.7 Å) |
  | the zone without its heavy-atom filter | passed | caught at debug speed (2 Å and 2.7 Å) |
  | the potential keeping the frozen-only terms | passed | caught: no all-frozen term, and the gap is their sum |
  | the turn without the current centroid | (new) | caught by the turned, translated ring |
  | the turn read from the first atom only | (new) | passed, then caught after a non-rigid closed form was added |

  **Counts.** `cargo test -p pantometry-forcefield -- --list` counts **265 tests, nineteen of them
  ignored**: 246 run by default. The 20 new ones are all in `the_complex_in_motion.rs`: thirteen
  run by default, in 2.8 s unoptimised, and seven are ignored. Of those seven, two assert
  (equipartition and the zone rule), and five are measurements: four trajectories, which assert
  the buffer's bits, and the NVE spread, which asserts nothing.
  245 + 20 is the 265, and 233 + 13 the 246.

- **`pantometry-forcefield` computes free energies: an alchemical Hamiltonian, a Boresch
  restraint, TI and BAR from λ windows, each checked against a closed form, and benzene decoupled
  from T4 lysozyme L99A in vacuum.** This is step 3c-1, the engine that 3c-2 (a faster OBC II) and
  3c-3 (benzene's absolute binding free energy by double decoupling) need. Three new modules:
  `alchemy`, `boresch` and `free_energy`.
  - **`Decoupling`** decouples one group of atoms from the rest of a `ForceField` with three
    couplings in a `Lambda`:
    - the group–environment Coulomb pairs scaled linearly, so `∂U/∂λ_e` is the coupled cross
      Coulomb energy exactly;
    - its van der Waals pairs through a soft core, `λ V_LJ(r_sc)` with
      `r_sc⁶ = α σ⁶ (1 − λ)^p + r⁶`, α = 0.5 and p = 1 as Mobley, Chodera and Dill 2006 used;
    - a `Boresch` restraint switched on linearly.

    The group's own terms are kept at every λ: decoupling, not annihilation, so the ligand's
    intramolecular energy cancels between the two legs of a cycle. Energy, forces and `∂U/∂λ` are
    analytic. A solvated force field, a group bonded across its boundary, and an empty or whole
    group are each refused by name (`AlchemyError`).
  - **`Boresch`** holds the six relative coordinates (one distance, two angles, three dihedrals)
    between three receptor and three ligand atoms. `Boresch::choose` picks the pair of candidate
    triples whose four hinge angles are furthest from straight, with `r₀` in a stated range.
    `release_free_energy` is the analytic free energy of releasing the restraint to 1 M, and
    `release_free_energy_extended` the same with its two approximations undone in closed form.
    `STANDARD_VOLUME` is `10⁻³ m³ / N_A`.
  - **`Windows`** runs BAOAB at each state of a schedule. Each window is deterministic from
    `window_seed(seed, k)`, samples on the dynamics' own step count, so it is the same bits however
    it is cut into calls, and holds its whole state, so a clone is a checkpoint and a stopped
    campaign resumes. A `Sample` is `∂U/∂λ` and the energy differences to both neighbours, from
    the λ-dependent part alone.
  - **TI**: the trapezoid on any path, Simpson's rule on a uniform line. **BAR**: Bennett's
    equation solved by bisection to the last bit, with Shirts et al.'s variance for independent
    samples (`bar`), and the delta method's for a window's correlated samples (`bar_correlated`)
    and for a chain of windows (`bennett_chain`). **Not MBAR**: it would need every sample's
    energy at every state, and BAR between neighbours is what Mobley et al. 2007 used on this
    system.

  **Read, and what was not.**
  - **Beutler et al. 1994**, **Boresch et al. 2003**, **Bennett 1976** and **Shirts et al.
    2003** are all behind paywalls, and none was opened. Each formula is taken from a source that
    was read, and the documentation says which:
    - the soft core from the GROMACS reference manual's free-energy section;
    - the Boresch coordinates and closed form from Clark et al., *J. Chem. Theory Comput.* 19,
      3686 (2023), PMC10308817, Figure 3 and eqs 6–8, read off the published equation images;
    - BAR and its variance from Shirts and Chodera 2008 (arXiv:0801.1426), eq 11 and Appendix E.
  - **Mobley et al. 2006 and 2007** were read through PMC summaries: α, p, the protocol's order
    and the force constant `K₀ = 10` in `K₀(ξ − ξ₀)²`. That is K = 20 in this crate's `½K`
    convention, which is the convention the closed form's `(2π k_BT)³` requires.

  **Found on the way: summing the intervals' BAR variances understates the error.** By 1.7 on
  the particle below, where the samples were thinned as well; by 1.26 on exact, independent
  samples along a five-state chain, which is the covariance alone.
  - The first version thinned each window to every `⌈2τ⌉`-th sample and summed Shirts's variance
    over the intervals, as common practice does.
  - On a Lennard-Jones particle decoupled from a fixed atom, over 48 seeds, the mean was right:
    −0.5000 ± 0.0050 `k_BT` against −0.5062 by quadrature. But σ̂ was 0.0198 against a spread of
    0.0345. At a stride five times longer, where the samples are close to independent, it was
    0.0102 against 0.0166.
  - So thinning was not the cause. Adjacent intervals share a window — its samples are the
    forward set of one interval and the reverse set of the other — and the covariance was missing.
  - The error is now the delta method's, one series per window, so both intervals a window feeds
    are in it, with each series' own autocorrelation time. Every sample is used.
  - TI was not affected: each window enters it once, and its σ̂ was 0.0299 against a spread of
    0.0333.
  - **BAR is solved by Newton's method inside a bisection bracket**: each evaluation moves one
    end of the bracket, and the search stops when no float is left strictly inside it. Every
    printed result is the bisection's to the digit, at a fraction of the evaluations, which is what
    paid for the calibrations below.

  **The benzene demonstration (reported, not asserted).** 181L in vacuum, with the 3b setup: a
  6 Å mobile zone in a 10 Å binding, hydrogens relaxed, the zone minimised and the buffer frozen.
  - **The restraint** chosen by the rule is on the backbone C, CA and N of one pocket residue
    (complex atoms 177, 176 and 175) and benzene's C3, C4 and C5. Its reference: r₀ 6.452 Å,
    θ_A 77.2°, θ_B 105.5°, and hinges between 77° and 120°.
  - **The schedule** has 22 windows: the restraint on in six, the charges off in four more and
    the van der Waals off in twelve. Each window runs at 0.5 fs in a 1 ps⁻¹ bath, discards 2 ps
    and samples 10 ps every 10 fs.

  | segment | BAR (kcal/mol) |
  | --- | --- |
  | restraint on, coupled | +3.053 ± 0.077 |
  | charges off | +1.477 ± 0.009 |
  | van der Waals off | +13.799 ± 0.127 |
  | **complex leg** | **+18.329 ± 0.200** (TI, trapezoid: +18.421 ± 0.202) |
  | restraint released to 1 M, analytic | −7.835 |

  - **Every interval's overlap is between 0.38 and 0.50**, so the schedule is denser than it
    needs to be.
  - **The restraint is not free to turn on.** At λ_r = 0, `⟨U_B⟩` is 12.2 kcal/mol, because the
    unrestrained ring turns away from the minimised pose the restraint holds it to (3b's finding).
  - **The cost** was 873 s for the 22 windows: 39.7 s per window of 24 000 steps, 1.65 ms per step.
  - **A vacuum leg has no solvent counterpart and cannot be compared with experiment.** It is here
    to exercise the engine end to end.

  **Checked against closed forms** (`tests/the_free_energy_against_closed_forms.rs`, fifteen
  default tests in 8.5 s unoptimised as libtest runs them, 18.7 s one at a time;
  `tests/benzene_decoupled_from_its_pocket.rs`, four in 0.4 s). The two harmonic wells below have
  `F₁ − F₀ = (3/2) k_BT ln(k₁/k₀)` exactly. **The calibration bands are set by the number of
  seeds**, so the exact-sample tests run many small campaigns: bands of ±0.14 to ±0.16 on the
  z-scores' variance, which a σ̂ wrong by about 8% leaves.
  - **BAR on exact, independent samples** of the two wells: 1600 seeds, `N_F = 400` and
    `N_R = 700`, so that `M ≠ 0` and its sign is seen. The z variance is 0.963 with Shirts's
    formula and 0.916 with the delta method, against 1 ± 0.14.
    - **The delta method is about 4% conservative on independent samples**: `Estimate::of`'s `τ`
      is a windowed sum of noisy correlations clipped at ½ from below, so it errs upward only.
  - **On correlated samples** (stationary AR(1) chains with φ = 0.8, 800 per state, 1600 seeds):
    the delta method's z variance is 1.001, against 1 ± 0.14; Shirts's formula, which assumes
    independence, gives 3.66. The 0.821 the first version of this test printed, at 200 seeds, was
    a fluctuation of a band four times as wide.
  - **Along a chain of five states** (150 samples per window, 1200 seeds): window by window, the
    z variance is 0.982; summing the intervals' variances gives 1.555, outside the band of ±0.16.
  - **Through `Windows` itself, over 600 small campaigns, on a schedule with a corner**: one atom
    whose stiffness goes `k₀ → 2k₀` under λ_r and then `→ 4k₀` under λ_e, so the corner window's
    trapezoid weight spans two components and both BAR intervals meet there.
    - Unbiased: TI's mean deviation from the trapezoid on the exact integrand is
      −0.0049 ± 0.0043 `k_BT`, and BAR's from the exact ΔF −0.0053 ± 0.0036.
    - Calibrated: the z variance is 1.047 for TI and 1.032 for BAR, against 1 ± 0.23.
    - The z means, −0.11, are not asserted: σ̂ grows with the estimate here, so an unbiased
      estimate's z leans negative; the deviations are what check bias.
    - With 100 steps of equilibration from rest the variance read 1.14 for both, because the first
      samples remembered the start. The test equilibrates for 500.
    - Simpson's rule refuses the schedule.
  - **TI and BAR from MD windows** on the two wells, eleven windows:
    - Simpson's rule gives 2.0779 ± 0.0167 `k_BT` against 2.0794. Its bias, computed on the
      exact integrand, is 3.4e-4, and is asserted below a tenth of σ̂.
    - The trapezoid is held to its own rule's value on the exact integrand: its bias, +0.0105, is
      not small.
    - BAR gives 2.0839 ± 0.0155.
    - Each window's `⟨∂U/∂λ⟩` is within 4σ̂ of `(3/2) k_BT (k₁ − k₀)/k(λ)`.
  - **Boresch, at 10 kcal mol⁻¹ Å⁻² and rad⁻²**:
    - eq 6 is integrated as `e^(−βU)` with `U` the restraint's own `Boresch::energy`, on
      configurations built from the six coordinates (and checked to give them back), with `8π²`
      from the same quadrature. That gives −7.140864 kcal/mol.
    - The excluded tails (r below the range, θ outside [0, π], |Δφ| > π) are integrated
      separately, and are 2.3e-5 of the integral. With them added, eq 6 is the extended closed
      form to 4e-13 `k_BT`, asserted at 1e-9.
    - The old bound, the tails' Gaussian estimate, was 2.5e-3, a hundred times the tails, and
      caught a missing `r²` correction by 2%. That correction is now six orders above the bound.
    - Eq 7 itself is −7.106823.
  - **The Jacobian of eq 6**: Cartesian grids over A's and B's positions, integrating
    `e^(−βU)` with `U` again `Boresch::energy`, agree with the one-dimensional product to 1e-13 at
    two spacings.
  - **The restraint's force** against central differences of its energy, at scale 0.7, and its
    energy `Σ ½ K (ξ − ξ₀)²` of its own coordinates to 4 ε.
  - **The anchor rule** passes over a hinge at 175°.
  - **The soft core**:
    - at λ = 1 it is UFF's pair to the bit, energy and radial derivative;
    - at λ = 0 it is exactly zero;
    - at `r = 0` it takes its closed form `λD(s² − 2s)` with `s = 2/(α(1 − λ))`;
    - `∂/∂λ` and `∂/∂r` match central differences.
  - **One LJ particle tethered inside a fixed atom's wall**, decoupled over 21 windows:
    - the reference is the configurational integral reduced to one dimension by the angular
      integral `2 sinh(a)/a`, and its quadrature reproduces `Z₀` to 1e-12;
    - BAR and Simpson's TI are each within 4σ̂ of the −0.506 `k_BT` it gives;
    - coupled again from the other end with another seed, BAR closes to +0.017 ± 0.038 and TI to
      −0.0005 ± 0.040.
  - **In 181L's 3 Å pocket**:
    - `∂U/∂λ_e` is the same bits at seven states, and equals both the test's own pair sum and
      `Binding::cross_terms_at`'s;
    - the soft core at λ_v = 1 is the binding's cross van der Waals to the bit;
    - fully coupled, the decoupling is the force field to a gap of 2.8e-14 kcal/mol, against a
      traced allowance of 2.0e-10;
    - decoupled, the forces are the rest's to the bit;
    - at λ = (0.6, 0.4, 0.7), with the Boresch restraint, the forces and `∂U/∂λ` match central
      differences.
  - **Determinism**: a campaign is the same bits in one call, in calls of seven steps, and
    stopped half-way, cloned and resumed. A window alone is the same as among the others.
  - **In release and ignored**, measured:
    - **the particle over 48 seeds.** BAR's mean is −0.5037 ± 0.0044 `k_BT` at 1 fs and
      −0.5069 ± 0.0046 at 0.5 fs, against −0.5062. Its σ̂ is 0.0282 against a spread of 0.0304
      between seeds; TI's is 0.0299 against 0.0333.
    - **the wells over 200 MD campaigns.** The z variance is 0.99 for BAR and 1.09 for TI where
      the samples decorrelate quickly. Where they are strongly correlated (γ = 0.2ω, stride 2) it
      is 1.29 and 1.35, both inside the band of ±0.40 but about 12% small in σ̂. That is the
      autocorrelation-time estimator's limit, and it affects both estimators alike.

  `ForceField` gained `split_across` (crate-internal).

  **Sabotage.** Twenty-six sabotages were run, each restored by copying the original back,
  `touch` and SHA-256. Twenty-four were caught the first time:
  - the soft core's α with the wrong sign;
  - `∂U/∂λ_v` without the chain rule through `r_sc`;
  - no exact branch at λ = 1;
  - the Coulomb term quadratic in λ_e;
  - the restraint's `∂U/∂λ_r` scaled by λ_r, and its coupling energy unscaled;
  - the cross pairs left in the rest as well;
  - BAR's `M` with the wrong sign (caught because the test now uses `N_F ≠ N_R`);
  - BAR's forward term with `Δf`'s sign flipped;
  - Shirts's variance without `−(1/N_F + 1/N_R)`;
  - the reverse set's sign, both in the chain and per interval;
  - the delta method without the autocorrelation time;
  - the chain's covariance term with the wrong sign;
  - the trapezoid without its ½, and Simpson's 4 and 2 swapped;
  - Boresch without `sin θ_A0`, without `8π²`, and with `(π k_BT)³`;
  - `V°` per m³ rather than per litre;
  - the extended correction's sign;
  - the angle gradient's sign, and the dihedral forces unscaled;
  - the anchor rule keeping the straightest hinge.

  Two needed a fix to the test, and were run before and after it:

  | sabotage | before | after |
  | --- | --- | --- |
  | one seed for every window | passed: the check compared windows at different states, which differ anyway | caught: two windows at the same state must differ |
  | samples keyed by the call's own count | hung: a campaign cut into calls never completed, and the loop waiting for it had no bound | caught: the loop is bounded, and incompleteness fails |

  **A review found the uncertainties unchecked**, and eight more sabotages passed. Each was run
  against the tests before the fix and after it, restored the same way. All twenty-six of the
  first set were run again after the fixes and are still caught.

  | sabotage | before | after |
  | --- | --- | --- |
  | `bennett_chain`'s error ×1.2 | passed: band ±0.33 | caught by the chain (±0.16) and the corner campaigns |
  | `bar_correlated`'s σ̂ ÷1.1 | passed: band ±0.40 | caught by the correlated samples (±0.14) |
  | Boresch without its ½, in the energy, the forces and the dihedrals alike (a 2K restraint released as K, 1.24 kcal/mol) | passed: every quadrature re-typed `½K` | caught: the quadratures integrate `Boresch::energy`, and the energy is held to `½K Δξ²` |
  | TI's and BAR's totals' errors ×3 | passed: nothing calibrated σ̂ through `Windows` | caught by the corner campaigns |
  | TI's error ×½ | passed | caught by the corner campaigns |
  | the extended correction without its `r²` term | caught by 2% of the bound | caught by six orders |
  | the trapezoid's weight on one component at a corner | passed: no schedule had a corner | caught by the corner campaigns |
  | Simpson without its straight-line check | passed | caught: the corner schedule is refused |

  **Counts.** `cargo test -p pantometry-forcefield -- --list` counts **287 tests, twenty-two of
  them ignored**: 265 run by default. The 22 new ones are seventeen in
  `the_free_energy_against_closed_forms.rs` (fifteen default, two ignored measurements) and five
  in `benzene_decoupled_from_its_pocket.rs` (four default, the demonstration ignored). 265 + 22 is
  the 287, and 246 + 19 the 265.

- **`pantometry-forcefield`: OBC II in a frozen buffer is twice as fast, to the bit, and the rest
  of its cost is the platform's `exp` and `ln`.** This is step 3c-2. Its target, generalized Born
  within a few times the vacuum step, **is not met**: the step is 28 times vacuum, down from 53.
  - **The profile came first** (release, one core, `x86_64-pc-windows-gnu`). The 987-atom complex
    of 3b (a 6 Å zone of 322 atoms in a 10 Å binding) cost 86.5 ms per evaluation:

    | part | ms |
    | --- | --- |
    | the Born radii: 973 000 pair integrals, a logarithm each | 31.2 |
    | the pair terms: 486 000, an exponential each | 23.0 |
    | the chain rule through the radii: the 973 000 integrals again | 30.8 |
    | UFF | 1.6 |

    On this toolchain `exp` costs 33 ns and `ln` 18.5 ns, against 1.1 ns for `sqrt`.
  - **Three changes, none of which changes a bit.** Every sum still takes the same terms in the
    same order, so the energy, every Born radius and every force are the earlier evaluation's.
    - **Each pair integral is computed once.** The chain rule reads the `dI/dr` that the radii
      were computed with, instead of computing the integral again.
    - **The frozen atoms' descreening is kept.** `GeneralizedBorn::with_frozen(mask, at)`
      computes the integral between each pair of frozen atoms once, and `Complex` gives its model
      the buffer. That is 442 000 of the 973 000 integrals.
      - Every evaluation first checks that each frozen atom is where it was, to the bit. If one has
        moved, every term is computed.
      - `reuses_frozen_terms` says which path an evaluation will take.
      - **Not kept, because they are not constant**: the radii, which every mobile atom changes,
        and so the pair terms between frozen atoms, which depend on them. The brief proposed
        keeping the frozen–frozen pair terms; their `f_ij` depends on `R_i R_j`.
    - **No `e^(−κf)` without salt.** It is `e^0`, exactly one.
  - **No cutoff and no threads.** A cutoff is an approximation, which the brief made opt-in and
    measured, after the exact changes. It was not built: the pair terms fall only as `r⁻¹`.
  - **Measured after:**

    | | before | after |
    | --- | --- | --- |
    | one evaluation | 86.5 ms | 43.5–44.3 ms (57.6 without the kept terms) |
    | the radii | 31.2 ms | 18.5 ms: 532 000 integrals, 35 ns each |
    | the pair terms | 23.0 ms | 23.4 ms: 486 000 pairs, 48 ns each |
    | the chain rule | 30.8 ms | 2.3 ms |
    | a step of dynamics | 87.2 ms | 44.8–45.5 ms, 28 times vacuum's 1.6 ms; 0.96 ns/day |

  - **What limits it is the platform's `exp` and `ln`**, one of which is in every pair term and
    every integral. With both replaced by arithmetic stand-ins of no accuracy, for the timing
    alone, the step took 12.2 ms, about 10.5 of it GB's divisions and other arithmetic. The
    other 33 ms of the 45 are `exp` and `ln`. Two ways past it are left for a decision, because
    each changes something this step was told to keep:
    - an `exp` and `ln` of the crate's own would change the bits of every GB result, though not
      its accuracy, and at about 4 ns a call would bring the step to about 16 ms;
    - a cutoff would change the model.
  - **For 3c-3:** 20 windows of 2 + 10 ps at 0.5 fs is 480 000 steps, about **6.0 h** for the
    complex leg (11.6 h before). The brief's 4.8 M steps and 116 h are ten times that schedule.
    For that many steps it would be 60 h. The solvent leg is the ligand alone, and costs
    microseconds a step.
  - **`GeneralizedBorn::accumulate_reference` and `born_radii_reference` are the code that was
    there before, kept** as what the evaluation is held to.

  **Checked, with no tolerance** (`tests/the_generalized_born_is_its_direct_sum.rs`, five default
  tests in 0.8 s unoptimised, and one ignored). The evaluation is held to the reference bit for
  bit: the energy, the energy alone, every radius, and every force added to a non-zero start.
  - **On aspirin**: OBC II, OBC I, HCT, ε = 4 and 0.15 M salt.
  - **On a cluster that reaches every branch of the pair integral**, with two atoms at one point.
    That puts an `r = 0` pair in the sums.
  - **On 181L's 3 Å pocket**, with its QEq charges.
  - **Under six masks**: no kept terms, none frozen, all frozen, every third, the first half, and
    all but the last. The small complex's own mask, with the backbones frozen, is added.
  - **At three kinds of position**: where the terms were kept; with the mobile atoms moved, where
    the kept terms must be used; and with a frozen atom moved by 0.1 Å and by one ulp, where they
    must not be.
  - **A complex keeps its buffer's terms** through minimising and 40 steps.
  - **Its run is the direct evaluation's run, bit for bit**: 60 steps of BAOAB on the potential,
    and on a `Potential` that adds `accumulate_reference` to the vacuum terms in the order
    `ForceField::evaluate` does. That means 3b's OBC II trajectory is the same trajectory.
  - **In release and ignored**, the same on the 987-atom complex, its run over 20 steps, and the
    cost. Two cost assertions: with the kept terms the evaluation is under 0.9 of its cost without
    them (measured 0.75), and without them under 0.85 of the reference's (measured 0.68).
  - **Unchanged and passing**: every GB closed-form and gradient test, and 3b's bit-for-bit force
    test of the potential against the whole complex.

  **Sabotage.** Twelve sabotages were run, each restored by copying the original back, `touch`
  and SHA-256. The tests ran with `--no-fail-fast`: without it, cargo stops at the first failing
  binary. In the first pass that made three sabotages look caught only by 3b's force test, because
  the new file never ran.
  - Caught by the default tests:
    - a mobile–frozen integral left out of a frozen atom's radius;
    - kept terms used without checking that the frozen atoms are where they were;
    - the positions compared to within 10⁻²⁰ m rather than to the bit (caught by the one-ulp move);
    - a kept slope with the wrong sign;
    - the kept terms summed apart and added at the end (a reassociation alone);
    - the kept matrix transposed;
    - the first frozen atom's kept terms computed from another atom's position;
    - the no-salt shortcut taken in salt;
    - the chain rule skipping `r = 0` before reading its slope (only the coincident pair sees it);
    - the complex's model not told about its buffer.
  - Caught only by the release cost assertions, because each gives the same bits slower:
    - the kept terms never used: 58.3 ms with them, against 57.8 without;
    - the chain rule computing every integral again: 88.4 ms against the reference's 85.7.

  **Counts.** `cargo test -p pantometry-forcefield -- --list` counts **293 tests, twenty-three of
  them ignored**: 270 run by default. The six new ones are all in
  `the_generalized_born_is_its_direct_sum.rs`: five default, and one ignored that asserts. 287 + 6
  is the 293, and 265 + 5 the 270.

- **`pantometry-forcefield` decouples a ligand under OBC II, and benzene's absolute binding free
  energy to T4 lysozyme L99A comes out at +6.1 ± 0.3 kcal/mol by double decoupling, against
  experiment's −5.19 ± 0.16.** This is step 3c-3. Nothing is asserted against experiment.
  - **`Decoupling` accepts a solvated force field.** 3c-1 refused one (`AlchemyError::Solvated`,
    now removed). Its solvation term is taken out of the rest and evaluated as
    `G(x; λ_e, λ_v)` by the new `GeneralizedBorn::decoupled`:
    - **the group's charges are `λ_e q` in every term of eq 2.** Its Born terms and its pairs
      among themselves go as `λ_e²`, its pairs across as `λ_e`. The intramolecular Coulomb term
      stays in the rest at full charge, so the decoupled ligand keeps its vacuum self-interaction
      and loses only its solvation;
    - **the descreening integral between atoms on opposite sides is scaled by `λ_v`**, and between
      atoms on the same side it is not. The ligand's volume leaves the environment's Born radii
      with its van der Waals;
    - `∂G/∂λ_e` and `∂G/∂λ_v` are analytic, and the chain rule through the radii takes a cross
      pair's `dI/dr` times `λ_v`;
    - **no nonpolar term in either leg**, because the surface has no useful gradient (2b). It is
      estimated at the end states and reported separately.

    **The decoupled end state is the environment alone in solvent plus the ligand alone in
    vacuum, exactly.** That is the end state openmmtools' `AbsoluteAlchemicalFactory` builds
    (`_alchemically_modify_GBSAOBCForce`, read at commit `f6ef22a8`): it multiplies an alchemical
    atom's charge, and its descreening of every other atom, by `lambda_electrostatics`. Its path
    differs in two ways, chosen differently here, neither of which changes ΔG:
    - it scales the Born self term linearly, not as the square of a scaled charge;
    - it removes the descreening with the charges, not with the van der Waals.

    Tying the descreening to `λ_v` keeps the radii fixed while the charges go, so that segment is
    exactly quadratic in `λ_e`, and a Born ion's TI is then exact on a trapezoid.
  - A whole molecule may now be the group of a solvated force field: that is the solvent leg.
    `Decoupling::with_solvation` gives a `from_pairs` model a GB term, for a Born ion.
    `Alchemical::couplings` evaluates several states at one configuration, each to the bit of its
    own `coupling`; `Decoupling`'s shares the descreening, and `Window::advance` uses it.
    `GeneralizedBorn::decoupled_energies` is the same for the solvation term alone.

  **Checked** (`tests/benzene_bound_in_generalized_born.rs`, ten default tests in 1.0 s
  unoptimised), on 181L's 3 Å pocket with QEq charges and on one and two ions:
  - **fully coupled, the solvation term is `accumulate`'s to the bit**, energy and every force,
    and the whole energy is the solvated force field's to a gap of 2.8e-14 kcal/mol against a
    traced allowance of 5.2e-10;
  - **with the charges off, the group's charges are in no term**: the solvation term is
    `accumulate` with benzene's charges zeroed, to the bit, and the coupling energy is the cross
    van der Waals plus that, with no cross Coulomb; `∂U/∂λ_e` there is the cross Coulomb plus the
    linear GB terms, rebuilt in the test from the radii and eq 3;
  - **fully decoupled, it is the pocket alone in solvent plus benzene alone in vacuum**: the
    solvation term and its forces are the pocket's own GB to the bit (zero on benzene), and the
    total is the pocket's solvated force field plus benzene's vacuum one to 2.8e-14 kcal/mol
    against 3.8e-10;
  - **the forces and `∂U/∂λ_e`, `∂U/∂λ_v`** against central differences at three states, with
    second-order one-sided differences at an end of [0, 1];
  - **the frozen buffer's kept descreening is used and changes no bit**: with only benzene mobile,
    the decoupled energy, both derivatives and every force at four states equal the same model's
    without the kept terms, at the start, with benzene moved (kept terms used) and with a pocket
    atom moved (not used); and `couplings` over six states, two of them equal, is each state's
    `coupling` to the bit;
  - **a Born ion decoupled through `Windows`** gives `−G = (q²/2ρ̃)(1 − 1/ε) × 332.0637` by the
    trapezoid, by Simpson and by BAR, to 8 ε, the rounding of five weighted terms or four roots of
    `G`'s size (measured: at most 1.69 ε): 101.836305512 kcal/mol for ρ 1.7 Å and q −1, and
    36.927128801 for ρ 1.2 Å and q +0.5. The integrand `2λ_e G` is linear, so the trapezoid is
    exact. **The totals see only the ends of the path**: a Born self term linear in `λ_e`, with a
    consistent derivative, integrates to the same `−G`. What tells the path is the per-sample check
    that each gradient is `2λ_e G`, to 4 ε;
  - **two ions follow eqs 2–8 written out at nine `(λ_e, λ_v)`**: each radius OBC II's with
    `I = λ_v ×` the other's pair integral, energy and `∂/∂λ_e` to 1e-13. `∂/∂λ_v` is held to a
    difference of the written-out energy at `h` = 1e-4, central inside [0, 1] and second-order
    one-sided at its ends, within `h² M₃/6` (or `/3`) plus `4 ε max|G|/h`, with `M₃` the largest
    third derivative on a 0.01 grid, doubled. The measured errors are a quarter to a ninth of that.
    **With both ions frozen**, so that the one kept descreening term crosses the partition, every
    result is the uncached model's to the bit;
  - **the measurements' interval BAR is `Windows::bennett`'s to the bit**: the Born ion's windows,
    read into the records the measurements analyse, give each interval's estimate, variance and
    overlap with the same bits. The helper that computes it, `interval_bar`, is now the one place
    both the insertion rule and the printed table take an interval from;
  - **the cycle**, now the function `binding_free_energy` with the cycle in its documentation, on
    hand-chosen numbers that a flipped release or swapped legs would change;
  - **the solvent leg's molecule**: benzene alone, the group every atom, is the solvated force
    field at full coupling and its vacuum force field to the bit at `λ_e = 0`. Its `∂U/∂λ_v` is
    exactly zero, which **follows from the construction** — a whole-molecule group has no cross
    pair and no cross descreening — and checks only that the code adds nothing there.

  **Sabotage.** Fifteen, each restored by copying the original back, `touch` and SHA-256, run with
  `--no-fail-fast` over this file, 3c-1's two and 3c-2's direct-sum file. All fifteen were caught:
  - by the two-ion closed form: the cross descreening left unscaled, scaled by `λ_v²`, `∂G/∂λ_v`
    left out, the Born self term linear in `λ_e` (openmmtools' form), a pair skipped when its
    scaled charge is zero, and `∂G/∂λ_e` of the Born term without its factor 2;
  - by the finite differences: most of those, and a cross pair's `dI/dr` unscaled in the force,
    a pair's `∂/∂λ_e` without the other charge's scale, and the solvation force not added at the
    dynamics state;
  - by the bit-exact end states: the cross descreening unscaled, the partition's two sums
    swapped, the rest keeping its solvation term (also by the coupled and solvent-leg tests), and
    the radii always taken from the two partial sums even at `λ_v = 1` — a reassociation alone,
    caught by the coupled test to the bit;
  - by the frozen-cache test: a kept frozen term left out of the two partial sums, and states
    deduplicated by `λ_e` alone;
  - by the Born ion and three of 3c-1's `Windows` tests: a window's neighbours evaluated in the
    wrong order.

  **Added after review** (`numerics-reviewer`, which also derived the cycle's signs independently
  and agreed). Four more sabotages, each run before its fix and after it, restored the same way:

  | sabotage | before | after |
  | --- | --- | --- |
  | the reverse set not negated, in the insertion rule and in the analysis | passed: only the ignored measurements reached it | caught by `the_interval_bar_is_the_windows_own` |
  | the release's sign flipped in the cycle | passed: the cycle was only printed | caught by `the_cycle_has_its_signs` |
  | every kept frozen term sent to the same-side sum | passed: the tests froze only the pocket, whose kept terms never cross | caught by the two ions frozen together |
  | the Born self term linear in `λ_e`, its derivative consistent | — | caught by the Born ion's per-sample check, the two ions and the charges-off test |

  The two-ion `∂/∂λ_v` allowances were 1e-4 and 1e-8 of the derivative, unsourced, against
  measured errors of 2.8e-7 and 9e-11; they are now the difference's own error, above. The Born
  ion's 64 ε is now 8 ε.

  **The measurement** (ignored, release, one core). Each window was written to a file outside the
  repository as it ran, so that a stopped run could resume; it ran once, to the end, and exited 0.
  - **The cycle**, in 3c-1's convention: `ΔG°_bind = ΔG_solvent − ΔG_complex − ΔG°_release`. Each
    `ΔG` is the free energy of decoupling in its leg, and `ΔG°_release` is the analytic release of
    the restraint to 1 M, which is negative. Equivalently `+ ΔG°_restrain`, with
    `ΔG°_restrain = −ΔG°_release` the free energy of restraining the decoupled ligand from 1 M.
  - **The setup** is 3b's and 3c-1's: a 6 Å mobile zone of 322 atoms in a 10 Å binding of 987,
    hydrogens relaxed, the zone minimised under OBC II (converged in 1274 steps to 1.8e-3
    kcal/mol/Å, 77 s), the buffer frozen and its descreening kept. The anchor rule chose backbone
    atoms 232, 231 and 230 and benzene's 978, 979 and 980: r₀ 6.941 Å, θ_A 103.1°, θ_B 65.2°.
    Each window runs at 0.5 fs in a 1 ps⁻¹ bath, discards 2 ps and samples 10 ps every 10 fs.
  - **The schedule was thinned from 3c-1's**, to 15 windows from 22: five for the restraint, two
    more for the charges and eight more for the van der Waals, because every overlap there was
    0.38–0.50. Each window recorded its energy at every candidate state between its neighbours,
    so that a window could be added between two that had run without running either again. One
    was to be added wherever an interval's overlap fell below 0.1, three times the 0.03 that
    Klimovich, Shirts and Mobley (*J. Comput.-Aided Mol. Des.* 29, 397 (2015), PMC4420631) find
    tolerable "with enough samples". **None fell below it**: the lowest is 0.315.

  | window | λ_r | λ_e | λ_v | → next: BAR (kcal/mol) | overlap | EXP fwd | EXP rev |
  | --- | --- | --- | --- | --- | --- | --- | --- |
  | 0 | 0 | 1 | 1 | +0.154 ± 0.043 | 0.481 | +0.185 | +0.075 |
  | 1 | 0.1 | 1 | 1 | +0.099 ± 0.008 | 0.498 | +0.100 | +0.096 |
  | 2 | 0.25 | 1 | 1 | +0.164 ± 0.018 | 0.492 | +0.143 | +0.230 |
  | 3 | 0.5 | 1 | 1 | +0.301 ± 0.031 | 0.478 | +0.282 | +0.321 |
  | 4 | 1 | 1 | 1 | −2.580 ± 0.006 | 0.494 | −2.575 | −2.585 |
  | 5 | 1 | 0.5 | 1 | −2.593 ± 0.006 | 0.494 | −2.592 | −2.591 |
  | 6 | 1 | 0 | 1 | +0.188 ± 0.025 | 0.386 | +0.135 | +0.192 |
  | 7 | 1 | 0 | 0.9 | +0.597 ± 0.042 | 0.315 | +0.554 | +0.645 |
  | 8 | 1 | 0 | 0.75 | +0.984 ± 0.040 | 0.331 | +1.116 | +1.005 |
  | 9 | 1 | 0 | 0.6 | +1.388 ± 0.047 | 0.366 | +1.499 | +1.394 |
  | 10 | 1 | 0 | 0.45 | +1.796 ± 0.109 | 0.396 | +1.807 | +1.827 |
  | 11 | 1 | 0 | 0.3 | +1.234 ± 0.072 | 0.419 | +1.429 | +1.168 |
  | 12 | 1 | 0 | 0.2 | +1.209 ± 0.106 | 0.417 | +1.142 | +1.283 |
  | 13 | 1 | 0 | 0.1 | +1.104 ± 0.137 | 0.379 | +1.015 | +1.077 |
  | 14 | 1 | 0 | 0 | | | | |

  | leg or segment | BAR (kcal/mol) | TI (trapezoid) |
  | --- | --- | --- |
  | solvent: charges off (4 seeds, spread 0.0003) | +2.3395 ± 0.0003 | +2.3403 |
  | solvent: van der Waals off (identically zero; see below) | exactly 0 | exactly 0 |
  | complex: restraint on | +0.718 ± 0.065 | +0.753 |
  | complex: charges off | −5.172 ± 0.011 | −5.171 |
  | complex: van der Waals off | +8.501 ± 0.301 | +8.491 |
  | **complex leg** | **+4.047 ± 0.309** | +4.073 ± 0.313 |
  | restraint released to 1 M, analytic (extended form −7.801) | −7.784 | |
  | **ΔG°_bind = +2.340 − 4.047 + 7.784** | **+6.08 ± 0.31** | |
  | nonpolar estimate, `0.005 × ΔA`, ΔA = −293.3 Å² at the start | −1.47 | |
  | **ΔG°_bind with it** | **+4.61 ± 0.31** | |
  | experiment, Mobley et al. 2007, Table 1 | −5.19 ± 0.16 | |

  - **Hysteresis.** The first half of every window's samples gives +3.63 ± 0.31 for the complex
    leg and the second half +4.42 ± 0.31: a drift of 0.79 ± 0.44, 1.8σ, which says that 2 ps of
    equilibration did not settle every window. EXP summed forward gives +4.24 and backward +4.14,
    against BAR's +4.05. A second seed of the complex leg, another 5 h, was not run. The
    autocorrelation times reach 37–48 samples at both ends of the path (λ_r = 0, and λ_v ≤ 0.3),
    so those windows hold 20–30 independent samples each and their σ̂ is itself uncertain by
    about a fifth.
  - **Eight of the solvent leg's eleven windows sample a segment that is identically zero**: with
    the whole molecule as the group nothing depends on `λ_v`. They are kept, so that the two legs
    run one schedule and the zero is seen on sampled configurations, at about 130 s of the leg's
    178.
  - **The solvent leg's cross-check.** Its charges' free energy, +2.3395, against benzene's
    ⟨ΔG_GB⟩ in the coupled state, −2.342; the linear-response average −(⟨G⟩₁ + ⟨G⟩₀)/2, +2.339;
    and Zwanzig's exponential average from the coupled end, +2.339. Benzene is nearly rigid and
    its GB energy barely moves, so all four agree to 3e-3 kcal/mol.
  - **The hydration free energy** under the model is −2.340 polar, and −1.133 with the nonpolar
    term (0.005 × 241.4 Å²). **FreeSolv's experimental value is −0.90 ± 0.20, not the −0.87 the
    brief gave**: `database.txt` v0.52 from `github.com/MobleyLab/FreeSolv` (master, fetched
    2026-10-05, SHA-256 `2d13f095…f260`) reads `mobley_3053621; c1ccccc1; benzene; -0.90; 0.20`,
    referenced to 10.1039/P29900000291. The model is 0.23 too negative, about the experiment's
    uncertainty.
  - **Against experiment the model is +11.3 kcal/mol off (+9.8 with the nonpolar term): it does
    not bind benzene.** Nothing is asserted. What limits the comparison, roughly by size:
    - **GB's electrostatic desolvation with QEq charges.** The charges cost +7.51 kcal/mol of
      binding (+2.34 in solvent, less −5.17 in the complex), where vacuum's charge segment gave
      −1.48. 2c-2 and A-2 found the same penalty at fixed pose (+12.21 in the 10 Å pocket at the
      crystal pose), with its cross terms unchecked against Poisson–Boltzmann or explicit water.
      QEq's charges are not the ones OBC II was fitted with.
    - **The cavity solvated as water when it is empty.** Removing benzene's volume lets GB
      solvate the cavity's walls, and L99A's real cavity is empty (A-2). The van der Waals segment
      is +8.50 against vacuum's +13.80. A-2 put this part of the polar term, the pocket
      desolvated by benzene's volume, at +4.31 kcal/mol in the 10 Å binding at the crystal pose.
    - **The frozen buffer and the cut.** The 10 Å binding's cut surface is solvated in every
      state (A-2: a cut pocket solvates its cut), and nothing outside 6 Å moves.
    - **UFF's van der Waals** for the cavity, and no nonpolar term in the sampled Hamiltonian.
    - **Sampling**: ±0.31 by BAR, a 1.8σ drift between halves, and one seed. That is small beside
      the rest.
  - **Wall time.** The solvent leg, 178 s. The complex leg, 5.11 h: the minimisation 77 s and the
    15 windows 18 318 s, 1157–1252 s each, 48.2–52.2 ms per step with the sampling, which
    evaluates up to five states at each of 1000 samples. 3c-2's cost test, run again after this
    change, reads 43.2 ms per evaluation and 44.6 ms per step without sampling: the partition's
    sums in the descreening cost nothing measurable.

  **Counts.** `cargo test -p pantometry-forcefield -- --list` counts **305 tests, twenty-five of
  them ignored**: 280 run by default. The twelve new ones are all in
  `benzene_bound_in_generalized_born.rs`: ten default, and the two measurements. 293 + 12 is the
  305, and 270 + 10 the 280. 3c-1's refusal test now holds that a solvated force field is accepted.

- **`pantometry_core::math::{exp, ln}`: an exponential and a logarithm of the workspace's own,
  the same bits on every platform, and generalized Born's step from 44.9 ms to 18.2.** 3c-2 left
  two ways past the platform's `exp` and `ln` for a decision; this is the first of them, taken.
  - **In the kernel**, as maths no domain owns, beside `vector` and `transform`: neither function
    knows anything about a physics, and a second domain should not have to depend on the
    forcefield to have them. Nothing in the kernel calls them, so no kernel result moved.
  - **Original code, under the crate's own `MIT OR Apache-2.0`, and no notice to carry.** The
    first version was ported from Arm's routines as musl has them and from fdlibm's `e_log.c`,
    and carried Sun's and the MIT notices; that would have made the kernel's licence a compound
    expression. It was **rewritten clean-room**: from the published method — P. T. P. Tang's
    table-driven `exp` and `ln` (*ACM TOMS* 15(2), 1989; 16(4), 1990), Cody–Waite reduction,
    Dekker's exact sum, Sterbenz's lemma, Taylor's series — without opening the port's code or
    tables, or any library's source. `math.rs` cites the mathematics and no code, and carries no
    third-party notice.
  - **The method.** `exp`: `x = (128k + j) ln 2/128 + r`, `n = 128k + j` rounded by adding
    `1.5 · 2⁵²`, `r` by a Cody–Waite split of `ln 2` whose head has 35 bits, so the head's product
    and difference are exact; `2^(j/128)` from a 128-entry double-double table; `eʳ − 1` by Taylor
    to `r⁵`; one rounding of `T_hi + (T_hi p + T_lo)`, then an exact `2ᵏ`. Below `2⁻¹⁰²²` the sum
    is added to 1 at `2¹⁰²²` times its size, so it is rounded **once** onto the subnormal grid.
    `ln`: `x = 2ᵏ m`, `m` in `[0.703125, 1.40625)`; `F` is `m` to eight fraction bits, `c` a
    nine-bit `1/F` from a 257-entry table (1 in the two cells about 1), and
    `r = m c − 1 = (F c − 1) + (m − F) c` — **every operation in it exact**, which the cells'
    width and `c`'s length are chosen for and the table's test checks cell by cell — then
    `k ln 2 − ln c` exact in its head, Dekker's sum, and Taylor to `r⁷`. **No division** in
    either. Only `+ − × ÷` and integer operations on the bits, in a written order, with no fused
    multiply–add. Every constant is written by `tools/math-reference/generate.py` (`--exp-table`,
    `--ln-table`, `--constants`) with mpmath at 50 digits.
  - **The reference tables are committed, and left out of the published kernel.** At 1.4 MB, they
    are the third exception to "nothing generated is committed", after `tools/presets` and
    `tools/parts`. `tools/math-reference/README.md` says why. Packaged with the first, 1.2 MB of
    them, they took `pantometry-core` from 147.4 KiB to 1.0 MiB compressed, for every user of the
    kernel. `Cargo.toml` excludes `tests/data/` and the one test that reads it, so a published
    crate holds no test that cannot compile, and CI and both gates still run the test from the
    repository. **The point recipe was extended for the rewrite's own seams** — `exp`'s table
    boundaries `(n + ½) ln 2/128` for `|n| ≤ 512`, `|x| = 708`, `ln`'s 256 cell edges at
    `k = −1, 0, 1` and `0.703125 · 2ᵉ` in every binade — appended after the old points, whose rows
    in the regenerated tables are byte for byte what they were: 64 993 and 76 308 points.
  - **fdlibm's `exp` and `ln` throughout were ported first, and measured worse on both counts**:
    0.8994 and 0.7526 ulp, and a step of 23.7 ms, against 4.4 and 3.4 ns a call in a loop. The
    gap between the loop and the step is what chose the table methods.
  - **Accuracy, against mpmath at 50 digits** (`tests/exp_and_ln_against_mpmath.rs`, against a
    committed table of the exact values' rounding at points the test draws itself by the script's
    recipe — whose input digest the table carries and the test checks first). The bounds are
    **derived from the method**, written beside them, and held at every point:

    | | bound | worst, 64 993 / 76 308 committed points | worst, ×20 random points (804 287 / 737 812), run once |
    | --- | --- | --- | --- |
    | `exp` | 0.520 | 0.5106 ulp | 0.5106 ulp at `x = 0.55506`, a boundary of `j`; 0.074% not correctly rounded |
    | `ln`, `c = 1` (`x` in `[1 − 3·2⁻¹⁰, 1 + 2⁻⁹)`) | 0.508 | 0.5000 ulp | 0.5010 ulp at `x = 1.0017` |
    | `ln`, `k = 0`, `c ≠ 1` | 0.508 | 0.5000 ulp | 0.5002 ulp at `x = 0.97946` |
    | `ln`, `k ≠ 0` | 0.501 | 0.5000 ulp | 0.5000 ulp; 6 points (0.001%) not correctly rounded, all of `ln` |

    On the committed points `ln` is correctly rounded at every one. The first version's were
    0.5058 and 0.5331 ulp against bounds of 0.52, 0.63 and fdlibm's 1. Each comparison also shows
    it can see one ulp: the values moved up one double are more than one ulp out somewhere.
  - **The constants are checked against closed forms summed in the test, not against the
    script**: `ln 2` by `Σ 1/(n 2ⁿ)` in double-double, its split, `128/ln 2`; every `2^(j/128)` by
    its Taylor series and by the exact products `T_j T_(128−j) = 2` and `T_j² = T_2j`; every
    `−ln c` by `2 atanh((c − 1)/(c + 1))`, each `c` the nearest nine-bit `1/F`, and per cell the
    two properties `ln` assumes: `|m c − 1| < 2⁻⁸` at both edges, and the exponent condition for
    Dekker's sum.
  - **Identities, each bound derived from the asserted ulp bounds**: `exp(ln x) = x` to
    `(0.508 |ln x| + 0.520) ε` relative, `ln(exp x) = x` to `(0.520 + 0.508 |x|) ε`,
    `ln(xy) = ln x + ln y` for `xy` exact; exactly `exp(±0) = 1` and `ln(1) = +0`; both monotone on
    consecutive doubles at every seam of the method, on 2¹⁷ doubles about 1 and on sorted grids
    across each domain.
  - **Special values**, as IEEE and `f64` have them, each tested by name: NaN, `±∞`, `±0`,
    negative `ln`, `exp`'s overflow at `0x40862e42fefa39ef` (finite there, `+∞` one double up)
    and underflow at `0xc0874910d52d3051` (the smallest subnormal there, `+0` one double down),
    subnormal arguments of both, subnormal results of `exp`, and `ln 2ᵉ` for every `e` from −1074
    to 1023.
  - **The bits are pinned, and it is a new digest**: FNV-1a of both at 8192 points each and their
    special values, `0x84532fb880d98a57`, the same here in debug and in release. The first
    version's, `0x150a0cfa04b62b07`, pinned another function over other points. CI runs it on
    Linux, macOS on arm64, Windows, release, and `wasm32-wasip1`; that matrix, not this machine,
    is what the claim rests on, and the rewrite has not yet been through it.
  - **The cost** (release, `x86_64-pc-windows-gnu`): 1.96–2.03 ns a call for `exp` against
    `f64::exp`'s 33.0–33.7, and 2.86–2.92 for `ln` against 17.4–17.6, in a loop over generalized
    Born's arguments. A first draft of `ln` that carried `(m − F)/F` as a product and its exactly
    computed remainder measured 3.67 ns; making `r` exact by construction instead removed the
    remainder. **The target of about 16 ms a step is not met**: on 3b's 987-atom complex one
    OBC II evaluation is 16.6 ms (43.1 before), 21.6 without the kept terms (57.0) and 33.6 by the
    direct sum (86.4); the radii 7.7 (18.1); a step of dynamics **18.2 ms against 44.9**, 11.4
    times vacuum's 1.6, 2.37 ns/day — all with the first version. The rewrite, timed in the same
    session as the first, is 18.48–18.49 ms against its 18.29–18.30: 1% slower a step, though
    `ln` alone is the same speed in a loop. With both functions replaced by arithmetic of no
    accuracy the step is 12.3 ms, so they still cost about 6 ms: 5.5–6 ns a call in place. Each
    sits in a chain of dependent divisions and square roots, and the step pays their latency, not
    their throughput.
  - **Adopted by generalized Born only**: `ln(U/L)` in `descreening`, `E` in Still's eq 3, and
    `e^(−κf)` in salt — the evaluation, the reference, and the alchemical evaluation alike, so
    3c-2's equivalence holds to the bit, with no test changed. **What moved**, the platform's
    against the first version, on aspirin (OBC II, and in 0.15 M salt) and on the 3b complex: no
    energy by a bit; no radius of aspirin's, and 20 of the complex's 987, by at most 2.95 ε;
    forces by at most 0.23 ε of aspirin's largest and 3.97 ε of the complex's. **And the rewrite
    against the first version**, the same systems: aspirin not a bit, energy, radius or force, in
    either model; the complex's energy not a bit, in OBC II or in salt; 18 of its radii by at
    most 2.95 ε of themselves; forces by at most 3.64 ε of its largest (3.96 in salt). Every
    forcefield test passes unchanged: 280, and 25 ignored.
    - **Not bit-identical yet**: GB's `tanh`, once per atom, is still the platform's, as are the
      QEq charges it is given.
    - **Not rerun**: 3b's 8 ps OBC II trajectory and 3c-3's free energy were measured with the
      platform's functions. A run now is another trajectory of the same model, since dynamics
      amplifies a few-ulp force, and the same within its statistics; `complex.rs` says so.
    - **The rest of the forcefield's `exp` and `ln` are not switched** — QEq, BAR, the Langevin
      coefficient, Boresch, UFF's bond-order term: none is hot, and every result they feed also
      passes through UFF's torsion `cos` (QEq through `powi`), so switching them would make
      nothing cross-platform and move their bits.
  - **Other crates that call `exp` or `ln` outside tests, for a later decision** (switching any
    would move its pinned and closed-form numbers): the kernel's `Rng::gaussian` (`ln`, with
    `cos`) and `Rng::poisson` (`exp`); `pantometry-em` `cavity.rs` (a Gaussian pulse);
    `pantometry-mechanics` (the damped overshoot); `pantometry-molecular` `fluid.rs` (a Langevin
    decay); `pantometry-optics` `diffraction.rs`, `propagation.rs` and `spectrum.rs` (Gaussians,
    and `exp_m1` for Planck); `pantometry-porous` (Arrhenius viscosity and kinetics, `puck.rs`);
    `pantometry-quantum` (a wave packet); `pantometry-thermal` `solid.rs` (`ln` in a view
    factor); `pantometry-view` `colour.rs` (`exp_m1` and a Gaussian). `pantometry-acoustic`,
    `-fluid` and `-pharmacokinetic` call them only in tests.
  - **Sabotage of the rewrite**: fourteen, each restored by copying the original back, `touch`,
    and SHA-256, with `--no-fail-fast`. Caught: `ln 2`'s tail off by 10⁻⁷ (the `exp` bound and
    the digest); `exp`'s `r⁵` term dropped; `ln`'s `r⁷` term dropped (only by the 0.508 bound, in
    `k = 0`); `−0` taken for a negative number by `ln`; `exp`'s Cody–Waite tail dropped; `ln`'s
    Dekker error dropped, and `r` computed as `m c − 1` (the 0.501 bound, the identities, the
    digest); subnormal results rounded twice (the 0.520 bound only); `exp`'s overflow branch
    scaled by `+∞`; `F` rounded down rather than to nearest (the 0.508 bound only); one bit of a
    `−ln c` flipped at run time (the bounds, monotonicity, the digest); and one low bit each of an
    `EXP_TABLE` and an `LN_TABLE` word — 2⁻¹⁰¹ relatively, invisible to any accuracy test, and
    caught by the table tests alone. **One passed and was understood**: bit 40 of `T_lo` for
    `j = 37` flipped at run time moves results by about 2⁻⁶⁵, a ten-thousandth of an ulp, and the
    table test reads the constant, not the run-time value; the same kind of flip in the constant
    is the one caught above.

  **Counts.** `cargo test -p pantometry-core -- --list` counts **139**, two ignored (the ×20 table
  and the timing): the seventeen new are eleven in `math`, four in `exp_and_ln_against_mpmath.rs`
  and two doc tests.
- **`pantometry-forcefield` in a periodic box: Ewald electrostatics, UFF's van der Waals cut off
  with its long-range correction, and a ligand decoupled inside the box.** This is step W1 of the
  four-step explicit-water track that +6.08 kcal/mol against experiment's −5.19 under OBC II
  opened: periodic boundaries and Ewald here, TIP3P with constraints (W2), benzene's hydration
  (W3) and the solvated complex (W4) to come. Nothing here is water yet. New modules `periodic`
  and `ewald`.
  - **The box is orthorhombic** (`PeriodicBox`): enough for a water box and a solvated protein,
    and its minimum image is one comparison per axis. No triclinic box.
  - **Molecules are made whole through their bonds at every evaluation**, breadth first from each
    molecule's lowest atom, and the bonded terms run on those positions unchanged. The other
    choice, the minimum image inside every bonded term, would have rewritten every angular term,
    which reads positions rather than displacements. The walk is `O(N)` and **changes no bit** when
    no bond crosses a face: aspirin's four bonded terms in a 40 Å box are the vacuum force field's
    to the bit. Molecular dynamics never wraps, so a trajectory stays continuous.
  - **Van der Waals is cut off at the Ewald cutoff, unshifted, with Allen and Tildesley's
    long-range correction** (eqs 2.144–2.145) for the energy and the pressure. UFF's
    `D[(x/r)¹² − 2(x/r)⁶]` is Lennard-Jones's `4ε[(σ/r)¹² − (σ/r)⁶]` exactly, with `ε = D` and
    `σ = x/2^⅙`, and **its geometric combination makes the double sum of the correction separable**,
    the square of `Σ √D_i x_i⁶` and of `Σ √D_i x_i³`, so it is `O(N)`. The virial is the model's
    own strain derivative, so its tail part is `E_tail`; A&T's `P_tail` differs from `E_tail/V` by
    the impulsive term of a potential that steps at `r_c`, `(2π/3V²) ΣΣ r_c³ u(r_c)`, which has
    the sign of `u(r_c)` and so makes `P_tail` the more negative for UFF; `pressure_correction`
    gives it for W2, and a test holds the difference to that term to 1e-12.
  - **Ewald**: real space by minimum image inside `r_c ≤ L/2`, reciprocal space over a sphere of
    wave vectors, the self term, the background `−π Q_net²/(2Vα²)`, and the 1-2 and 1-3 pairs'
    `erf(αr)/r` taken back out wherever they are. UFF keeps 1-4 pairs whole, so no pair is
    scaled. Analytic forces, and the virial `W_ab = −∂U/∂ε_ab` with the reciprocal stress
    `δ_ab − 2 k_a k_b (1/k² + 1/4α²)`. **α and `k_c` are chosen by Kolafa and Perram's (1992)
    estimates**, read secondarily in Saffar Shamshirgar, Hess and Tornberg (arXiv:1712.04718):
    the accuracy `δ` is each RMS estimate in units of `k_e Q/r_c`, so it does not depend on the
    charges. A thousand waters in a 31 Å box at `r_c` = 9 Å and `δ = 10⁻⁵` get α = 0.301 Å⁻¹ and
    1 051 wave vectors.
  - **Found: the estimate leaves out the larger error.** The reciprocal cut also drops the diagonal
    of `|S(k)|²`, which tends to `Q` at large `k`, so the energy is low on average by
    `k_e Q (α/√π) erfc(k_c/2α)` — derived here, `reciprocal_bias`. On a disordered box of 96 charges
    the measured error was that bias, 0.88–1.10 of it at every δ from 10⁻⁶ to 10⁻¹¹, and ten times
    Kolafa and Perram's spread. It is a constant per atom at fixed charges, so it moves no force,
    and it has no cross term, so it is not in a decoupling's `∂U/∂λ_e`. **Kolafa and Perram's
    real-space estimate is right**: with the reciprocal sum made exact, the RMS error over sixteen
    disordered boxes was 1.022, 1.061 and 1.054 of it at δ = 10⁻⁴, 10⁻⁶ and 10⁻⁸.
  - **`erfc` is this crate's own**, clean-room: Taylor's series of erf below 0.4375, and
    `exp(−x²) erfcx(x)` above, with `erfcx` a Chebyshev series on nine intervals to 6 and in
    `1/x²` beyond, every coefficient fitted with mpmath by `tools/erfc-reference/generate.py`, no
    library's code or table read; `x²` carried exactly by Dekker's product, and the kernel's `exp`.
    **Against mpmath at 50 digits, at most 3.21 ulp over 992 240 points** (2.74 from 0.4375 to 6,
    3.21 above 6, 1.22 below 0.4375, 1.05 for negative arguments), run once; the 438 points the
    test carries are held to 3.5 ulp, and their worst is 2.67. It is the same function on every
    platform. The reciprocal sum's phase factors still start from the platform's `sin_cos`, 3N
    calls an evaluation.
  - **`PeriodicForceField`** is UFF in the box — the bonded terms on whole molecules, every
    non-bonded pair by a cell list of cells at least half the cutoff wide, the correction and
    Ewald — and a `Potential`, so `MolecularDynamics` integrates it unchanged. It never builds the
    `N²` pair list `ForceField::new` keeps; `ForceField` gained a crate-private bonded-only build
    for it.
  - **`PeriodicDecoupling`** is `Decoupling` in the box, an `Alchemical` for `Windows`: the
    group's van der Waals pairs with the rest through the same soft core, and **its Coulomb pairs
    with the rest scaled by `λ_e` as the cross terms of the Ewald sum** — the real pairs across,
    `2 Re(S_rest* S_group)` in reciprocal space and the background's cross part, the self terms
    cancelling. **The group's own non-bonded terms are computed in vacuum, at every λ**, with no
    interaction with its own images, so the decoupled state is the rest in its box and the group
    alone in vacuum: `Decoupling`'s decoupled state exactly, whatever box each leg of a cycle used.
    The coupled state leaves the group's image interaction out too, a finite-size artefact: for
    two UFF waters in a 9.3 Å box, −0.25 kcal/mol, the long-range correction's share included. The cross part of the long-range correction is
    scaled linearly by `λ_v`.
  - **Checked against closed forms** (`tests/the_ewald_sum_against_closed_forms.rs`, point charges;
    `tests/a_molecule_in_a_periodic_box.rs`, UFF):
    - **Madelung constants** from periodic cells, at δ from 10⁻³ to 10⁻¹¹:

      | δ | NaCl, 64 ions | CsCl, 54 ions |
      | --- | --- | --- |
      | 10⁻³ | 1.749368 (1.8e-3) | 1.762345 (3.3e-4) |
      | 10⁻⁵ | 1.747635 (7.0e-5) | 1.762947 (2.7e-4) |
      | 10⁻⁷ | 1.7475663 (1.7e-6) | 1.7626744 (3.8e-7) |
      | 10⁻⁹ | 1.74756461 (2.0e-8) | 1.76267480 (2.4e-8) |
      | 10⁻¹¹ | 1.7475645949 (2.6e-10) | 1.7626747733 (2.5e-10) |
      | exact | 1.747564594633 | 1.762674773071 |

      **Kolafa and Perram's estimate does not hold for a crystal** — the error is 1.3 to 19 times
      it on rock salt and 0.2 to 17 on CsCl, because an ordered lattice puts its charge in Bragg
      peaks and shells — and is not asserted there. What is asserted is exact: every term the two
      cutoffs leave out, summed by brute force over the image pairs out to `r_c + 7/α` and the
      wave vectors to `3 k_c` and added back, gives both constants to 2.4e-15 or better at every δ.
    - **A lone charge in a cube is the Wigner lattice**, `k_e q² ξ/(2L)`, `ξ = −2.837297479480620`:
      with the omitted shell of wave vectors added back, to 8e-14 at three cutoffs and two
      accuracies — the background is what makes that so — and the bias estimate is 0.78–1.05 of
      that shell.
    - **The total does not depend on α**, from 0.23 to 0.99 Å⁻¹ over nine cutoffs and accuracies:
      each is the tightest's less the difference of their reciprocal biases, within 4σ of
      Kolafa and Perram's estimate (measured −3.4σ to +1.1σ) — the signed claim, not a bound on
      the distance.
    - **The non-periodic limit at its rate**: a neutral cluster in boxes of 12 to 32 Å is its own
      Coulomb energy less the tinfoil term `2π k_e μ²/(3V)`, and what is left times `L⁵` is
      constant to 2.8% from 16 to 32 Å (5.7% from 12), held to 5%. Aspirin with neutral charges in 40, 60 and 80 Å boxes: the energy's remainder
      falls by 7.618 and 4.219 against `L⁻⁵`'s 7.594 and 4.214, and the forces', after the dipole
      field `(4π k_e/3V) q_i μ`, by 7.630 and 4.221.
    - **Exclusions**: a cluster with six pairs excluded is the non-periodic energy without them to
      8e-7 of it, against a smallest excluded pair 4.5e5 times larger.
    - Forces against central differences and the virial against the strain derivative, each to a
      tolerance from the step's truncation and the energy's rounding — on point charges with a net
      charge and exclusions, and on 27 and 64 UFF waters, by brute force and by the cell list; the
      energy and forces unmoved, to 5e-16 of the parts, by a translation, by wrapping every atom on
      its own with 27 molecules cut by a face, and by scattering every atom over images up to four
      boxes away; the correction against its integral by quadrature to 1e-10, and for one type
      against A&T's printed forms to 1e-13; the soft-core decoupling's `∂U/∂λ_e` the three Ewald
      sums' cross terms to 1e-12, its end states the systems they say to 1e-12, and its λ
      derivatives and forces against differences; Langevin dynamics in the box the same bits cut
      into 24 calls, two or three.
  - **Sabotage: twenty-three, every one caught**, each restored by copying the original back,
    `touch` and SHA-256, the two test files run with `--no-fail-fast`: the self term's sign (six
    tests); the exclusion correction removed (three); the wave vectors without their 2π (nine);
    the pair search's minimum image never wrapping up (five), and the box's by `floor` (four);
    the background's sign (the Wigner lattice alone); the reciprocal force's sign (four); the
    reciprocal virial without `1/4α²` (both virial tests); the correction's 9 made 3 (its
    integral alone — the virial is the model's own derivative and stays consistent); the
    correction left out of the virial; `erfc` without the exact square, and its Clenshaw sum's
    first coefficient a part in 10¹⁵ off (both by the mpmath points alone); molecules not made
    whole (the wrapping test); the cell stencil shifted one cell (α-independence and the
    estimate's own test); the decoupling's background cross term, its soft-core λ gradient, the
    long-range cross term unscaled, and 1-3 pairs not excluded; Kolafa and Perram's real estimate
    without its `(α r_c)⁻²`; the real-space force without its Gaussian (six); the reciprocal
    Gaussian 1% wide (eight). **One passed first, and never reached its guard**: dropping the
    group's own pairs, because the decoupled group was one water, whose every pair is 1-2 or 1-3.
    With two waters it is caught — **and the two waters found a defect**: the group's own pairs
    were taken by minimum image, which for a group spanning more than half the box measures a
    pair through an image, 1.76 kcal/mol off the vacuum in the 9.3 Å test box. They are now taken
    on the whole molecules as given, and the minimum image put back is the twenty-third sabotage,
    caught.
  - **A review found five more that passed, and each now has the test that fails it** (each
    sabotage run before the test was added, passing, and after, caught; restored by copy, `touch`
    and SHA-256):
    - **The cell list was never compared with brute force**, and every box it ran on had `2L/r_c`
      whole, where a cell count rounded up comes out right: `floor` made `ceil` passed. Now 120
      charges in a 17.1 × 15.3 × 14.0 Å box at `r_c` = 5 Å, the real-space energy against every
      pair by minimum image, agree to 2.3e-17 of the pairs' magnitudes, and the sabotage is caught.
      The 27-cell branch of the cell list was unreachable — three cutoff-wide cells across is six
      half-cells — and is removed.
    - **No box was not a cube**, so the wave-vector extent, or the cell count, taken from `L_x` for
      every axis passed. Now rock salt's 2×2×3 and 3×2×2 supercells give its Madelung constant,
      completed by brute force to 1.0e-14, and both are caught. Their own truncation error is up
      to 103 times the estimate, the tightest 1.2e-9: `k_c` is set from `V^⅓`.
    - **The bonded virial was never compared with molecules cut by a face**: taking `Σ r ⊗ F` on
      the raw positions passed. The image test now compares the virial, to 4.5e-15 of its largest
      with 27 molecules cut, and catches it.
    - **No excluded pair was beyond `r_c`**: skipping their correction there passed, in the sum and
      in the force field. Now the energy with the exclusions less the energy without them is
      `−k_e Σ q_i q_j / r` over the excluded pairs (`erf(αr)/r` alone for a pair past `r_c`) to
      rounding, in a 40 Å box and in an 8 Å one at `r_c` = 1.5 Å with two pairs past it; and the
      force field's Ewald energy is the sum's at `r_c` = 1.4 Å, inside a water's H–H. That
      replaces a bound on the large-box remainder with a margin of 1.25.
    - **The off-diagonal virial was checked only for symmetry**: its sign flipped passed. Now
      `W_xy`, `W_xz` and `W_yz` are each the central difference of the energy under a simple
      shear, written in the test — the same pairs and integer wave vectors, `k_b → k_b − ε k_a`,
      the volume and every phase unchanged — in a non-cubic box with exclusions and a net
      charge, agreeing to 2.3e-8 of each or better.
    - A&T's `P_tail` is now held to `E_tail/V` plus the impulsive term written in the test, to
      1e-12, where it was only required to be the lower.
  - **The cost** (release, one core, `x86_64-pc-windows-gnu`; water's density, `r_c` = 9 Å,
    point charges with the waters' exclusions):

    | atoms | box (Å) | δ = 10⁻⁵: waves, ms (real, reciprocal) | δ = 10⁻⁶: waves, ms (real, reciprocal) |
    | --- | --- | --- | --- |
    | 375 | 15.5 | 257, 2.0 (1.5, 0.4) | 522, 2.2 (1.5, 0.8) |
    | 1 029 | 21.7 | 423, 8.5 (6.8, 1.7) | 895, 10.1 (6.7, 3.4) |
    | 3 000 | 31.0 | 1 051, 37.7 (25.8, 11.9) | 2 192, 51.1 (25.4, 25.7) |
    | 6 591 | 40.4 | 1 955, 104 (56, 48) | 4 300, 160 (55, 105) |
    | 12 288 | 49.7 | 3 309, 249 (97, 152) | 7 385, 434 (95, 338) |
    | 24 000 | 62.1 | 5 533, 697 (198, 499) | 12 838, 1 347 (198, 1 149) |

    The whole force field on the 3 000-atom box is 40.2 ms at 10⁻⁵ and 52.7 at 10⁻⁶. Real space is
    linear in `N`; the reciprocal sum goes as `N^1.8` here and as `N²` in principle at a fixed
    cutoff. A first cell list, cells a whole cutoff wide and the minimum image by `round`, cost
    78 ms of real space at 3 000 atoms where this one costs 26.
  - **PME is needed for W4, not for W3.** At W4's 20 000–30 000 atoms one evaluation is 0.7–1.3 s,
    most of it reciprocal and growing as `N²`: a nanosecond at a 2 fs step would be four to eight
    days on one core. W3's few thousand atoms cost 40–50 ms an evaluation, a nanosecond in six to
    seven hours, where real space is half and a neighbour list with a skin, not PME, is the next
    saving. Smooth PME (Essmann et al. 1995) is left for later, and this sum is what it will be
    checked against.
  - **Not here**: PME, a neighbour list, a triclinic box, a barostat, constraints, QEq or
    generalized Born in a box.

  **Counts.** `cargo test -p pantometry-forcefield` passes **305**, 27 ignored (280 and 25 before):
  fourteen in `the_ewald_sum_against_closed_forms.rs`, with the million-point dump for
  `generate.py --measure` ignored, and eleven in `a_molecule_in_a_periodic_box.rs`, with the timing
  ignored. Thirty sabotages in all, every one caught. `tools/erfc-reference/README.md` says what
  the script prints and where it is committed.

- **`pantometry-forcefield` has water: rigid TIP3P held by SETTLE inside BAOAB, a box of it, and
  the liquid at 298 K set beside the literature.** This is step W2 of the four-step explicit-water
  track: W1 built the box and Ewald, W3 is benzene's hydration free energy in this water and W4
  the solvated complex. New module `water`.
  - **The model is TIP3P as Jorgensen et al. 1983 published it** (*J. Chem. Phys.* 79, 926, Table
    I): r(OH) 0.9572 Å, ∠HOH 104.52°, q(O) −0.834, q(H) +0.417, and an O–O Lennard-Jones term
    `A/r¹² − C/r⁶` with A = 582.0 × 10³ kcal Å¹² mol⁻¹ and C = 595.0 kcal Å⁶ mol⁻¹, so ε =
    0.152 073 kcal/mol and σ = 3.150 656 Å. **The hydrogens have no Lennard-Jones term**: CHARMM's
    modified TIP3P, which gives them one, is another model. **Read secondarily**: the paper could
    not be opened (the publisher and the one copy found refused). The row is as Table 1 of Izadi,
    Anandakrishnan and Onufriev (*J. Phys. Chem. Lett.* 5, 3863 (2014), read in arXiv:1408.1679)
    and Wikipedia's "Water model" article reproduce it, which agree to every digit, and LAMMPS's
    TIP3P page gives the same geometry, charges and ε and σ.
  - **The charges sum to zero exactly**, not to rounding, and the dipole is the closed form
    `2 q_H r_OH cos(θ/2)` = 0.488 56 e Å = 2.347 D (Izadi et al. tabulate 2.348).
  - **Water and a UFF solute mix by UFF's own rule**, geometric in the distance and the depth: the
    oxygen is a UFF atom with `x = 2^⅙ σ` and `D = ε`, a hydrogen one with `D = 0`. **A judgement and
    an approximation**: neither parameter set was fitted with the other. For it: it is the rule the
    solute's own pairs use, and it keeps W1's long-range correction separable, `O(N)`. Against:
    GAFF, AMBER and Mobley's benchmark combine σ arithmetically. For benzene's aromatic carbon
    against a water oxygen the two rules differ by 0.09% in the distance and not at all in the
    depth. `PeriodicForceField::tip3p` is water alone, `PeriodicForceField::solvated` a UFF solute
    followed by waters, with `with_solute_charges` and `rigid_waters`. `PeriodicForceField::bonded`
    is now an `Option`, `None` for water alone (an unreleased W1 signature).
  - **SETTLE** (Miyamoto and Kollman, *J. Comput. Chem.* 13, 952 (1992)) is **derived here from the
    conditions it solves, not transcribed**: the constraint impulses are a combination of the old
    bond vectors, so they lie in the old plane, sum to zero and exert no torque. That fixes the
    two out-of-plane angles from the out-of-plane coordinates, the centre of mass, and the in-plane
    turn from `α cos θ + β sin θ = γ`, of whose two roots the smaller turn. The velocities are made
    RATTLE-consistent by the mass-weighted projection onto each water's rigid motions, a 3×3 solve.
    SETTLE's cube roots and the box's are by bisection in `+ − × ÷`, not the platform's `cbrt`.
  - **Inside BAOAB** (`MolecularDynamics::with_constraints`): SETTLE after each drift, with the
    velocity given the same correction `Δq/h`; the projection after each kick and after `O`, which
    takes an isotropic Gaussian to Maxwell–Boltzmann on what is left; the bath's work booked after
    the projection. Velocity Verlet becomes RATTLE. **Six degrees of freedom a water.** Frozen
    atoms still work: a water frozen whole is not SETTLEd and keeps its bits, one frozen in part
    is refused by name. `thermalised` projects its draw. `half_step_velocities` gives the
    velocities after `O`, from which translation and rotation are read. **With no constraints
    every new step is skipped and a run is the bits it was**: the whole forcefield suite passes
    unchanged.
  - **The box** (`WaterBox::lattice`): `k³` waters on a simple cubic lattice at 0.997 g/cm³, each
    turned by a uniform rotation — the unit quaternion of four normals of `Rng::for_index(seed,
    water)` — and `WaterBox::equilibrate` melts it: 0.5 ps at 0.5 fs in a 50 ps⁻¹ bath, then at 2 fs
    in a 5 ps⁻¹ one.

  **Checked against closed forms and invariants** (`tests/a_rigid_water_against_closed_forms.rs`,
  thirteen tests, 17 s unoptimised):
  - Table I as it stands; ε and σ against LAMMPS's printed 0.1521 kcal/mol and 3.1507 Å to their
    last digit, a third source for the transcription (C = 600 in both places, which the
    derivation alone cannot see, now fails); `4εσ¹²` and `4εσ⁶` give back A and C to 1e-14; the oxygen as a UFF pair is
    `A/r¹² − C/r⁶` at six distances; the charges' sum is `0.0`; the triangle's lengths, angle and
    centre of mass; the dipole against its closed form to 1e-14 and against 2.35 D.
  - The box's force field written out: the O–O pairs by hand to 1e-13, A&T's one-type tail to
    1e-13, the electrostatics the Ewald sum with each water's three pairs excluded, no bonded term;
    and methane among eight waters, its pairs with each oxygen by the geometric rule to 1e-12.
  - **SETTLE against SHAKE**, iterated here to its fixed point, on 512 waters given a step's
    worth of random displacement: the three lengths to 1.56 ε × the water's largest coordinate
    (bound 16), the centre of mass kept, and SHAKE's answer to 16.0 of the same unit (bound 64),
    against constraint corrections at least 2.6 × 10¹¹ of it. **The bound is in ε × the
    coordinate, not ulps of the bond**: a position is stored to half an ulp of itself, and the
    first bound, 8 ulps of the bond, failed at 8.5 for nothing but the waters' distance from the
    origin.
  - Projected velocities: nothing along a constraint to 8.1e-16, idempotent, rigid motions kept.
  - **Constrained NVE's RMS energy error falls as `h²`**, on eight waters in vacuum with no cutoff,
    written out in the test so that no cutoff's steps enter: measured over sixteen starts, the 32
    ratios per halving from 1 to 0.5 to 0.25 fs were 3.911–4.057, and the band is 3.8–4.2. **From
    2 fs to 1 fs they were 3.36–4.09**: one start in sixteen reads 3.36 there, 2 fs being past where
    `h²` alone describes it, so the test halves from 1 fs. The largest departure rather than the
    RMS was tried first and spread more (3.20–4.59 on the same start).
  - **Equipartition for rigid bodies**: in a 20 ps⁻¹ bath the 27 waters' translation read
    1.001 ± 0.026 and their rotation 1.022 ± 0.027 of `(3/2) N k_BT` each, each within 4σ, and
    `degrees_of_freedom` is `6N`. **The first version equilibrated 0.1 ps and read 1.057**
    (z = +2.4): measured over 20 ps in release, 0.1 ps left translation at 1.0143 ± 0.0068, and 1 ps
    brought four starts to 0.985–1.002, so the test equilibrates 1 ps.
  - **Equipartition to a part in a thousand**: 27 free rigid waters under no potential, in a bath
    with `γh = 1`, so that the samples are nearly independent (τ = 0.65 over 40 000): translation
    0.9985 ± 0.0009 and rotation 1.0005 ± 0.0009 of `(3/2) N k_BT`, held to 4σ.
  - The books balance under a bath with constraints: 5.1e-3 kcal/mol over 800 steps of 0.5 fs,
    against a bound of ten times constrained NVE's RMS error from the same start at the same step
    (5.06e-3), stated as a constant.
  - NVE's degrees of freedom are `6N − 6` with the rigid motion removed, asserted in every NVE
    run.
  - A run is the same bits in one call or 24, two or three; every water rigid after it; and the
    whole-step velocities have nothing along a constraint, to 2.2e-16, and so have the velocities
    after `O` at the positions they were projected at (`half_step_positions`, new), to 5.2e-16.
  - Frozen waters keep their bits with `6 N_free` degrees of freedom; a water frozen in part
    panics.
  - The lattice is at its density to 8ε, every centre of mass on its site, its orientations the
    seed's, and uniform: the bisectors' mean and fourth moment over 4 096 waters within 4σ of 0 and
    1/5.
  - W1's guarantees with TIP3P: forces against central differences, and the energy and forces
    unchanged when every atom is moved by its own whole number of box lengths.

  **Sabotage: seventeen, every one caught**, each restored by copying the original back, `touch` and
  SHA-256, the file run with `--no-fail-fast`: SETTLE's `sin ψ` sign and its other root (each by
  four tests, SHAKE among them); the degrees of freedom without the constraints (two); **no
  projection after `O` (the books alone: the noise along the bonds at 1 fs is 4% of the
  rotation's share, inside its 4σ)**; the charges swapped (five); no velocity correction after
  SETTLE (two); the work booked before the projection (the books); rotations from a quaternion
  uniform in a cube (the fourth moment); the molecule's mass one hydrogen short (five); the
  hydrogens given the oxygen's Lennard-Jones (four); the projection with the oxygen's mass on a
  hydrogen (four); `thermalised` not projecting (two); the drift never SETTLEd (three); the frozen
  mask not rechecked (two); the oxygen's `x` taken as σ (two); one O–H not excluded (one). **One
  passed first: no projection after the second kick.** It follows the same trajectory — the next
  step's first kick projects at the same positions, and projection is linear and idempotent — and
  only the whole-step velocities are wrong, whose kinetic energy error is `O(h²)` too. The cut
  test now holds those velocities to the constraints, and catches it.

  **A review found six gaps, and each now has the test that closes it**, every sabotage run
  against the test file before the fix and after, restored by copy, `touch` and SHA-256:

  | sabotage | before | after |
  | --- | --- | --- |
  | the constrained bath's noise 4% wide | passed (rotation z = +3.6 in the reviewer's run) | caught by the free-water equipartition |
  | the noise 2% wide, 4% in the energy | passed | caught: translation 1.0389 ± 0.0009, z = +42 |
  | NVE counting nine degrees of freedom a water | passed all twelve | caught by the `6N − 6` assertion |
  | the projection after `O` at the positions before the drift | caught by the books at 1e-3 of traffic, 2–4× margin | caught: 0.67 kcal/mol against the new 0.051 bound |
  | no `Δq/h` velocity correction | caught | caught: 1.36 kcal/mol, and by the `h²` test |
  | `Δq/h` halved | caught | caught: 0.68 kcal/mol, and by the `h²` test |
  | no projection after `O` | caught by the books alone | also caught by the velocities after `O` and by the free-water equipartition (rotation 1.86) |
  | `C = 600` in the code and in the test's own Table I | passed | caught by LAMMPS's ε and σ |

  - **The books' first new bound grew with the integrator's error**, ten times NVE's RMS error
    measured again in the same test, and a missing `Δq/h` passed it at 0.2 of itself because
    its NVE error grew too. The bound is now a stated constant from the correct integrator's
    measurement.
  - **A sub-ulp tolerance**: the UFF form of the O–O pair was held to `1e-14 A/r¹²`, a tenth of an
    ulp of the dominant `C/r⁶` at 9 Å, and passed by rounding's luck; the scale is now
    `A/r¹² + C/r⁶`.

  **The liquid** (`tests/liquid_water_against_the_literature.rs`, ignored, release). 216 waters in
  an 18.64 Å box and 512 in a 24.86 Å one, `r_c` = 9 Å, δ = 10⁻⁵ (309 and 606 wave vectors);
  melted 20.5 ps, then NVT at 2 fs in a 1 ps⁻¹ bath (100 and 50 ps) and 100 ps of NVE from its
  last state. **Reported, not asserted**; errors from each series' own autocorrelation, or from
  independent 25 ps blocks:

  | | 216 | 512 | literature |
  | --- | --- | --- | --- |
  | ⟨T⟩ after `O` (K) | 298.3 ± 1.2 | 298.3 ± 1.2 | 298 |
  | U/N (kcal/mol) | −9.610 ± 0.014 | −9.612 ± 0.015 | −9.86, Jorgensen 1983 (Monte Carlo, NPT, spherical cutoff), **as quoted, not read** |
  | g_OO first peak (Å) | 2.770 ± 0.004 | 2.782 ± 0.006 | 2.77, Izadi et al. 2014 Table 3 |
  | its height | 2.72 ± 0.02 | 2.73 ± 0.02 | not tabulated there |
  | P (bar) | +180 ± 37 | +139 ± 30 | −2.9 ± 2.5 at 33.00 nm⁻³, N = 256 (Yeh and Hummer 2004); here 33.33 |
  | NVE ⟨T⟩ (K) | 304.0 ± 0.4 | 293.3 ± 0.2 | |
  | D, NVE (10⁻⁵ cm²/s) | 5.80 ± 0.27 | 5.04 ± 0.13 | D_PBC 5.123 ± 0.027 (N = 256) and 5.315 ± 0.014 (512), 298 K (Yeh and Hummer) |
  | D under the bath | 5.67 ± 0.46 | 4.79 (two blocks) | 5.06 ± 0.09 NVT, 5.19 ± 0.08 NPT (Mahoney and Jorgensen 2001, 9 Å cutoff, 267 waters) |
  | k_BTξ/(6πηL) at the NVE run's ⟨T⟩ | +1.10 | +0.80 | |
  | D∞ from NVE | 6.90 ± 0.27 at 304.0 K | 5.84 ± 0.13 at 293.3 K | 6.05 by fit, 6.11 corrected, 298 K (Yeh and Hummer); 5.5 (Izadi et al.) |
  | D∞ under the bath | 6.75 ± 0.46 at 298.3 K | 5.60 (two blocks) at 298.3 K | |
  | experiment | | | 2.30 |

  - **The correction** uses Yeh and Hummer's own TIP3P viscosity, 0.308 mPa s, and ξ = 2.837297,
    with `k_BT` at each run's own mean temperature. **η's own temperature dependence is not
    included**: η is their 298 K value throughout, which understates the correction of the warmer
    run and overstates the cooler one's. (The table first used 298 K for both NVE runs; those
    D∞ were 6.88 and 5.85. The values here are the same runs, recomputed by the formula the test
    now uses; the runs were not repeated.)
    Yeh and Hummer and Mahoney and Jorgensen were read in primary; Izadi et al. in primary for
    their table, which cites others for its TIP3P column. The g_OO height is not compared with
    anything: no readable source that tabulates it was found, and a search snippet giving 2.71
    from an unidentified paper is not a citation.
  - **U/N against −9.86**: Jorgensen's is a cutoff Monte Carlo at the model's own density (0.980
    in Izadi et al.'s table); this is Ewald at 0.997 and its reciprocal bias, 4.9e-3 and 7.1e-3
    kcal/mol a water, is not put back. Izadi et al.'s ΔH_vap of 10.26 kcal/mol implies about
    −9.67 with `RT` taken off. The two box sizes agree to 0.002.
  - **The pressure is positive because the density is the experiment's**: TIP3P's own at 1 atm is
    0.980 (Izadi et al.), and the box is 1% denser than Yeh and Hummer's, which TIP3P's
    compressibility (57.4 × 10⁻⁶ bar⁻¹, Izadi et al.) makes about +175 bar.
  - **Each NVE run is at its own temperature**, 304.0 K and 293.3 K from baths at 298.3: an NVE run
    starts from one draw of the canonical total energy, whose spread, `k_BT √(N C_v/k_B)` with
    TIP3P's `C_p` of 18.74 cal/(K mol) (Izadi et al.), moves the mean temperature by 6.7 K at 216
    waters and 4.4 K at 512, one standard deviation. Neither drifts (−0.002 and +0.003 kcal/mol/ps
    over the box). The NVE temperature is now the library's own count, `temperature()` with
    `6N − 3` after `without_net_momentum` (new), not one written in the test. Water's D rises about
    2% a kelvin, so each NVE D belongs to its own temperature, and **the two sizes' D∞, 6.90 at
    304 K and 5.84 at 293 K, bracket Yeh and Hummer's 6.05–6.11 at 298 K** rather than test the
    size correction. That would need both at one temperature, longer.
  - TIP3P's excess over experiment's 2.30 is there, two and a half times: the model, not this code.

  **The cost** (release, one core, `x86_64-pc-windows-gnu`; rigid TIP3P at 2 fs, classical Ewald at
  δ = 10⁻⁵, `r_c` = 9 Å):

  | waters | atoms | wave vectors | ms a step | one evaluation | SETTLE + one projection |
  | --- | --- | --- | --- | --- | --- |
  | 216 | 648 | 309 | 7.03 | 6.77 | 0.029 |
  | 512 | 1 536 | 606 | 26.2 | 25.3 | 0.072 |
  | 1 000 | 3 000 | 1 051 | 55.2 | 53.5 | 0.136 |

  The constraints are a quarter of a percent of a step; the force field is the rest.
  - **W3's projection.** Benzene in 500 to 1 000 waters, about 20 windows of 300 ps at 2 fs, is 3
    million steps: **about 22 h at 512 waters and 46 h at 1 000**, before the vacuum leg, which is
    negligible. **But benzene's C–H bonds are UFF's and free**, and the step they allow is 0.5 fs
    (3a's measurement): at that step the same campaign is 87 and 184 h. So W3 needs one of a
    constrained solute X–H (SHAKE on the solute, which SETTLE does not do), a multiple time step,
    or a rigid benzene; and `Windows` does not yet take constraints, which is W3's to add. A
    neighbour list with a skin, not PME, is the next saving at this size (W1).
  - **Not here**: SHAKE for a solute, a barostat, PME, H–H or O–H Lennard-Jones (CHARMM's
    TIP3P), other water models, a measured viscosity.

  **Counts.** `cargo test -p pantometry-forcefield` passes **319**, 30 ignored (305 and 27
  before): thirteen in `a_rigid_water_against_closed_forms.rs`, one in `water`'s unit tests, and three
  ignored in `liquid_water_against_the_literature.rs` (two liquids and the timing).

- **`pantometry-forcefield` holds a solute's bonds to hydrogen, decouples it from rigid water in
  the box, and benzene's hydration free energy in TIP3P comes out at −1.57 ± 0.14 kcal/mol,
  against experiment's −0.90 ± 0.20.** This is step W3 of the four-step explicit-water track: W1
  built the box and Ewald, W2 the water, and W4 is the solvated complex. New module `shake`.
  Nothing is asserted against experiment.
  - **`Shake`: SHAKE for the positions and RATTLE for the velocities**, for any set of bond
    lengths (Ryckaert, Ciccotti and Berendsen 1977; Andersen 1983 — cited for the method, not
    opened; each condition derived here). The impulse on a bond lies along its old vector, and
    **each bond's quadratic is solved exactly**, by its smaller root written without cancellation,
    then swept Gauss–Seidel until every bond is within the stated tolerance, `|r² − d²| ≤ 2 tol d²`
    with `tol` = 10⁻¹², two hundred times above the rounding a 25 Å box allows. **Bonds that share
    no atom — benzene's six C–H — are solved in one sweep, to rounding.** The velocities are
    RATTLE's projection, `κ = r·u/(w|r|²)` per bond, whose fixed point is the mass-weighted
    orthogonal projection SETTLE makes in closed form. A frozen atom has inverse mass zero, so a
    bond to it holds its partner at the length; a bond frozen at both ends is skipped.
    `Shake::to_hydrogen` holds every bond to a hydrogen at UFF's natural length `r_IJ`, 1.0814 Å
    for benzene's C–H. **For bonds that share no atom the constrained distribution needs no Fixman
    correction**: the mass-weighted Gram matrix of their gradients is a constant diagonal.
  - **`MolecularDynamics::with_bond_constraints`** composes it with BAOAB where SETTLE is: SHAKE
    and `Δq/h` after each drift, each atom corrected once by its whole displacement; the projection
    after each kick and after `O`. It combines with SETTLE as long as no atom is in both, which is
    refused by name. **One degree of freedom a held bond** with an atom free to move. Without bond
    constraints every new step is skipped: the whole suite's earlier tests pass unchanged.
  - **`Windows::with_dynamics` and `Window::with_dynamics`** build windows from a template
    `MolecularDynamics` — masses, frozen atoms, rigid water, held bonds — and `Windows::new` is
    that with an unconstrained template, to the bit.
  - **`PeriodicDecoupling`'s couplings come from the cross terms alone**: the group's pairs with
    the rest inside the cutoff, by the evaluation's own wrapped minimum image, the reciprocal sum's
    cross terms, the background's and the long-range correction's, then each state's soft core.
    A sample at 29 states costs 4.3 ms against 25.4 for one evaluation of the 1 533-atom box;
    `coupling` is `couplings`' entry to the bit, and a unit test holds both to the evaluation's
    own λ-dependent terms to 64 ε. `PeriodicDecoupling::cross_dispersion_correction` gives the
    correction's cross part. **W1's treatment of the group's own images is kept**: benzene is
    neutral, so there is no Wigner term, and the image interaction the coupled state leaves out
    measured −0.008 kcal/mol at the start, with a dipole of 0.021 e Å whose tinfoil term is
    2 × 10⁻⁵ kcal/mol. The decoupled state is the water in its box and benzene alone in vacuum
    with its intramolecular terms, so `ΔG_hyd = −ΔG_decouple` with no vacuum leg.
  - `WaterBox::without_overlaps` takes out the waters a solute overlaps, by minimum image.

  **Checked against exact identities** (`tests/constrained_bonds_against_closed_forms.rs`, eleven
  default tests in 1.8 s unoptimised, and the band's measurement ignored):
  - benzene's C–H held to 0.67 ε of the larger coordinate (bound 16) after every step at 2 fs
    beside rigid water, the whole-step velocities and those after `O` tangent to 0.75 ε (bound the
    larger of 16 ε and RATTLE's 10⁻¹²);
    the six bonds are the bonds to hydrogen, at the stretches' natural lengths, and independent;
  - twelve bonds sharing atoms — the ring as well — within 1.01 tol of their lengths and of
    tangent after every step, and their NVE energy error falling by 4.06 per halving from 1 fs;
  - the degrees of freedom `3N − 3N_w − N_b`, less six in NVE, with a bond to a frozen atom
    costing one and a bond frozen at both ends none; **a bath fills exactly that many**: with no
    potential, at `γh = 1`, each C–H dumbbell's translation 0.998 ± 0.004 of `(3/2) k_BT`, its
    rotation 1.000 ± 0.004 of `k_BT`, the waters 1.000 ± 0.001 and the temperature through
    `degrees_of_freedom` 1.000 ± 0.001; and held bonds alone 1.001 ± 0.003 of `15 k_BT`;
  - **constrained NVE at 2 fs falls as `h²`**: benzene and the ten waters nearest it, rigid, in a
    30 Å box with `r_c` = 14 Å so that nothing crosses the cutoff. Over sixteen starts the ratio
    per halving was 4.229–4.508 from 2 fs, 4.050–4.109 from 1 fs and 4.012–4.027 from 0.5 fs,
    above 4 because the next term is `h⁴` and positive; the bands are 3.9–4.8, 3.8–4.2 and
    3.97–4.05. **And the excess over 4 shrinks by four per halving**, as the `h⁴` term's must:
    `(r₃ − 4)/(r₂ − 4)` read 0.2365–0.2593, held to 0.15–0.35. The ratios alone could not see an
    error of lower order with a small coefficient, whose excess grows instead (below). The RMS
    error at 2 fs was 0.066–0.48 kcal/mol over 0.2 ps. **This test is the whole guard of SHAKE's
    `Δq/h`**: removing it passes every Langevin and equipartition test and the hydration file;
  - a frozen carbon holds its hydrogen at the bond's length, keeping its bits; a run with SETTLE and
    SHAKE together the same bits in one call or cut three ways; an atom held twice refused.

  (`tests/benzene_hydrated_in_tip3p.rs`, ten default tests in 7.3 s unoptimised, three ignored):
  - **the end states are the systems they say** at a configuration constrained dynamics left:
    coupled, the solvated field less benzene alone in the box plus benzene in vacuum; decoupled,
    the water alone plus benzene in vacuum, energy and every force, to 10⁻¹² of the parts; and
    the cross correction written out pair by pair to 10⁻¹³;
  - `∂U/∂λ_e` and `∂U/∂λ_v` against central differences of the whole energy at four states;
  - every state of `couplings` its own `coupling` to the bit, and their differences the whole
    energy's to 64 ε; and the same **with a Boresch restraint** (benzene to two waters, 6.62
    kcal/mol stretched after the dynamics): `∂U/∂λ_r` the restraint's energy to the bit and the
    whole energy's central difference in `λ_r`, the branch W4's complex leg needs;
  - **a charge held at a distance has a closed form.** An ion held by SHAKE 3.3 Å from a frozen
    oxygen carrying the opposite charge, in a 9 Å cube: `ψ − 1/r − 2πr²/3V` is harmonic in the
    ball of radius `L`, so its sphere average is its value at the centre, Wigner's `ξ/L`, and
    `⟨E⟩ = k_e qQ (1/d + ξ/L + 2πd²/3L³)`. A 32 × 64 product rule gives −6.329688804 kcal/mol
    against −6.329688808, inside four times Kolafa and Perram's estimate; the directions span
    −8.47 to −5.05. **The windows decouple it to its quadrature**: BAR +7.090 ± 0.071 and TI
    +7.087 ± 0.070 against the exact +7.094 and the trapezoid's +7.091, where `−⟨E⟩` alone would
    be +6.330 — 0.764 apart, against the 8σ (0.564) that two four-σ windows need not to touch;
    the decoupled window's `⟨∂U/∂λ_e⟩` −6.336 ± 0.081 against the uniform −6.330;
  - a solute takes out exactly the waters brute force over 27 images says it overlaps;
  - constrained windows the same bits however they are cut, `Windows::new` the unconstrained
    `with_dynamics` to the bit, and the measurement's own sampler the windows' to the bit.

  **Sabotage: twenty-one at first, every one caught**, each restored by copying the original back, `touch`
  and SHA-256, the lib's unit tests and both files run with `--no-fail-fast`. SHAKE along the new
  bond (the `h²` test alone); RATTLE's sign on one atom (six); a frozen atom's inverse mass kept
  (two); the tolerance without `d²` (nine); RATTLE stopping after a sweep (the shared bonds);
  `to_hydrogen` holding every bond (six); the degrees of freedom without the bonds (three), and a
  bond to a frozen atom not counted (one); the projection without RATTLE (four); no `Δq/h` for
  held bonds (the `h²` test alone); the couplings without the background (the sphere's closed form
  alone), without the correction in `∂U/∂λ_v` (two), with the soft core at `λ_e` (three), without
  the minimum image (three) and past the cutoff (three); the cross correction halved (the end
  states); `with_dynamics` dropping the template's constraints (three) and seeding every window as
  window 0 (the sampler test alone); the overlaps without the minimum image (one). **Two passed
  first**, and each now has the test that fails it: `Δq/h` applied twice to an atom in two bonds
  (caught by the shared bonds' energy, `h²`), and a bath with held bonds and no rigid water skipping
  the projection after `O` (caught by held bonds alone in a bath).

  **A review found four more that passed, and each now has the test that fails it** (each run
  before the fix and after, restored by copy, `touch` and SHA-256):

  | sabotage | before | after |
  | --- | --- | --- |
  | SHAKE's `Δq/h` × 0.99 | passed: ratios 4.375 and 4.177, inside both bands | caught by the excess's shrinking (0.491–1.322 over sixteen starts against the band 0.15–0.35) |
  | `Δq/h` × 0.95 | passed the 2 → 1 fs band | caught by the same, and by the shared bonds' energy |
  | `l.electrostatics × U_B` for `l.restraint × U_B` in the couplings | passed: no test built a restraint in the box | caught by the restraint test |
  | `∂U/∂λ_r` halved | passed | caught by the restraint test |
  | `DEFAULT_TOLERANCE` 10⁻⁶ | passed: the shared bonds were held to the tolerance read back from the code | caught by four tests: the bound is now the literal 10⁻¹² |
  | a velocity that is not finite iterated to the sweeps' end | panicked "did not converge" | refused by name, a unit test |

  The tangency bounds were 16 ε, while RATTLE skips a bond already within `10⁻¹² |r||u|` of
  tangent: they are now the larger of the two, and the measured worst stays 0.75 ε. The "test can
  tell" margin of the held charge was 10σ with 8% to spare and no reason given; it is now 8σ, the
  separation two four-σ windows need, with a margin of 1.35.

  **The measurement** (ignored, release, one core; written window by window outside the
  repository, run once, exit 0 after 33 564 s).
  - **Before it ran**, the cost: 25.7 ms a step for benzene in 507 waters, 25.4 of it the force
    field and 0.015 the projections, and 4.3 ms for the couplings at all 29 candidates; a window
    0.79 h, thirteen 10.3 h, below the 30 h at which the run would have been referred back. It took
    9.3 h: the steps ran at 22–27 ms.
  - **The box**: 512 lattice waters at 33.00 nm⁻³ (0.9872 g/cm³), where Yeh and Hummer measure
    TIP3P's pressure under Ewald at −2.9 ± 2.5 bar, less the five within 1.4 Å of benzene, about
    its volume. W2's 0.997 g/cm³ read +139 to +180 bar, and a hydration free energy carries the
    solute's partial molar volume times the pressure, a few tenths of a kcal/mol there. 24.94 Å, past
    `2 r_c` plus benzene's 5.0 Å. Benzene grown in through the soft core, then 20 ps at 2 fs.
  - **Constrained NVE in that box**, from the start, 1 ps: RMS energy error 0.264, 0.126 and 0.117
    kcal/mol at 2, 1 and 0.5 fs. The floor below 1 fs is the plain cutoff's steps, which no step
    size removes; above it, 2 fs is `h²`'s.
  - **The schedule**: charges off at λ_e 1, 0.5, 0, then van der Waals at λ_v 0.9 to 0 by 0.1;
    29 candidates recorded per sample, one window to be inserted between any pair below 0.1
    overlap. **None was**: the lowest was 0.225, between λ_v 0.4 and 0.3. Each window 20 ps
    discarded and 2000 samples every 100 fs, 2 fs, 1 ps⁻¹, 298.15 K, benzene's QEq charges from
    181L as 3c-3's.

  | window | λ_e | λ_v | TI term | τ | → next: BAR | overlap | EXP fwd | EXP rev |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | 0 | 1 | 1 | +0.775 ± 0.023 | 3.2 | +1.011 ± 0.027 | 0.348 | +1.029 | +0.968 |
  | 1 | 0.5 | 1 | +0.537 ± 0.029 | 1.9 | +0.216 ± 0.017 | 0.395 | +0.208 | +0.247 |
  | 2 | 0 | 1 | +0.378 ± 0.011 | 1.1 | +0.802 ± 0.008 | 0.449 | +0.816 | +0.808 |
  | 3 | 0 | 0.9 | +0.774 ± 0.013 | 0.7 | +0.753 ± 0.009 | 0.443 | +0.333 | +0.750 |
  | 4 | 0 | 0.8 | +0.724 ± 0.013 | 0.9 | +0.685 ± 0.010 | 0.435 | +0.700 | +0.680 |
  | 5 | 0 | 0.7 | +0.644 ± 0.015 | 0.9 | +0.589 ± 0.012 | 0.423 | +0.607 | +0.585 |
  | 6 | 0 | 0.6 | +0.526 ± 0.019 | 1.1 | +0.461 ± 0.017 | 0.403 | +0.425 | +0.473 |
  | 7 | 0 | 0.5 | +0.376 ± 0.029 | 2.0 | +0.240 ± 0.027 | 0.369 | +0.234 | +0.242 |
  | 8 | 0 | 0.4 | +0.055 ± 0.055 | 4.2 | −0.240 ± 0.049 | 0.298 | −0.308 | −0.238 |
  | 9 | 0 | 0.3 | −0.608 ± 0.076 | 5.3 | −1.166 ± 0.065 | 0.225 | −0.859 | −1.234 |
  | 10 | 0 | 0.2 | −1.587 ± 0.083 | 6.5 | −1.326 ± 0.041 | 0.363 | −1.367 | −1.339 |
  | 11 | 0 | 0.1 | −0.911 ± 0.033 | 5.2 | −0.454 ± 0.014 | 0.463 | −0.436 | −0.466 |
  | 12 | 0 | 0 | −0.005 ± 0.006 | 2.5 | | | | |

  kcal/mol; τ in samples.

  | | BAR (kcal/mol) | TI (trapezoid) |
  | --- | --- | --- |
  | charges off | +1.227 ± 0.038 | +1.278 ± 0.038 |
  | van der Waals off, the long-range correction's +0.915 included | +0.343 ± 0.135 | +0.400 ± 0.136 |
  | **decoupling** | **+1.570 ± 0.140** | +1.678 ± 0.142 |
  | **ΔG_hydration = −ΔG_decouple** | **−1.57 ± 0.14** | −1.68 ± 0.14 |
  | GAFF/AM1-BCC in TIP3P, FreeSolv v0.52's calculated value (10.1021/acs.jced.7b00104) | −0.81 ± 0.02 | |
  | OBC II, the same charges (3c-3): polar; with the nonpolar term | −2.34; −1.13 | |
  | experiment, FreeSolv v0.52 `mobley_3053621` | −0.90 ± 0.20 | |

  - **0.67 kcal/mol too negative, 2.7 of the combined errors**, and 0.76 below GAFF in the same
    water. By part: **electrostatic −1.23 ± 0.04**, about half of OBC II's −2.34 for the same
    charges, so GB's over-solvation of QEq benzene was most of 3c-3's polar number; and
    **nonpolar −0.34 ± 0.14**, of which the long-range correction is −0.91 and what lies inside the
    cutoff +0.57.
  - **Hysteresis.** The first half of every window gives +1.377 ± 0.205, the second +1.753 ± 0.184:
    0.38 apart, 1.4σ. EXP forward +1.382 and reverse +1.478 against BAR's +1.570. TI is 0.11
    above BAR, mostly the trapezoid's bias where `⟨∂U/∂λ_v⟩` turns over between λ_v 0.4 and 0.1,
    the windows whose autocorrelation times are longest (4–7 samples).
  - **Window 3's forward EXP, +0.333 where BAR gives +0.753, is one sample.** At step 80 200 of the
    λ_v = 0.9 window, moving to 0.8 lowers the energy by 3.72 kcal/mol, with `∂U/∂λ_v` = +55.5
    against a median of about +8: a water against benzene's repulsive wall, which the softer state
    relieves. Its weight `e^{6.3}` ≈ 540 among 2000 samples is the whole average; without it the
    forward EXP is +0.702. The sample 100 fs before it is the next most extreme (−1.57), so it is
    one event of a few hundred femtoseconds. That is EXP's known failure in the direction that
    removes a repulsion; BAR uses both sides and is not moved, and the interval's overlap is 0.443.
  - **What limits the comparison**, roughly by size:
    - **UFF's van der Waals against TIP3P's oxygen by UFF's geometric rule**, which neither was
      fitted for: a UFF hydrogen's well is 0.044 kcal/mol, aromatic carbon's 0.105 at 3.851 Å, and
      twelve of those against water make the nonpolar part. GAFF's parameters were not read here,
      so how much of the 0.76 against GAFF is this, rather than the charges, is not separated;
    - **QEq's charges**, ±0.098 e, not the AM1-BCC charges GAFF's −0.81 used;
    - **sampling**: ±0.14 by BAR, halves 0.38 apart, one seed, 200 ps a window;
    - **the box**: 507 waters, NVT at the model's one-atmosphere density rather than under a
      barostat; for a neutral solute the finite-size terms are measured as small (above).
  - The steps ran at 22–27 ms rather than the measured 25.7: the machine's load moved during the
    run, and the windows' wall times are in the log.

  **Counts.** `cargo test -p pantometry-forcefield -- --list` counts **378 tests, thirty-four of
  them ignored**: 344 run by default (319 and 30 before). The new: three unit tests in `shake`, one
  in `periodic`, eleven default and one ignored in `constrained_bonds_against_closed_forms.rs`, and
  ten default and three ignored in `benzene_hydrated_in_tip3p.rs`. 319 + 25 is the 344, and
  30 + 4 the 34.

### Changed

- **`Molecule` declares `max force`, `rms force`, `converged` and `minimiser steps` as
  diagnostics.** They describe the minimisation rather than the molecule, and the scene walk's pin
  on the labels every scene emits is what asked: it refuses a new label until somebody decides
  whether a sweep should compare it as an answer. 74 `(label, unit)` pairs now, 63 answers and 11
  diagnostics, and 12 `(domain, label)` pairs; `pantometry verify` prints the four as solver
  diagnostics and leaves them out of the window sweep's verdict.

- **A conjugated ether, ester, acid or thioether oxygen or sulfur is `O_R` or `S_R`, not `O_3` or
  `S_3+2` — `pantometry-forcefield` had been typing it differently from the paper, and the 1d
  measurements showed it.** `uff::assign` now types a divalent O or S with only single bonds,
  bonded to at least one sp² or resonant atom, as resonant. **What was wrong**: with the oxygen of
  anisole and the sulfur of thioanisole typed `O_3` and `S_3+2`, the paper's group-6 sp³–sp²
  torsion row applied and their relaxed barriers were 19.87 and 14.46 kcal/mol against Table II's
  3.6 and 1.7; methyl vinyl ether's C–O–C was 107.84° against Fig 6's 118.3°. The paper fitted
  `O_R`'s radius to methyl vinyl ether's O–CH₃ bond (p. 10025) and its angle to that molecule's
  C–O–C (p. 10028), names anisole and thioanisole "tests of eq 17" (p. 10029) — the sp²–sp² row —
  and types a salicylidene ring oxygen "resonating" (p. 10034). **Now**: anisole 3.6280, thioanisole
  1.6655 — **asserted in Table II, the first of the paper's tests of its method there to agree**
  (the eight rows asserted before are all rows it fitted); methyl vinyl ether's C–O–C 118.05°,
  O–CH₃ 1.4284 Å against 1.428 (asserted) and C–O 1.4137 against 1.413; methyl formate's C–O–C
  115.56° against 113.6° (107.66° before). The paper's for an O on a C=C or a ring; **choices**,
  recorded in `uff::assign`: an O–H (phenol, carboxylic acid), an ester's O, an O between two sp²
  atoms, an O on an sp² nitrogen and a thioether or thioester S are resonant too; an O bonded only
  to sp and sp³ atoms (a cyanate's) is not; neighbours are read before the rule, so it does not
  propagate (a peroxy ester's second O stays `O_3`); the bond keeps the dictionary's order, 1.
  Unchanged: acetaldehyde 0.1728, isoprene 6.6275, ethylbenzene 3.8381 and cis HO–OH 6.5191 —
  none has such an atom; the sign of `r_EN` still wins all three ways on all 24 heteronuclear
  bonds, now with methyl vinyl ether and methyl formate typed by the crate itself. **Aspirin**: its
  acid O1 and ester O3 are `O_R`, and the typing test says so with the reason. At the dictionary's
  ideal geometry its energy falls from 223.095 to 185.672 kcal/mol (torsion 30.673 → 1.334, bond
  2.139 → 5.026, angle 17.007 → 6.035); relaxed, from 18.573 to 29.583 kcal/mol — higher, because
  the ester and acid are now held planar instead of turned perpendicular (C2–O3–C8–O4 −89.0° →
  −17.4°, C3–C7–O1–HO1 90.4° → 179.9°); O4–H1 4.242 → 2.666 Å; heavy-atom RMSD to the crystal
  1.113 → 0.894 Å, now below the ideal start's 0.905. Its minimum takes 454 steps instead of 443,
  FNV-1a `2455b1d7dc555a40` instead of `dac2edd5b0a72aa3`; the facade's ten frames run 185.672 →
  35.084 kcal/mol instead of 223.095 → 72.679. **`Variant::group6_on_sp2` and
  `Group6OnSp2` are removed**: the group-6 sp³–sp² row they switched is now reached only by a
  three-coordinate oxonium oxygen on an sp² atom, which a test builds, and its numbers are held
  through `torsion_parameters` directly. `ElectronegativitySign` stays, as the evidence for the
  sign. Each of the four sabotages of the new rule was caught.

### Fixed

- **The viewer's scale bar did not measure what its label said.** `scale_bar` assumed half the
  screen was half the subject, `(nice / widest) * 0.5` of the half-width, and that ignores the
  focal length `Camera::fit` chooses and the distance the camera stands at. Framed on scene 14's
  20 mm bar, which is 1018 px on screen, the `10 MM` ruler under it was 137 px — 3.7 times short.
  Its test, `the_bar_covers_the_metres_it_names`, checked `across * widest * 2 == nice`, which is
  the same identity, so it could not disagree. The bar now comes from the camera:
  `viewer_core::Camera::across_per_metre` projects the framing's centre and a point one span from
  it along the screen's right direction, through the same `Camera::project` every vertex takes,
  and the label is still rounded up a 1-2-5 ladder to between 0.167 and 0.417 of the half-width.
  It follows a zoom now, which the old one could not. `the_bar_measures_what_the_geometry_measures`
  draws a 30 x 10 x 5 mm box's diagonal (as opened, turned, zoomed in and out) and scene 14's bar,
  and requires the ruler's metres-to-screen to equal the drawn geometry's to 8 `f32` epsilons over
  the shorter length, perspective included in closed form. It agrees to `3.6e-7` at worst. Put the
  old formula back and the box fails at 0.359x. The figures, measured in pixels: `bench-app.png`
  goes from `10 MM` at 136 px to `20 MM` at 165 px. The committed file had 147 px, because it was
  already stale against today's renderer everywhere, not only at the bar. `protein-app.gif` goes
  from `2 NM` at 108 px to `2 NM` at 152 px. Scene 14 framed on its bar goes from `10 MM` at
  140 px to `2 MM` at 128 px. These lengths include the end tick, which is why 137 reads as 140.
  The chooser tiles are byte-identical, since `--thumbnail` draws no legend. `docs/README.md`'s
  bench command read `bench.json` from `app/`, where the example does not write it, and failed. It
  reads `../bench.json` now. **The editor's bar was wrong too, by another route.** It probed along
  world +x through the geometry's projector, which is foreshortened by the view's turn and off the
  centre's depth: at the camera every scene opens on it measured 0.65 of the true scale, a ruler
  drawn at 65% of its label. It takes `across_per_metre` now, and the bracket figure's label moves
  from `20 mm` to `10 mm`; both editor figures are retaken.
- **Scene 14's bar could not be looked at, in the editor or the viewer.** The scene holds a 20 mm
  bar, a 4.4 m room and orbits spanning 2.7e11 m, and both shells framed everything at once, so it
  opened on three planets' markers. Zoom cannot help — `Camera::zoom` clamps the distance to
  1.2..9 — and the editor's **Frame selection** could not either: it kept the world's framing and
  moved only the focal length, on the argument that re-centring "would move everything else
  relative to it". Under that framing the bar is `7.5e-14` of a span, `Framing::local` put both its
  ends on one `f32` (`[0.031563334, 0.06726968, -0.010378995]`), and `Camera::fit`, measuring from
  the framing's centre 2.0e10 m away, moved that one point to `(0.258, 0.850)`, the top edge of the
  viewport, with its ends `2.4e-13` of the viewport apart. The argument was not true of this
  viewport: every reader takes the one framing. **Frame selection now makes the selection's box
  the framing**, settled once per paint above every reader — handles, shaded batches, flat painter,
  labels, scale bar, probe — and the bar's outer end lands at 0.85 of the half-frame with the bar
  spanning 0.76 of the viewport, on `project` and through the `f32` matrix alike. **Fit view**
  returns to the world. Two things the old framing had hidden came with it:
  - **The shaded meshes' cache key did not hold the framing**, though every vertex is built in it,
    so a framing that moved with nothing else drew the old meshes under the new camera. It was
    reachable before this change: drag a domain, then press Fit view.
  - **The CPU paths drew what a GPU clips.** `Camera::project` can only clamp, and framed on the
    bar one planet is behind the eye and projects to `(-2.8e14, -8.8e14)`; box edges were drawn
    from there and labels were clamped back onto the screen from there. `Camera::sees` is the
    depth half of the clip test, and the flat painter, the labels, the probe and the viewer drop
    what fails it.
  The outliner's run rows also carried each panel's own box where `Node::bounds` promises the world
  box, which a re-framing would have aimed at empty space for any placed panel.
- **The viewer can frame one panel**: `--frame-panel NAME`, and `F` in the window to step through
  the whole run and each panel. Framed whole, the bar's 61 samples reached the GPU at one position;
  framed on it, at 61. A `--snapshot` of scene 14 lights **13 318** pixels framed whole and
  **400 616** framed on the bar — most of them the room the bar lies along, which touches it. A
  whole-run snapshot is byte-identical to the previous binary's on the three committed fixtures and
  on scene 14. Every new test was run against the old behaviour and fails on it; the cull, the
  cache key and Fit view's return were each sabotaged separately and each was caught.
- **A molecule's scale read `0 um`.** The editor's scale bar and its inspector's extents had no
  unit below the micron, so three decimals of one rounded aspirin — 0.81 × 0.40 × 0.94 nm — to
  `0.001 x 0 x 0.001 um`, and its scale bar to `0 um`. Below a micron they say nanometres now. The
  glTF and USD writers' note about the radius they draw bodies at had the same shape: `{:.6} m`
  printed aspirin's as `0.000000 m`, a size of nothing, and is in scientific notation now.
- **The viewer coloured every field on the selected panel's scale.** `panel_vertices` worked out
  each panel's own span, as its comment said, and handed the field branch `self.span` instead.
  Scene 25 has a temperature in kelvin and a displacement in metres on one grid, and the
  displacement, some 1e-6 m, was painted on a scale from 300 K: its tile was a solid block of the
  coldest colour, which no test noticed. `a_field_takes_its_own_scale_and_not_the_selected_panels`
  renders one field with each of two panels selected and fails on the old line with every vertex
  the bottom colour. Found by regenerating the chooser's tiles, which adding scene 33 had found
  stale: `06`, `07`, `14`, `25` and `31` had been drawn before the renderer moved, and `18` moved
  with this fix. They are the script's output again. `14-a-world` falls from 170 lit pixels to 14
  and is pinned as the third sparse tile: a 20 mm bar beside orbits 1.5e11 m across, framed
  together, opens on three planets' markers, and the tile now says what the scene shows.

- **Xenon's reason for being out of the catalogue was stated and not measured.** The docs said
  its 10.8% miss on liquid density was an unsourced pair. A sourced one misses too: Beattie,
  Barriault and Brierley's fit to xenon's own virial coefficients, *J. Chem. Phys.* 19, 1222
  (1951), is 3.5% low in triple-point temperature and 8.4% low in density. So it is not only a
  transcription. `a_gas_this_model_cannot_describe` now measures that pair beside the table's,
  and xenon stays out.
- **The protein and the bracket opened from New project as a missing file.** A preset opens as an
  unsaved `scene.json`, which has no directory beside it, so `structures/1CRN.pdb` and
  `parts/l-bracket.stl` were looked for wherever the editor had been started. The protein showed
  an error and a viewport saying it had no geometry. The chooser was supposed to warn about the
  bracket, but its flag came from a search for `"stl"`, so the protein got no warning at all.
  Presets now carry the files they name, the editor reads those first, and saving writes them
  beside the scene without overwriting a file already there. Nothing had opened a preset through
  the editor. `a_preset_opens_with_its_files` now clicks all thirty-two tiles, and it failed on
  exactly these two before the change. It also runs the two, saves one and reverts one: Revert
  after a preset had kept the preset's files attached to whatever `scene.json` it then loaded.
- **An empty viewport said the same sentence whatever the reason.** It said "sources, lumps and
  networks are readings, not places" to a scene that had not checked and to a protein that simply
  had not run yet. `editor_core::nothing_to_draw` now says which of the three it is, and both
  editors call it.
- **Watching a run serialised it quadratically.** `editor_core::run_streaming` sent the whole run
  so far after every frame. For the cavity scene that was about 31 GB to reach a 38.5 MB run, on
  the editor's run thread. It now sends each frame once as `Streamed::Frame` and the settled run
  as `Streamed::Whole`, which the editor appends and then takes. The test that streams every
  scene went from 1565 s to 107 s. `editor-core` is not published, so its new signature changes
  nothing outside the tree.

## [0.22.0] — 2026-09-26

Thirty-six commits, and the theme is making a run say what it holds. A protein was forty-six points
and forty-six numbers; it carries its residues and its bonds now, and is drawn as the surface it is.
A fluid of real argon was written in reduced units under a field that promised metres; it names its
gas and reports SI. The editor gained a library of parts and a panel to drop them into a scene, and
the run file gained a panel shape — which is the one break below, and the reason its format number
moved to 3.

Most of what else is recorded here is what doing that found: a recommended grid that its own build
refused, thirty-nine strings that had lost their line continuations, and thirty sentences of prose
that had stopped being true.

### Changed — breaking

- **`pantometry_scene::PanelData` gained a variant, and `Points` gained two fields.** `Surface` is
  new — vertices, triangles over them and a value per vertex — and `Points` carries `labels` and
  `bonds`. The enum is not `#[non_exhaustive]`, so a `match` over it that names every variant with
  no `_` arm stops compiling, and so does a `Points { .. }` pattern that lists its fields without
  `..`. A consumer builds neither by hand — a scene makes them from its domains — so reading is
  the only way to be affected.

  Nothing else in the published crates moved. Every public enum, struct and trait was compared
  against `v0.21.0`, and all **997** public functions present in both have the signature they
  had. The comparison was first run on a change that did happen — `ligand_from_pdb`'s return type
  between two unreleased commits — so its zero is a result and not a silence.

### Changed

- **A set of bodies with no periodic cell declares the box it occupies, not a cube about the
  origin.** It was `max|coordinate| × 1.2` on every axis: right for an orbit, and wrong for
  anything read out of a crystal file. Crambin, 2.1 × 1.7 × 2.5 nm and about 10 Å off the origin,
  declared a **5.37 nm cube** — seventeen times its volume — and a viewport framing that drew the
  protein at a sixth of the width. Per axis now, with the margin taken from the longest side:
  2.71 × 2.33 × 3.13 nm. Every scene with bodies tightens, and the library gate did not move.
- **A `Fluid` in reduced units reports its speeds as `sigma/tau`.** It said `m/s`, and a wire
  format whose positions are metres wrote 108 atoms spanning 1.7 nm as a box `5.0388 x 5.0388 x
  5.0388 m`; `pantometry view` drew a scale bar reading `2 M` under it, correct about the numbers
  it was given. A fluid built from a real gas still says `m/s` — see `substance` below.

### Added

- **`pantometry_optics::geometry::profile`** — the curve a surface *is*, as points in space, from
  one edge through the vertex to the other. `cap_intersect` and `conic_intersect` say where a ray
  meets a surface; nothing said where the surface was, so `optical_bench` drew its rays and no
  glass at all: light bending in mid-air at angles the picture did not explain.

  The doublet is drawn now, and the claim that makes it a drawing of the thing rather than a second
  description of it is checked point by point: every drawn vertex is where a ray aimed at it lands,
  inside a floor derived from the intersection’s own cancellation. One `Surface` list is both what
  the rays are traced against and what the picture is drawn from.

  The shape is checked against the definition of a conic, `(1+k)z² − 2Rz + h² = 0`, which owes
  nothing to either routine — because the ray check cannot carry it. An axial ray holds `h` fixed,
  so `conic_intersect` inverts the same sag the profile was drawn from and agrees with it whatever
  shape it has: with `k` multiplied by 1.3 inside `conic_sag_si`, every test in the crate passed.
- **`Bodies::label` and `Bodies::bonds`**, both defaulted, so no domain changes by a line. A run
  holding a protein was forty-six positions and forty-six numbers and nothing else, and the editor
  drew it as forty-six unconnected spheres — the picture being honest about the file. `Protein`
  says `chain:number:name` for each body now, and joins two only when they share a chain, are
  numbered one apart and sit within `PEPTIDE_REACH` (4.5 Å, Cα to Cα); `Protein::with_residues`
  carries the names across. The run format writes `labels` and `bonds` only when there are any.
- **`pantometry_molecular::substance`: argon and krypton**, as `σ`, `ε/k` and mass, so a scene that
  names its gas runs in SI. The pair is checked against each gas's measured **triple point**, which
  neither fit saw: argon lands 0.8% out on temperature and 0.4% on density, krypton 1.7% and 2.2%.
  Neon and xenon are kept in `UNSHIPPED` with the measurement that put them there, and the de Boer
  parameter says why they differ: neon's is 0.594 and the potential is wrong for it (8.66% on
  density), xenon's is 0.063 — the most classical of the four — so its 10.77% is a pair that could
  not be sourced. Each shipped residual is also **pinned**, because the 5% band that holds the model
  lets a transposed digit through: `3.405` written as `3.450` is 4.24%. `Fluid::is_reduced`.
- **`Structure::ligand_from_pdb`**, a `HETATM` ligand as atoms with their elements.
  `Structure::from_pdb` refuses `HETATM` on purpose, so one crystal structure could not give both
  halves of the complex it holds. The element is read from columns 77–78 and never guessed from
  the atom name, which is ambiguous exactly where it matters: `CA` is an alpha carbon in one file
  and a calcium ion in the next.
- **`PanelData::Surface`**, a solid a renderer can light. `pantometry view` built triangles only for
  an isosurface through a grid, so a doublet meshed with care arrived as lines and left as lines.
  It crosses the run format, glTF, USD, the filmstrip, the HTML report, the viewer, the editor and
  the GPU shell. The first mesh had **4675** edges not shared by exactly two triangles, for two real
  reasons: every azimuth's own copy of an apex, and the cemented face drawn by both elements.
- **`ligand_binding`**: adenylate kinase open and empty (`4AKE`) and closed on AP5A (`1AKE`). Shown
  only the open structure, the enzyme's softest mode points at the closure with an overlap of
  **0.799**, against 0.039 for a direction in 642 dimensions that knows nothing. Walked along it,
  the distance from the closed form falls 7.13 Å to 4.29 and climbs back — a harmonic direction is
  a straight line through a curved path — and the example writes that as a 48-frame run. The
  backbone is a swept solid of 7680 triangles, held closed and right-side-out against the exact
  area and volume of a prism and an icosahedron rather than against a second sweep.
- **Figures a reader can see the instrument in.** `optical_bench` projects its bench to an SVG
  through the same rotate, tilt and divide as the report's viewer. `docs/bench-app.png`,
  `docs/protein-app.gif` and `docs/editor-protein.png` are the three things only a GPU can draw,
  so each carries a `.txt` of the geometry it is of, and a check compares it on every CI run.
  `pantometry view --all-frames` writes a GIF from one load: 37 s where forty-eight calls to
  `--snapshot` measured sixteen minutes.
- **The viewer draws the run's readings and says which keys do anything.** A run carries the
  numbers it measured frame by frame and `view.rs` did not contain the word, so a window showing a
  protein closing could not say how far it had closed; the controls went to stdout, where a person
  looking at the window is not.
- **A parts library, and an asset panel to put a part in a scene.** `tools/parts/make.py` writes six
  STL solids — plate, rod, wedge, L-bracket, heat sink, pipe — and each is held against the exact
  volume and area of its own profile. The editor lists every `.stl` and `.pdb` beside the open
  scene; drag one onto the view or press `+` and it becomes a domain, a part on the grid
  `pantometry fit` would choose for it. `--ui-dump --drag` drives the drag, with a control released
  over the outliner that must add nothing.
- **`"grid_origin": "parts"` on a block** starts its grid at the parts' lowest corner, so geometry
  modelled about its own centre can be used without moving it. A word rather than a coordinate:
  a millimetre coordinate goes through a decimal and back, and 12.65% of such round trips land on
  a different `f64` — enough, on a grid that fills its part exactly, to refuse it by one ulp.

### Fixed

- **A flat cap rejected a ray aimed at its own rim, and a curved one at the same aperture did
  not.** `cap_intersect`’s flat branch compared `rho2 <= semi²` outright while its spherical branch
  compared `rho2 > semi² + 1e-12`, so an aperture’s boundary depended on how the surface behind it
  curves. Found by drawing: the flint’s flat back face is drawn out to its 9.5 mm rim, and a ray
  fired at that rim met the crown’s front surface 4 mm upstream — the drawn point was on a surface
  that said it was not there.

  The first regression test for it fired from `(9.5, 0, 0)`, where `x² + y²` is `semi²` to the bit
  and the comparison holds either way; reverting the fix left it green. A rim is drawn at
  sixty-four azimuths and `9.5 cos t` squared plus `9.5 sin t` squared is *not* `9.5²` to the bit,
  which is the input the two versions disagree on. The test fires where the drawing draws now.
- **A conic clamped past its own edge returned `R` where the sag there is `R/(1+k)`.** Only a
  sphere and an ellipsoid can reach that line at all — for `k <= -1` the root is never imaginary —
  and for the sphere the two are the same number, which is why nothing noticed. An oblate spheroid
  was clamped `1+k` times too far: at `R = 50 mm, k = 3` it drew the surface at 50 mm instead of
  12.5. The test written for it reached the clamp only by accident and then not at all — a profile
  stops where `inner` is exactly zero, which is the other branch — so it now calls `conic_sag`
  past the edge, where there is no such luck.
- **"Every drawn point" meant one meridian of six, and none of the three rims** — 147 of the
  1383 vertices the picture is made of. The sections are drawn at six azimuths and the rims at
  sixty-four, so an error in the azimuth would have been invisible at `y` and everywhere else in
  the frame; nothing in `pantometry-optics` fired a check off the `y` axis either. Every drawn
  vertex is checked now, and that is what found the flat cap above.
- **A doc said a sphere and its paraboloid are "tens of micrometres" apart at the rim.** Measured,
  174 — a size written where a measurement belonged, and three times out from the liveness floor
  sitting under it.
- **The floor on that point-by-point check dropped a factor of `S/|R|`.** For an axial ray the
  rounding that matters is the one in `R² − h²`, which the square root divides by `2√(R² − h²)`;
  a floor of `4·eps·S` passes only while the standoff and the radius are the same size. Measured:
  at a metre of standoff the gap is `1.430e-15 m` against that floor’s `9.326e-16` and fails. It is
  `4·eps·(S + S²/|R|)` now, per surface rather than one built from the largest radius.
- **`app/viewer-core/README.md` quoted `1046 pixels ... in 3 shades` for the bench snapshot**, with
  "one per field angle" as the check. Lines are depth-shaded now, so the unchanged committed
  fixture reports 77506 in 100 shades — the drift is the renderer, not the run. Nothing asserted
  either number, so nothing said it had moved.
- **The front page counted thirteen domains and named twelve.** `pantometry-protein` was never
  written into the list when the count moved. A correct number beside a short list is invisible to
  a guard on numbers, so the list's own length is held against the crate arithmetic now. Five more
  claims on that page had stopped being true and were measured again.
- **The lens on the front page was three vertical lines**, with the rays drawn straight through
  the glass under a comment saying the bending was too small to draw. `lens_spots` draws the
  surfaces it traces now.
- **A caption claimed an achromat the picture could not show.** Its F and C foci are 21.9 µm apart
  in a frame 118 mm wide, so the three colours coincided and the last one painted was all a reader
  saw. A second panel sets it beside the singlet of the same power, whose colours spread 1546 µm.
- **The citation block restated three numbers nothing compared** — the version, its DOI and the
  concept DOI — and had been wrong four times. Each is held to the thing that owns it; sabotaged
  three ways, including a version DOI where the concept one belongs, and each fires.
- **The legend could not spell its own units.** The glyph table held ten letters and a missing one
  drew nothing, so `refractive index` rendered as `_E__ACT_VE ___E_`; the scale bar knew metres and
  millimetres to three decimals, so a 4.5 nm protein read `0 MM`; and the unit label ran off the
  canvas. The alphabet is complete, an unknown character draws a box, and no two characters may
  draw the same strokes.
- **A scale bar changed length between frames under an unchanged label.** It divided by the
  frame's box while the camera frames the whole run: `2 NM` ran 108 px at the ends of the swing and
  119 in the middle, the ratio of the two frames' boxes. It takes the run's size now.
- **What `pantometry fit` recommended for a part off the origin was refused by the build it was
  for.** It measured every candidate from the parts' lowest corner and left the corner out of the
  fragment, while the builder put the grid at the origin. The fragment says `grid_origin` now.
- **Five of twelve ordinary bricks were refused on the grid `fit` recommended.** The table measured
  the exact `thinnest / across` and the fragment printed it rounded to four decimals, so the grid
  written was never measured, and on the thinnest axis — which has no slack — a cell rounded down
  left the part a sliver outside. The cell is rounded up before it is measured, the counts come
  from the builder's own arithmetic, and all 48 placements tried build and agree with the table to
  the cell.
- **Two of the six parts could not be used by any scene, and the one part that shipped before them
  had no check at all.** `rod` and `pipe` were built about the origin and reached outside a grid
  that started there. Volume, area and the closed-edge count are all translation-invariant; it took
  building a scene around one to see it.
- **The editor's note about hidden panels named the first to go and silenced the status bar.** At
  500 points the assets, the inspector and the scene text had all gone and it said only the first;
  and it was rewritten every paint, so nothing else the editor said survived its frame. It names
  every panel now, when the set changes. `--layout-at` held its own copy of the panel widths and
  would have reported a viewport 190 points wider than the one laid out.
- **The run format gained a panel kind and two keys without its number moving.** `surface` and the
  `labels` and `bonds` of `points` all shipped under format 2, so a format-2 reader — which is
  `deny_unknown_fields` on purpose — met them as an unknown variant and an unknown field: right to
  refuse, wrong about why, and the one confusion the `format` key exists to prevent. The writer and
  the viewer are at **3**, the doc that said the number "has always been one" while it read 2 is a
  table of what each number added, and every key the writer emits is held against that table —
  a key added under an unmoved number now fails with the instruction to bump it.
- **Thirty-nine strings had lost their line continuation**, so their words stood an
  indentation apart: the refusal a person gets for asking for a GPU read `this binary has no
  device.` and then twenty-six spaces, and `pantometry-view` wrote twenty-six into every surface
  panel. Each was a Rust string edited through a Python one, where the backslash and newline are
  Python's continuation and go first; rustfmt does not reflow strings, so nothing said so. They
  came from twenty-nine commits reaching back to 2026-08-08. Each is restored as it was written,
  and a test reads every Rust file in the tree for a string holding a run of nine spaces on a line
  over a hundred columns — the length is what separates a join from the fourteen lines that hold
  such a run on purpose.
- **The wasm CI job asked a script for "the latest wasmtime" and was twice told `{`.** The script
  printed its error and exited 0, so the failure surfaced two steps later as a missing file. It
  installs a named release checked against its hash.

### Found

- **A drawn shape makes a check the prescription cannot.** The bent doublet’s crown is 4.000 mm on
  the axis and **1.801 mm** at its edge; the paraxial arithmetic that chose the curvatures has no
  opinion on whether the two surfaces cross, and a prescription whose elements cross is one no
  workshop can cut. `optical_bench` measures it now.
- **A report draws one card per panel**, so putting the glass in a panel of its own made two
  pictures of one bench — the instrument in one and the light in the other. The glass shares the
  rays’ panel and their scale, and the panel’s unit reads `deg field, 4 = glass` rather than
  leaving a colour bar to imply a lens is a four-degree ray. **FRICTION 49**: `PanelData::Paths`
  colours every path by one quantity, and a drawing of an instrument has structure in it that is
  not a measurement. Recorded rather than actioned — one drawing is not enough to say whether the
  answer is an optional value, a second run list, or a convention.
- **FRICTION’s source table stopped at finding 47** while the sentence above it counted
  forty-eight, and `CLAUDE.md` described the open findings as "the same underlying decision" where
  the file itself gives a different reason for each. The report’s closing section also still said
  eleven domains have scenes, which has been thirteen since `32-a-dose-distributing-and-leaving`.
- **A fan triangulation gets the volume exactly right and the area wrong.** The signed areas of a
  fan telescope to the shoelace sum whatever the polygon, and a fan is still a valid combinatorial
  triangulation, so the solid is closed too. Only the area sees triangles laid across a gap: the
  heat sink fanned is **43%** high on area and exact on volume. `make.py` clips ears, and the test
  builds a U to show all three.
- **Six significant figures is `5e-6` a coordinate, not `5e-7`.** The first tolerance on the parts
  was derived from the wrong end of the mantissa and sat below what the format can produce; the
  pipe, whose cross-section is a difference of two nearly equal areas, amplifies it 2.92 times.
  Moving the parts changed nothing but which coordinates were rounded, and took the pipe from
  `2.0e-7` to `6.8e-6` of a bound advertised as sevenfold. It is `1e-4` and derived.
- **The Lennard-Jones residuals are not an offset a constant could correct.** Argon's triple-point
  density sits below the measurement and krypton's above, and the constant each would need differs
  by 2.58%. The best single value takes argon from 0.39% to 1.27% to take krypton from 2.21% to
  1.31%: a redistribution, recorded beside the residuals so the next reader gets the arithmetic.
- **A window capture loses ninety-six per cent of a solid render's shading.** `PrintWindow` on a
  hardware-accelerated surface kept 247 distinct colours of the 6964 the headless render has, so
  the figures are rendered and not photographed.

## [0.21.0] — 2026-09-15

Fifteen commits, and the theme is the editor: it gained a way in, a way to choose what to simulate,
and one home per command. Most of what is recorded here is what looking at it found — the editor
could not be *looked at* except by opening it, and the hook that fixed that is the reason the rest
of this section exists.

The exception is the twelfth domain, which is not about the editor at all: compartmental
pharmacokinetics, and the count of commits above does not include it.

### Added

- **A structure refines with the block it follows**, so `25-what-140-kelvin-does-to-the-solder`
  stopped skipping its own resolution sweep. A structure's elements have to be that block's cells,
  so refining one without the other is a scene that will not build — which `DomainSpec::refined` is
  right to refuse, seeing one domain. `Scene::refined` sees all of them: it doubles every other
  domain, then every structure that follows one, and refuses when the named partner did not double.

  What the skip cost: the scene's strain energy moves **4.09%** between 512 elements and 4096 —
  0.0274238 J against 0.0263010 — and the strain it is about moves 1.48% in `z`. A first pass at
  that number said 5.31%, from a hand-written refinement that left the block's `contact` faces and
  `dissipation` boxes where they were; the number that matters is the one the sweep makes.

- **A field could be asked for its mean, its peak and its coldest, and not for a place.**
  `Scene::probes` names points, and reports each in every frame beside that domain's own scalars:

  ```json
  "probes": { "junction": { "in": "module", "at_mm": [6.0, 6.0, 12.0] } }
  ```

  **In millimetres, which is the decision the key turns on.** A point in cells moves when the grid
  is refined, so a resolution sweep would compare two places and call the difference
  discretisation. In millimetres the same point reads 44.301778 °C at 1.5 mm and 44.286211 °C at
  0.75 mm — a convergence pair for one design number.

  `24-a-power-module-junction-to-ambient` names its junction and its baseplate now, where it read
  the block's `peak` — the hottest cell anywhere, which is the junction there only because the die
  spans the whole cross-section. The two agree to the bit, which is worth having on its own: `peak`
  comes from the cells the domain holds and a probe from `ScalarField::at`, which interpolates.

  A probe outside the part, on a domain with no field, or misspelled is refused at build with what
  the caller needs; the library needed no change at all.

- **A scene could say what its drive is and not when it changes.** Every scene in this repository
  ran at a constant one, and the questions a design actually asks are not constant: the junction
  temperature of a module under a duty cycle, a motor at start-up, an espresso pulled with a
  pre-infusion. `Scene::stages` is a load profile, keyed by the domain it drives:

  ```json
  "stages": { "losses": [ { "at_s": 300.0, "watts": 12.0 } ] }
  ```

  It reaches the three kinds that have a drive and can be reached today — a `heater`'s `watts`, a
  `conductor`'s `volts`, a `puck`'s `bar`. A stage naming a domain that has no drive, spelling one
  its kind does not take, setting two at once, repeating a time, going backwards, or landing past
  the end is **refused at build** with the domain named; a misspelled drive does not parse.

  **The rule that could have been wrong quietly is step alignment.** A run advances in whole steps,
  so a stage asked for at 10.5 s on a 1 s step lands at 11: the scene heats for half a second
  longer than it says, the conservation audit closes because the joules that were paid were taken,
  and every reading is a correct run of a different experiment. Refused, with a `window_s` that
  lands every stage — **searched rather than derived**, because the first version derived
  `duration / ceil(duration / at_s)` = 10 s and `steps` is raised to `frames`, so the step stayed
  1 s and 10.5 still landed nowhere. A refusal naming a number that does not work is worse than one
  naming none.

  **It hangs off `advance`, not off `run`.** Written the other way first, which would have given
  the batch path one experiment and the editor's streaming path another — the divergence
  `a_streamed_run_reads_back` already exists for, one level up. `World` carries the clock now and
  every path goes through one applier. `World::steps` also stopped keeping its own copy of the step
  count and calls `Scene::steps`, because the check needs that number before a world exists and two
  copies of it would drift into exactly the failure the check is for.

- **`11-motor-thermal-network` starts under load, and shows an overshoot the steady run could not.**
  It ran at a constant 12 W and climbed monotonically to 55.04 °C, because a first-order network
  under a constant source has no other shape. A motor's losses go as `I²R` and its starting current
  is several times its running current, so the winding sees its worst temperature in the first few
  minutes. At 36 W for 300 s and 12 W after it peaks at **72.82 °C** and settles to 59.64 — the
  number the scene used to report is **17.8 K low** for choosing an insulation class. The peak being
  above the end is asserted as a shape rather than a value, because under a constant drive that
  difference is exactly zero.

### Fixed

### Added

- **`Domain::diagnostics`: which of a domain's readings describe the solve rather than the world.**
  A residual, a divergence a projection removes, a norm a unitary scheme preserves, a cell Reynolds
  number — each is worth a column and none is an answer, and a consumer comparing readings between
  two runs has no way to tell. `pantometry-world`'s battery could not, and printed
  `residual 0.000000 -> 0.000000 (281492.026%)` under a heading that says the grid did it.

  It worked around that with a list of five labels it kept itself. This is the same statement made
  where the knowledge is: five crates declare six labels, the battery asks a built world, and the
  answer is a `(domain, label)` pair rather than a label — which removes the ambiguity the list
  needed a pin for, where a *different* domain reporting `norm` as its answer would have been
  swallowed. Additive: `Domain` has thirteen other methods with defaults.

  FRICTION 40 proposed a flag on `Reading` and then argued, correctly, that the fields are public
  so a flag is a breaking change. The argument was sound and about the wrong object.

### Fixed

- **Two findings were fixed and the count did not know.** `FRICTION.md` says how many of its
  findings are resolved, and `friction_counts.rs` pins that sentence to the number of lines opening
  with `**Fixed`. Finding 10 opened with the retraction it was correcting — "this said 'not a defect
  anywhere', and that was wrong. Fixed 2026-09-07" — and finding 36 with what the battery does now,
  above three paragraphs each beginning "the X is fixed". Neither carried the marker, so the file
  said **forty and seven** for a tree that was **forty-two and five**, and the test agreed with it.

  A test that counts a convention cannot see a resolution that did not follow it. What it can see
  is a body that *claims* a fix and carries no marker, which is the only way the two have ever
  disagreed here — measured against the file, that rule flags exactly those two and none of the
  five that are genuinely open. It is a failure now rather than a silent subtraction.

- **The verification battery verified a run nobody wrote.** `run_measured`'s loop advanced
  `world.sim`, the simulation inside the world, which does **one** of the five things
  `World::advance` is: it steps the domains. It does not apply a `stages` profile, hand a structure
  the strain its block's temperature implies, or re-solve it. On a motor with a start-up load the
  battery reported its winding at **44.80 °C** where the run ends at **59.64** — 33% apart, and the
  battery's number entirely self-consistent: with no profile the element ran at 36 W throughout,
  emptied its tank at 800 s and cooled, with every audit margin healthy and no finding raised.

  Its doc has said "the loop is `World::run`'s" throughout, and this is the third time it drifted —
  it computed its own `duration / frames` until `window_s` existed. A loop that *is* the other one
  is the only version of that sentence which cannot go stale, so it calls `World::advance` now and
  `the_battery_measures_the_run_the_world_performs` compares the two reading by reading.

- **The editor ran a different experiment from the CLI, and its own doc said it did not.**
  `editor_core::run_streaming` promised "the last payload is byte-identical to what `run` returns
  for the same text, which the tests pin", and nothing pinned it — the only test on that path
  checked that its JSON *reads back*. It took `duration_s / frames` as its step, which ignores
  `window_s`, the key that exists so the step can be shorter than a frame. On a scene whose window
  is a quarter of its frame the two paths produced **15459 bytes against 15463**, and the streaming
  answer was identical to the same scene with no window at all — the step had never been shortened.

  No shipped scene could show it: all thirty have `steps == frames`, measured, so thirty scenes
  through that path in CI on every commit said nothing. The streaming loop takes `World::run`'s
  schedule step for step now. `the_two_run_paths_are_one_run` compares the bytes on three scenes,
  and the third was added because a sabotage passed: 80 steps over 20 frames is an exact four, so
  swapping `div_ceil` for plain division changed nothing there. A test that pins the step is not a
  test that pins the schedule.

- **A scene named for an espresso shot pulled a quarter of one.** `18-an-espresso-shot` ran for
  eight seconds and delivered 3.28 g from a 4.29 g dose — **4.81% extraction**, where a shot is
  pulled to 18–22%. Nothing failed, and the reason is the shape of every check it had: all of them
  compare the two baskets, and Darcy at a fixed pressure through a fixed bed gives a constant flow,
  so every one holds at any length. It is 25 s now: **19.90%** at 8.31% TDS, 10.26 g from 4.29 g, a
  2.4:1 ratio. The comparison is unchanged — the flow ratio measures 1.786775 either way — and the
  channelling reads harder for it, the gapped basket pouring nearly twice the liquid at 58% of the
  strength. The check asserts the yield lands in 18–22% as a band, because what makes it a shot is
  that it is in the range a barista pulls to.

  The editor's **New project** chooser offered the same quarter-shot, and the gate caught it: a
  `puck` template suggesting 8 s in 9 frames "which no scene that uses it runs". It is 25 s now.
  `the_schedule_a_template_suggests_is_one_a_scene_uses` holds every template's duration against a
  scene that runs it, so a table maintained by hand cannot sit wrong between the change and the
  noticing.

- **A notch's resolution sweep refused to run, and the reason was true of the wrong refinement.**
  `DomainSpec::refined` would not double a conductor with blocked cells, on the grounds that "a
  blocked cell is a one-cell notch, so refining shrinks it — a different geometry, not a finer
  one". That holds only for a refinement that keeps the *indices*: a cell at `(6, 0, 0)` on a 1 mm
  grid occupies `x ∈ [6, 7] mm`, and at 0.5 mm those millimetres are indices 12 and 13. Each
  blocked cell becomes its eight children and the notch is the same notch, which a test now asserts
  in millimetres rather than in indices — a count of eight is not the claim, and the sabotage that
  showed why was one that passed.

  What the refusal cost was the measurement. `17-a-busbar-with-a-notch` is titled "the resistance
  the shape actually has" and every check it had was an inequality or an identity, all of which
  hold at any grid. It is **4.65% above** the `h → 0` answer, and not because the grid is coarse:
  current crowds into a re-entrant corner, the conducting wedge there subtends `3π/2`, so
  `λ = π/ω = 2/3` and a quadratic functional of the field converges as `h^(2λ) = h^(4/3)`.
  Measured 1.33534, 1.33296, 1.33318 and 1.33367 across six grids, bracketing a closed-form
  **1.33333**. Refining is cheap once `R·t` is seen to be exactly thickness-independent — twelve
  digits at `nz` = 1, 2, 3, 5, 8 — but the scene *refuses* 0.0156 mm on its own 1e-9 drift budget,
  which reads 1.094e-9 over 246 k cells, so the finest thing it can say about itself still has
  0.046% in it. The scene keeps its 1 mm cells and gains a check against the Richardson limit, so
  the error it ships with is stated rather than discovered.

  **The reason first written for this was a theorem about a different method.** "A resistance is an
  energy-norm quantity" is the conforming-Galerkin identity, and by Dirichlet's principle it makes a
  voltage-driven resistor read *low* — every grid measured reads high. `Conductor` is a cell-centred
  conductance network with harmonic-mean faces and `Σ g Δφ²`, whose face currents are conservative
  per cell, so the argument is the dual one and Thomson's principle makes every grid an *upper*
  bound. The rate, the exponent, the limit and every percentage survived the check; the sentence
  explaining them predicted the opposite sign from the one every measurement showed.

- **The resolution sweep compared a solver residual as though it converged to something**, and on
  the notch printed `residual 0.000000 -> 0.000000 (281492.026%)` under a heading that says "what
  moved is discretisation". The values are 7.793e-13 and 6.644e-13, both converged and both printed
  as zero by a format chosen for temperatures; the denominator is the **4.08e-17** the residual
  wobbled by between two frames of the base run, which is conjugate gradients stopping at a
  different iterate. It also measured the residual "converging at order **-0.77**" beside four real
  ones.

  Flooring the denominator does not fix it: the wobble is 5.2e-5 of the residual's own magnitude,
  far above any rounding floor, because it is real variation in a quantity that converges to
  nothing. `verify::DIAGNOSTICS` names the five labels the scenes emit that describe the solve
  rather than the world — a residual, a divergence, a `div B`, a norm, a cell Reynolds number — and
  they now print their values with no percentage and stay out of `Sweep::worst`, which the window
  sweep raises a finding on. The scene walk pins all **47** labels the thirty scenes emit against
  the classification, so a domain that gains a residual cannot leave the list quietly stale.
  `Reading` cannot carry this itself: its fields are public, so a flag is a breaking change to a
  published crate. FRICTION.md finding 40.

- **CI's `the app` job ran none of the gate's `--locked`.** Four of its six steps take the flag —
  `fmt` and `deny check` have none — and `CLAUDE.md`'s app block passes it to all four while that
  job passed it to none. Not a policy: the seventh step in the job, added later, has it. The
  divergence points the wrong way: with a lockfile present `--locked`
  *errors* when the manifests would change it and the bare form silently regenerates one, so a stale
  `app/Cargo.lock` stops the local gate and passes CI. That is the tenth failure in `CLAUDE.md`'s
  list with the roles swapped. Its comment also said the suite covers "the twenty-eight scenes";
  there are thirty.

- **The elastic model knew where it stops applying, and nothing checked a scene against it.**
  `pantometry-elastic` is linear, has no plasticity, and documents what that costs: past yield it
  "returns a displacement that is arithmetically correct and physically meaningless, and nothing in
  the answer says which". `Elastic::from_substance` drops the yield strength on the way in, so by
  the time a body is solved the number that would say has gone.

  `25-what-140-kelvin-does-to-the-solder` ships **5.201×** past it — SAC305 assembled at its 217 °C
  reflow and sitting at 40 takes a free strain of 0.3806% against a yield strain of 0.0732%, worst
  in the solder layer at `z = 5`.

  `World` keeps each element's yield strain beside the expansion coefficients it already kept for
  the same reason, tracks the worst `|free strain| / yield strain` over the run, and `verify`
  reports it and raises it above one. **One is a physical boundary and not a chosen threshold**,
  which is what separates it from the arrival measurement above: there is no corpus to calibrate
  against and none is needed, because the model says where it stops.

  The scene is kept — the coupling is what it demonstrates and the coupling is right — and it now
  says what range it is in.

- **A structure that follows a block was passed through the resolution sweep unrefined**, with a
  comment saying it was "reported as unswept". It was not: the block doubled, the structure did
  not, and `World::build` refused the pair because an element and a cell have to be the same box —
  so the sweep reported the whole refined scene as **refused**, and
  `25-what-140-kelvin-does-to-the-solder` had been failing its own resolution sweep since it
  shipped. It skips with a reason now. Refining the two in step is a statement about two domains at
  once, which `DomainSpec::refined` cannot make; left undone rather than done wrongly.

- **A part inside a housing could not lose heat to the air in it.** `cooling` reaches a block's six
  outer faces, so a part rasterised inside one, or a bar with a clearance beside it, had surfaces
  no film could reach — it shed heat by radiating across the gap and by nothing else, and a scene
  describing a housing was describing an evacuated one.

  It failed by producing nothing: state a film on a face made of void and `cells_on` counts no
  solid cell there, so the film is charged to nobody. A copper bar walled in that way sat at its
  initial 200 °C for a whole run, and the run completed and reported four figures.

  `23-a-part-radiating-to-its-lid` is that scene. Its part settles at **536.22 K** as shipped and
  **471.22 K** with still air in the cavity, over the same 600 s — **65 K**, on a scene whose title
  says housing. Its closed form is the same pair of ODEs with one term added to each body, and it
  agrees to **0.12%** where it agreed to 0.20%.

  Air is not a replacement for a stated area, and the busbars are the counter-example:
  `30-two-phases-crossing-at-a-clearance` was tried this way and came out 33.0 K against 22.0,
  because attaching a bar's side area to a block face carries convection *and* radiation over it
  while air carries only convection. Air models a closed cavity; a stated area models a surface
  open to the room. FRICTION.md finding 37.

- **Two busbars four millimetres apart could not see each other.**
  `30-two-phases-crossing-at-a-clearance` states two blackened bars, ε = 0.9, crossing at a
  clearance. Two surfaces at 335.6 K and 319.0 K facing across that gap exchange **6.5711 mW** over
  the 8 × 8 mm patch where they cross — 7.64% of what the cooler bar dissipates — and two `Solid3D`
  domains exchange nothing but bus totals, which carry an amount and no location. The arrangement
  the scene is named for was not modelled.

  As one block it is `find_gaps`, the pairing `23-a-part-radiating-to-its-lid` is already checked
  on. Its closed form is a **pair** now — each bar sheds convectively and radiatively *and* trades
  with the other — giving 22.0246 and 6.2523 K against a measured 22.0218 and 6.2493. Uncoupled
  they sit at 22.4333 and 5.8139, so the exchange is forty times the tolerance the balance is
  asserted to. The `poses` entry it was written for is on the assembly.

- **A `void` region could not be filled back in, and the comment beside it said otherwise.**
  `Solid3D::empty` was one-way: `fill` set a cell's material and left it empty. The scene format's
  own words are *"applied in the same order as any other region — later wins where they overlap —
  so a gap cut into a part reads the way a coating on a layer does"*, and that is the seventh
  comment in this workspace to guarantee something the code did not do.

  It surfaced by writing a scene the readable way: two crossing bars are most simply stated as
  "void the block, then put the bars back", and that was refused for dissipating heat in cells it
  had just filled. The loud version. The quiet one is a part with a hole nobody put there.

- **The flat-field measurement could be silenced by merging two domains.** A domain's `peak` and
  `coldest` are the extremes of everything in it, so one block holding two objects at different
  temperatures reads as a large spread. Giving the busbars their radiative path would have made the
  finding stop naming them with neither bar gaining a field.

  It is **per connected body** now, six-connected over the cells a field panel reports as finite.
  That found a scene it had been blind to: `23-a-part-radiating-to-its-lid` read 0.79 as one domain
  and passed, and its part and lid are two isothermal objects at 0.052% and 0.049%.

- **A bracket rasterised from a shape had no route along it.** `29-a-designed-bracket-becomes-cells`
  is the one scene whose geometry comes from a file, and the route heat takes along that geometry
  is the only reason the shape is in it. It was a uniformly hot bracket shedding over its whole
  footprint, so every cell shed where it stood: 4 100 cells, all holding the same number, `Bi =
  8.6e-4` and a field spanning **0.067 K**.

  It carries a module's 20 W from the tip of one arm into a bolt pad at the tip of the other now —
  **52.4420 °C at the module against 30.5449 at the bolts**, 21.90 K across a 32.44 K rise. Its
  `conservation_tolerance` moves 1e-9 → 1e-8: it has 4 000 J passing through it rather than none,
  and the refined sweep needs more of the budget than the base does.

  It starts at **44 °C** rather than 20. From cold it took 900 s to settle, nearly all of it
  charging the bracket’s mass rather than establishing the gradient the scene is about — and nine
  test binaries walk every shipped scene, so that cost was multiplied by nine. 200 s from near the
  answer is within 0.09% of the settled peak and seven times cheaper.

- **Three true claims about eight hundredths of a kelvin.** `19-a-coating-stops-the-heat` asserted
  that the metal levels, that the largest step along z lands on the interface, and that the glass
  stopped it — all true, and all about a 0.08 K range, because the pulse was one cell at +60 K
  spread over 1 458 cells. It is the whole heated face now: the interface step is **9.4656 K of a
  19.1152 K rise**, half the profile on one face.

  As a `regions` entry rather than a `hot_spot` it can also be refined — `verify` refuses to halve
  a one-cell feature and says to state the initial condition as a region instead — so the scene has
  a resolution measurement for the first time, at 0.106% on the peak.

- **The resolution sweep did not double a cooled patch.** A patch is stated in cells, so at twice
  the grid the bracket's bolt pad came out a quarter of its own area, and the sweep reported
  **5.593%** on the peak: a changed boundary condition reading as discretisation. It is 1.409% now,
  which is the STL rasterising differently at 1 mm.

- **`verify` could not say that a scene's grid was doing nothing.** Every finding the battery
  raised was about arithmetic — determinism, a sweep, a drift, a rasterisation loss — and none of
  them asked whether the scene needed the arithmetic. A block whose cells all hold the same number
  passes every check it has, converges perfectly, conserves to twelve digits, and answers a
  question one ordinary differential equation answers; the report reads exactly like the report of
  a scene with a real gradient.

  It measures the last frame's `peak − coldest` over the range that domain's readings covered
  across the run, against a threshold with its corpus written beside it. Of thirty shipped scenes,
  eleven domains report both, and the gap is a factor of five:

      0.00000  30-two-phases-crossing         phase_a and phase_b
      0.00134  19-a-coating-stops-the-heat    joint
      0.02108  29-a-designed-bracket          bracket
      ----------------------------------------------- 0.05
      0.10645  24-a-power-module              module
      0.79150  23-a-part-radiating-to-its-lid housing

  `30-two-phases-crossing-at-a-clearance` states two busbars of thirty-two cells apiece whose
  spread is **exactly zero**, and `29-a-designed-bracket-becomes-cells` rasterises 4 100 cells from
  an STL to hold 0.067 K at `Bi = 8.6e-4`. `pantometry verify` exits 1 on six shipped scenes now,
  and `scene.rs` pins that set so a seventh is a failure rather than a quieter list. FRICTION.md
  finding 36.

- **A closed form agreed to `1.1e-4` with a power module whose junction temperature was 32.5%
  wrong.** `24-a-power-module-junction-to-ambient` is checked against a resistance stack written
  out term by term from the geometry, and it had agreed to `1.1e-4` of a 141 K rise since the day
  it was written. Both the check and the solver were faithful to the stack in the file; the stack
  was not the module.

  Its 100 µm solder was a 1.5 mm region, fifteen times its own resistance, because 1.5 mm is one
  cell. Its 0.63 mm DBC ceramic was another. And the **largest** resistance in a junction-to-ambient
  path — the interface material under the baseplate — was absent altogether, because a mounting is
  on the one face with no cell on the other side of it.

      as checked, to 1.1e-4    3.1379 K/W    junction 181.19 °C
      with the joints stated   4.1594 K/W    junction 227.15 °C

  The two errors pushed opposite ways and left something plausible. A closed form computed from the
  file cannot see the file: a premise shared by a model and its check is not tested by their
  agreement, however tight — and tightness is what makes it look tested. Recorded as EVIDENCE.md's
  twelfth section and FRICTION.md's finding 35.

- **A prose-count test hard-coded the total it was checking, twice more.** `counts_in_prose.rs`
  already carried two comments recording that mistake — a hard-coded half refuses a correct
  document the day a finding is actioned — and two of its templates still spelled "thirty-four" by
  hand. Finding 35 turned the gate red on eight documents that were all correct. Both halves come
  from the count now, everywhere.

- **A designed part received 61% of the cooling its scene stated.** `losing_from` takes an area and
  the block divides it among the cells on that face so each carries its share — divided by the
  cells the *grid* has there, and a part rasterised from an STL is mostly not there. A void cell
  took a share and lost nothing with it.

  Measured on `29-a-designed-bracket-becomes-cells`: 410 of the 676 cells on `z-min` are solid, so
  16.5 cm² acted as **10.0**, and its 120 s drop was 1.908 K against the 2.950 K its own lumped
  balance predicts. It is 3.151 K now, a ratio of **1.068** to that same balance — and the 6.8%
  over is the radiation the lumped form omits: 6061's ε = 0.09 gives a secant `h_rad` of
  0.842 W/m²K against 12 for the film, which predicts 1.070.

  `cells_on_where`, one line above `cells_on` in the same file, has filtered voids since it was
  written.

- **Asking for more pictures changed the answer by 58%.** A scene's `frames` set the coupling
  window as well as the capture rate, and each domain subdivided that window into whole substeps
  no longer than its own stability limit — so the step, and with it the accuracy, was a function of
  the frame count. `15-a-hot-spot-in-a-block` reads 5.177, 7.027, 8.065 and 8.186 K at 2, 11, 101
  and 1601 frames, textbook first order, and shipped at 11: **14% below its own grid's converged
  answer**. A stability limit is not an accuracy limit, and nothing in the format said what would
  be right.

  `window_s` separates them. With one stated the run takes `ceil(duration_s / window_s)` whole
  steps and `frames` chooses only which are photographed, asserted bit-for-bit. Absent, the old
  behaviour stands, so every scene written before the key is unchanged.

- **`verify` measured this all along and reported it as 0.212%.** The window sweep divided the
  shift by the reading's own magnitude, and for a celsius temperature that is dominated by the
  273.15 the scale carries. Against what the run actually did the same shift is 0.966%. The
  denominator is now the range the readings of one unit covered in one domain over the run —
  per domain and unit rather than per reading, because a block's `peak`, `mean` and `coldest` are
  three views of one field and judging each against its own travel made `mean` read 9.2% beside two
  printed values that were identical.

  And it is a finding now rather than a row, above half a percent, carrying the exit code. Three
  scenes tripped it and were corrected: `15` to 88 frames, `05` to 52, `27` to 1600.

- **The battery marched a different scene from the one the CLI runs.** `run_measured` computed its
  own `duration / frames`, which was the same number until `window_s` existed. Its own doc said the
  loop was `World::run`'s. Both step by `World::steps` now, and the window sweep halves the window
  rather than multiplying the frames — which on a scene that states one had become a knob that
  moved nothing.

- **Two scenes called themselves a crystal and a liquid and were the same frozen lattice.**
  `08-atoms-crystal` and `09-atoms-liquid` asked for `duration_s: 6.0e-12` against a domain built
  with `LennardJones::reduced()`, where the only time scale is `τ = σ*sqrt(m/ε)` and all three
  are 1 in SI — **one second**. The whole run was 6e-12 of it, ten orders short of the 0.01 the
  domain suggests for a single step. Six picoseconds is where *argon* melts, and argon's is
  2.16 ps.

  Measured: mean square displacement grew as exactly `t^2` — a ratio of **4.000** between the last
  frame and the half-way one, which is free flight with not one collision — and the "liquid" sat
  on the crystal's own fcc sites to 1e-11. At `duration_s: 6.0` the crystal saturates at 0.075 of
  a neighbour spacing (Lindemann melts near 0.1) and the liquid reaches 1.174 with a ratio of
  2.42, which is diffusion.

  **Every check those scenes had went on passing**, because the speeds were right and only the
  clock was wrong. And the 5e-2 conservation tolerance passed with 11.9 digits in hand on a run in
  which nothing moved; at the corrected duration it has 1.5. It had never been tested. `scene.rs`
  now asserts the displacement, which is what equipartition cannot see.

- **The chooser's timescale span said fourteen orders of magnitude and it is sixteen.** The
  extremes are a quantum `well` at 2e-13 s and an `orbit` at 7200 — 3.6e16 — and they always were,
  so the count was stale whatever `atoms` held. The example pair named in three doc comments was
  `atoms` against a thermal `network`, whose "3e14 times apart" was a product of the units error
  above; corrected, those two are 300 apart and the warning does not fire for them.

- **Every icon in the outliner was a missing glyph.** Five codepoints, all absent from the bundled
  face, so every row drew two hollow squares — a fold arrow and an eye — for as long as the
  outliner has had rows. Painted now.

  `--ui-dump` could not see this and no test could: it reports the *string*, and the string is
  fine. `unwritable=N` asks `epaint`'s own `has_glyphs` of the family each string was laid out in,
  which is the question the renderer will ask, and marks the line with `?`. Zero across seven
  screens at three widths.

- **Going back from the kinds screen to the pictures dropped every tick.** The set is taken out of
  `choosing` with `mem::take` at the top of the frame and that arm neither restored nor consumed
  it — harmless while `Back` meant "leave the chooser", wrong once it meant "up one level".
  Seeing it needs two clicks and `--ui-dump --click` took one point, which is why no state change
  across screens had ever been checked; it takes several now, in order.

- **`thermal::network`'s `Node` documented a guarantee that does not exist.** It said
  `Simulation` "already refuses two domains with one name", and `Simulation::with` pushes onto a
  `Vec` with no check: two domains sharing a name advance without complaint and `domain` returns
  the first. The handle's identity is a hash of that name, so the doc was resting a safety
  argument on it. The real limit is written down instead.

### Added

- **A block's steady state, solved rather than marched to.** `Solid3D::steady_state` finds the
  field at which every cell balances and `settle` puts a block at it. The residual is `flux_at` —
  the function the sweep itself marches with — plus the source and the clearance pairs the sweep
  applies beside it, so it finds the march's own fixed point rather than a second opinion about
  where it is.

  Successive over-relaxation with Newton on the diagonal, matrix-free.
  `ThermalNetwork::steady_state` builds a dense Jacobian and factors it, which is right for a
  handful of nodes and impossible for eight thousand cells.

  Checked against the closed forms and against the march: a driven bar at `293.579509 K` solved
  against `293.579508` from the series and `293.579507` marched; a radiating lump at `588.520214`
  against the root of its own quartic to six decimals; and a pair across a clearance at
  `388.5512/310.7172` against a march that needed **32 686 s** of simulated time to get there.

  On the shipped bracket the solve moves the marched answer `0.0465 K` further, in **1.15 s**
  against the 470 s that march cost each of nine test binaries.

- **`verify` reports how far short of its own steady state a run stopped.** "Did it arrive" was a
  judgement — look at the last two frames and decide — and it is a number now.

  It is a **measurement rather than a finding**, and that took measuring. Written as a finding at
  a percent, it named `19-a-coating-stops-the-heat`, `15-a-hot-spot-in-a-block` and
  `23-a-part-radiating-to-its-lid`, all three of which are transient on purpose: `19`'s claim is
  the interface step before the front reaches the glass, and running it to steady state moves the
  steepest step off the interface and breaks it. What separates a deliberate transient from a run
  cut short is the question the scene is asking, and this format does not carry one.

- **A block's void can be air.** `Solid3D::air_in` fills it at a stated temperature with a film,
  and every solid face touching a void cell sheds `h · dx² · (T − T∞)` to it. A scene spells it
  `air`. The area is the **grid's** rather than the caller's: a `cooling` entry states what an
  outer face exposes because a rasterised part covers less of it than the grid does, and the faces
  touching an internal void are exactly the ones the grid has.

  **Convection only.** A transparent gas does not radiate; what a surface facing a clearance
  exchanges with the surface across it is the pairing `find_gaps` already computes, and the two run
  beside each other because both paths are real and in parallel.

  Checked against the lumped exponential to the step's own first order, with the ledger agreeing to
  the bit — `201.105948 J` shed and `201.105948` counted — and against two conductances in parallel
  at steady state, `36.07013 K` against `36.07013`. Absent is a vacuum, so no scene changes unless
  it asks.

- **A face can be cooled only where it is bolted.** `Solid3D::losing_from_within` exposes a box on
  a face rather than the whole of it, and a scene spells it as `cooling`'s `from`/`to` — the
  face's two axes in x-y-z order, the reading `contact` already gives its own.

  A part is bolted at **pads**, and a face could only be cooled entire. Stating a smaller
  `area_cm2` does not stand in for it: the area is divided among the cells on the face, so a tenth
  of the area is a tenth of the conductance spread over all of it. Measured one step from uniform,
  a whole face and a 3×3 pad of the same stated area lost `0.043328955 J` each — the same joules to
  the bit — and left a block with `0.000000 K` corner-to-centre against one with `0.248756 K`.

  Checked against the two closed forms a bar has: bolted at one end, `345.9747 K` against
  `P·((n−½)·dx/(kA) + 1/(hA))` of `345.9747`; the same conductance spread along its length,
  `263.3732 K` against the fin equation's `263.2563`, agreeing to 0.044%. That second one corrects
  the claim it was written to make — spreading the conductance does not give a lump, it gives a
  fin, and the lumped balance is 222.22 K.

- **A face can carry a stated contact resistance.** `Solid3D::joined` puts a conductance on an
  interior face, in series with the two half cells already there — `1/k_face = 1/k_series +
  1/(dx·h)` — and `Solid3D::mounted_on` puts one on an outer face, between the half cell and the
  film. A scene spells them `contact` and `cooling`'s `contact_w_per_m2_k`.

  Every real thermal design is a chain of these: a die soldered to a substrate, a substrate bonded
  to a baseplate, a baseplate bolted to a heatsink through grease. All are tens of microns thick in
  a part discretised at a millimetre, and none of them could be stated. That is a **floor**, not a
  discretisation error: a `regions` entry is a box of cells, so the thinnest resistance it can
  express is `dx/k`, and refining the mesh lowers the floor rather than approaching an answer. No
  grid a real part can afford reaches 100 µm — on a 12 mm module that is 120 cells a side.

  A contact has a resistance and no thickness. Checked against `1/(A·h)` exactly across five
  decades of `h`, and against a marched steady state to four figures: 226.7028 K measured against a
  closed form of 226.7028, with the joint itself 78.1250 K against `P/(hA)` of 78.1250. `h = +∞`
  returns the harmonic mean bit for bit; `h = 0` is a clearance; a negative one is refused by name.

  It only ever lowers a conductance, so the explicit stability limit gets looser and a joint can
  never make a stable step unstable. The same conductance is read from both sides of a face, so
  conservation is untouched.

- **`pantometry-pharmacokinetic`, the twelfth domain: where a drug goes in a body.**
  `CompartmentModel` is *n* well-stirred compartments joined by intercompartmental clearances,
  with elimination out of the ones that eliminate — the same algebra as `ThermalNetwork`, whose
  module docs it copies the reasoning of, plus the two things that are not in a thermal network:
  drug leaving the system, and a dose arriving from outside it.

  **The kernel is unchanged, which is the twelfth time that claim has held.** The domain names no
  other domain and never touches the `Exchange` — it neither produces nor consumes anything
  another domain could want, and publishing on a channel nobody consumes is refused by the bus's
  own audit, correctly.

  Checked against the closed forms and nothing else: `A(t) = A₀e^{−kt}` for a one-compartment
  bolus; `(R/k)(1 − e^{−kt})` for an infusion, whose steady state `R·V/CL` explicit Euler
  reproduces *exactly* because the fixed point of the discrete map is the fixed point of the
  differential equation; the bi-exponential for two compartments, with `α` and `β` the roots of
  `λ² − (k10 + k12 + k21)λ + k10k21 = 0`, which is the only one of these that would notice a
  transposed index in a link — a link contributes `+q` and `−q` to the same sum, so the
  conservation audit is blind to it by construction; mass balance to `1e-12` of the dose over
  14 400 steps; and the **rate**, four step sizes each half the last, with the error required to
  halve within 3%.

  `max_stable_dt` is derived rather than chosen. Gershgorin on the concentration-basis matrix
  puts every eigenvalue within `Σ Q_ij/V_i` of `−(ΣQ + CL)/V_i`, so `|λ|max ≤ 2 max_i r_i` and
  `h ≤ 1/max_i r_i` implies the `2/|λ|max` explicit Euler needs. The same expression is the
  **positivity** limit — `A ← A(1 − h r) + …` — and that is the half worth having, because the
  scheme is still *stable* at twice this and a model run in between holds negative amounts of
  drug without ever diverging. Verified against the true `α` of the two-compartment model, which
  is available in closed form, rather than trusted.

  Five deliberate absences, in the module docs: no absorption compartment (`kₐ` acts on an amount,
  a clearance acts on a concentration — different mechanisms), no Michaelis–Menten elimination
  (which is exactly where the closed forms stop), no pharmacodynamics (a response model has no
  conserved quantity), no population variability (a sampling problem on a different execution
  axis from the one this kernel has), and no AUC as a typed quantity.

- **`VolumetricFlow` in `pantometry-units`, m³·s⁻¹ — `MassFlow`'s missing sibling.** New public
  API in a published crate, so it is here rather than only in the domain that wanted it. A
  pharmacokinetic clearance is a volume of plasma emptied per unit time and there was no name for
  that dimension; `pantometry-pharmacokinetic::Clearance` is an alias of this one. It is in the
  units crate rather than the domain crate for a reason the orphan rule decides: an inherent
  `impl` on a foreign type is not allowed, so a local alias could carry no `l_per_h` or
  `ml_per_min` — and a clearance whose only constructor is `from_si` is a factor of 3.6 million
  waiting to happen. `product!(Concentration, VolumetricFlow => MassFlow)` and
  `product!(VolumetricFlow, Time => Volume)` come with it, so `CL·C` landing in kg/s is checked
  by the compiler rather than claimed by a comment.

- **`--ui-dump`, one frame of the editor's interface as text.** `App::update` split into an eframe
  entry point and `App::ui(&egui::Context)` — the `eframe::Frame` argument was already unused — so a
  frame can be built, laid out and read with **no window and no GPU**: the shaded viewport arrives
  as a `Shape::Callback` that a headless run records and never executes. Two `ctx.run` passes,
  because egui settles.

  It reports every string with its rect, the paint callback's own rect, `callbacks`, and four
  counts over the strings: `texts`, `elided`, `cut` and `overlap`. Doors for the states a frame built from nothing cannot reach —
  `--click x,y` sends a press and a release, which is the smallest input that opens a menu;
  `--ran` runs the scene through `editor_core::run` so everything that only exists after a run
  does; `--iso`, `--solo` and `--new` set the three pieces of state a click cannot reach from
  outside.

  Twenty tests over it. Five of the defects below were found by looking at what it printed — the
  zero-width viewport, `fit view` laid out at x=704 of a 700-point window, `Domain` clipped at 160,
  the status bar that stopped wrapping, and the time label cut at 260 — and **three of its own
  columns were measuring nothing when written**: `elided` counted the ellipsis that means "opens a
  dialog", `cut` watched one edge after the strings that ran off the other had moved, and `overlap`
  compared clip rects, which is what told the colliding pair apart.

- **A start screen.** `pantometry edit` with no file opened the built-in scene, which is an answer
  to a question nobody asked. **New project**, **Open a scene** and what you had open before. The
  recent list is a JSON array beside the platform's other per-user configuration, written with
  `serde_json` rather than by turning on eframe's `persistence`, which would add `ron` and `home`
  to carry eight strings. A path that has gone is shown greyed and refused rather
  than dropped: a list that silently shortens itself looks like a list that forgot.

- **New project asks what you are simulating.** Every kind the scene format defines, with a line
  about each. Ticking one is a starting point; ticking several is a scene of several domains, which
  is what a custom simulation is here — there is no second flow for that case, because it is the
  same list.

  **It says what a checker cannot.** A scene carries one `duration_s` and these kinds do not share
  a timescale: `atoms` settles in 6e-12 s and a thermal `network` in 1800 s, fourteen orders of
  magnitude between them, and sixteen across the whole table — `well` at 2e-13 s to `orbit` at
  7200 s. Every combination of them is a *well-formed* scene, so nothing
  downstream refuses it — what happens is that the slow domain advances by a ten-trillionth of what
  it needs while the fast one finishes, and the picture looks like a bug in the physics. The
  chooser computes the span from the same numbers the scene is built with and puts it on the
  screen.

  The `duration_s` and `frames` a new scene starts from are not chosen: each comes from a shipped
  scene that uses that kind, and a test holds all nineteen against the thirty on disk. A starting
  point whose duration nothing has ever run is a starting point that opens on a refusal.

- **`File > Open…` and `File > Save as…`**, through the platform's own dialog. `rfd` is two new
  crates — everything it needs here was already in the tree for the window — and on Linux it talks
  to the desktop portal rather than linking GTK.

- **An eighth reviewer, `unearned-pass-hunter`**, for the failure the other seven are downstream
  of: did the check *run*, and against the thing it claims to test. Its table of disguises is drawn
  from instances this repository has shipped.

- **`EVIDENCE.md` and `EXAMPLES.md`**, and figures under `docs/`.

### Changed

- **`README.md` is 169 lines and was 883.** Eighty of them said what this is; the rest was
  explanation. What it is, install, run, and a table of where everything else went: ten sections
  of verification into `EVIDENCE.md`, the crate table and *What is not here* into
  `ARCHITECTURE.md`, why those CI jobs exist into `CONTRIBUTING.md`, and the API sample and the
  example table into `EXAMPLES.md`.

- **Seven of the toolbar's nine controls were also in a menu**, and the screen did not say which
  was the real one. `open`, `revert` and `save` are in File; `deep` is in Run; `watch file` and
  `run on change` are in Watch; `fit view` is in View. What is left is **run**, **verify**,
  **deep** and **fit view**, on a rule worth writing down: *state that changes what a button does
  stays beside that button; state that changes what happens when nobody presses anything goes to
  the status bar.* The strip is one row at every width down to 420 points; it was two at 500 and
  three at 420.

  The path field is a label. It was 220 points of editable absolute path, and it was a control
  because typing a path was the only way into a file on disk.

- **Seven menus, each named for one subject.** `View` held eleven items doing three jobs — how the
  viewport draws, which panels are on screen, and where the camera goes — and was the only menu
  whose name did not predict its contents. The three panel toggles are in a `Window` menu now. A
  seventh name to read is the price; a menu you can guess the contents of is what it buys.

- **The isosurface level left the View menu.** A menu is where a thing is turned on and where a
  command is issued; it is not where a *continuous* value is set, and this one was inside a menu
  drawn over the viewport whose surface it moves. The toggle stayed; the slider and the note that
  goes with it — what this frame reaches, and whether the level is inside it — are on the strip
  along the bottom of the viewport, under the colour bar whose scale they share.

- **`solo` is offered once.** It acts on the selection, the outliner is where a selection is made,
  and its header carries the checkbox with the state visible — so the View menu's copy went. Hiding
  the outliner would then have left the viewport drawing one domain out of five with nothing saying
  why, so the status bar says `solo — drawing only the selection` exactly when the control is not on
  screen to say it itself.

- **`File > Load` is `File > Revert`.** It re-reads the path already open, which is a revert; called
  `Load` it was the only item on that menu that looked like a way to open a different file.

- **The gate is 29 steps.** The two wasm targets were added after a test that reads documents off a
  disk shipped without `#![cfg(not(target_family = "wasm"))]`, passed twenty-seven local checks and
  turned CI's `test (wasm32-wasip1, wasmtime)` job red.

### Fixed

- **The editor's viewport was zero points wide at four window sizes, under a guard that passed.**
  `--layout-at 900` printed `view=290`; the rect the paint callback was handed was `0 x 870`.
  `panels_that_fit` reserved each side panel's *minimum* — 170, 240, 200 — and the panels were then
  created at their *defaults* — 260, 430, 320 — so the sum it checked against the 280-point floor
  was 400 points short of what the layout spent. Zero at 1000, 950, 900 and 700 points.

  That is the same defect `the_viewport_always_has_room` was written about, still there, because
  that guard reads the decision function and the decision function was not the layout. It passes at
  177 widths — 200 to 3 192 in steps of 17 — and **not one of the four is among them**: 897 and
  914 straddle 900. The panels are held to the budget the arithmetic
  reserves now, and both sums count the separator egui puts beside a panel — with three panels up
  the viewport settled at 264 against a floor of 280, and 16 points is two of them.

- **A parse error hid every report the editor made about its own actions.** The status bar was
  `match error { Some => error, None => status }`, and thirty-three places assign that status:
  *loaded X*, *saved X*, *ran: 11 frames*, *still busy with the last job*. Saving a scene with a
  missing field wrote the file and said nothing. Both now, left to right in reading order.

- **The editor could not read back the JSON it had written one second earlier.** JSON has no NaN. A
  `block` filled from an STL marks every cell outside the part not-a-number and
  `pantometry_view::data::compact` writes any non-finite value as `null`, which is the only thing
  JSON offers; `viewer_core::Run` declared those values `Vec<f64>` and refused the file. **Two of
  the thirty scenes** — 2 772 nulls in `23-a-part-radiating-to-its-lid` and 56 168 in
  `29-a-designed-bracket-becomes-cells` — so neither `pantometry view` nor the editor could open
  what `pantometry run` had just written. `null` reads back as `f64::NAN` now, which is the value
  the writer had.

  CI ran all thirty scenes through the binary and read none of them back through `viewer_core`,
  which is the reader the viewer and the editor both use. Thirty of thirty were green while two
  could not be opened. `a_streamed_run_reads_back` closes that.

- **Four things took the viewport's bottom and none of them knew about the others.** The scale bar,
  the colour bar, the count of what the shaded pass drew and the frame transport each anchored to
  `rect.bottom()`, and with a run open the slider was laid out through the middle of the count. A
  band each — 26 points — with the transport's rows reserved above them. The time label moved above
  the frame slider rather than beside it, because side by side they are one row that does not fit a
  260-point viewport.

- **A stale count in `CLAUDE.md`**: "the twenty-eight scenes" when there are thirty. The scene count
  was guarded in **seven** places, five of them sentences in one README, which is how a number can
  be guarded seven times and wrong in an eighth — and in a ninth, `CONTRIBUTING.md`, which this
  changelog's own audit found still saying twenty-eight. Ten places now. The menu bar and the
  toolbar both wrap as well, so neither runs off a narrow window.

- **`README.md`'s BibTeX and the version DOI.** 0.20.0's Zenodo record is
  [10.5281/zenodo.22233493](https://doi.org/10.5281/zenodo.22233493), and the citation block was a
  ninth place a version lives that no list covered.

### Added

- **`pantometry-protein`, the thirteenth domain: a coarse-grained elastic network model.** One node
  per residue at its alpha carbon, one spring between every pair inside a cutoff, and the
  eigenvectors of the Hessian are the collective motions the fold has. `Structure::from_pdb` reads
  a Protein Data Bank entry, `Network` builds the springs, `Modes` diagonalises, and `Protein` moves
  the structure along those modes at a temperature — written down in closed form rather than
  integrated, so there is no integrator error and the energy is exactly `(3n − 6) k_BT`.

  **No new dependency.** The workspace still resolves twelve external crates: a dense symmetric
  eigensolver is two hundred lines of cyclic Jacobi with a closed form for every part of it, and it
  is checked against `4 sin²(jπ/2N)` for a free chain, the Dirichlet chain's own eigenvectors, a
  planted spectrum with a triple root and a zero, `Σλ = tr A`, `Σλ² = ‖A‖²_F` and `‖Av − λv‖`.

  **What it answers, and what it does not.** Given only the *open* structure of adenylate kinase,
  its softest mode overlaps the 7.1 Å motion the enzyme actually performs on closing by **0.799**,
  where a direction picked without looking scores `1/√642 = 0.039`. Binding a ligand cannot make any
  residue more mobile — a theorem about the Schur complement, not an observation — and the residue
  the ligand sits on loses 38.8% of its motion. It says nothing about binding affinity, it is not
  docking, it has no side chains, and it cannot travel along a large transition.

- **A thirty-first scene, and the first whose geometry is a measurement.** `29` names an STL
  somebody drew; `31-a-protein-shaking-at-body-temperature` names coordinates somebody refined.
  It is also the first scene with no resolution to sweep — a structure *is* its own resolution — so
  `Scene::refined` refuses it and the **cutoff** is swept instead.

### Fixed

- **A check the model passes that one line of geometry also passes.** Predicted fluctuations against
  deposited B-factors: `0.659` for lysozyme, `0.571` for ubiquitin. The distance from the centroid
  scores `0.699` and `0.804` on the same data, and the contact number `0.697` and `0.768`. **Neither
  geometric predictor loses to the model on any structure in the crate.** The test written first
  asserted that the model won; it failed, and the repair was to assert the ordering that was
  measured. A B-factor is one scalar per residue and cannot see a claim about direction, which is
  why the adenylate kinase comparison exists.

- **Two theorems asserted about the wrong object.** "Binding adds springs, so by Weyl every
  eigenvalue rises" is not Weyl's inequality — a ligand adds *nodes*, and mode 5 of a test helix
  falls from `0.204` to `0.178` N/m. And the protein's own statistics are not the `pp` block of the
  complex Hessian's pseudo-inverse: residue 4 came out `1.003` times as mobile after binding, which
  the theorem forbids. Through the Schur complement it is `0.828`.

- **A shipped scene reporting that nothing moved.** The CLI formatted body and field values with
  `{:.4}`, so anything under `5e-5` printed as `0.0000`. `15-a-hot-spot-in-a-block` was reporting an
  elastic displacement of 17 to 50 µm as zero. Measured across all thirty-one scenes, before and
  after: exactly three lines change and every one was showing a real value as nothing.

- **Three things the gate on one machine cannot see**, each found by CI after a green local run:
  `app/Cargo.lock` goes stale when the facade gains a dependency and the library gate never enters
  that directory; a closed form evaluated with `sin` is not the same number on two platforms, so a
  floor counting only the solver's rounding held here at exactly zero and failed on Windows at
  `8.88e-16` against `4.44e-16`; and `catch_unwind` compiles for `wasm32` and catches nothing there,
  because that target is `panic="abort"` — which `--no-run` cannot discover.

- **A gate log with two authors.** A stopped run's child kept going, a new run truncated the log
  under it, and a twenty-four-step gate produced forty-four step markers and both `the gate passed`
  and `THE GATE DID NOT PASS` in one file. None of the existing guards touch that: every individual
  check ran and reported its own exit code.

### Changed

- **The tally of unearned passes, which two documents were keeping separately.** `CONTRIBUTING.md`
  said eight while holding six rows, `CLAUDE.md` narrated a sixth through a tenth of its own and
  said that file "has all eight", and the instance that reached `main` — `4a3654f`, a `grep` over a
  test log still being written — was in neither list. **Sixteen**, indexed in one table, held
  against its own rows by a test, and `CLAUDE.md` no longer carries ordinals.

- **`AGENTS.md`'s "What is in the box" listed twelve crates out of eighteen.** Its first line says
  every name below it is re-exported through the prelude; `elastic`, `em`, `fluid`, `porous`,
  `shape` and `protein` were re-exported and not below it. A missing row states nothing, so no count
  guard could see it — the new test compares the table against `crates/` in both directions, and
  found the sixth on its first run after five had been counted by reading.

- **FRICTION is forty-eight findings**, forty-three fixed and five recorded. The new one is a
  constructor that panics on data the caller read rather than wrote: `Protein::new` refuses a
  disconnected network, and the scene format has to pre-check `components()` to report it as an
  error instead of dying.

## [0.20.0] — 2026-09-01

The same day as 0.19.0, which is unusual and deliberate: `Designed` gained a required field hours
after it was published, and a break is better corrected in the next version than left to accumulate
behind one.

### Added

- **`Placed::axis_angle`**, and the HTML report says which frame each panel is drawn in.

  The report draws one panel per card, so two domains never share a picture and each picture is
  right in its own frame. What was missing is any statement of *which* frame: every spatial axis is
  labelled in metres, and for a placed domain those were the domain's own metres presented as the
  world's. Two busbars 20 mm apart both read `0 .. 32 mm`.

  Stated rather than baked, for the reason the run format states it — an extent is a **box**, and
  the axis-aligned box around a rotated one is bigger than it. `at 20, -12, 12 mm, turned 90° about
  (0, 0, 1)`, in the scene format's own vocabulary, recovered through the quaternion the run
  carries. A domain at the origin emits nothing, so every report ever written is byte-identical.

  **`atan2` and not `acos(w)`, and the reason written beside it was backwards.** The editor, where
  this arithmetic was worked out, said `acos` loses precision near a half turn. At a half turn `w`
  is zero and `|d acos/dw|` is 1 — the best-conditioned point on the curve. It is the small angle
  that hurts: 1.1e5 at a thousandth of a degree, which is a nudge of a rotation ring. Measured
  4.6e-5 relative at 0.0001°, against 0 for `atan2`. The choice was right and the account of it was
  not, and the account is what sends the next person to the wrong end.

- **A schedule no domain could survive is refused when the scene is built.** `FRICTION.md`'s
  finding 8. `World::build` asks every evolving domain for `Domain::max_stable_dt` and refuses any
  non-subcycling schedule whose frame window is longer, so the editor — which checks as you type
  through the same `build_with` — says so while somebody is still writing the file.

  Half of it was already there: `verify::stability_hazard` asked the same question of the same
  built, unrun world. The two are one now, in `build`, with the battery's wording kept. Merging
  them exposed a gap in the new one — it checked only `staggered`, and `one-way` does not subcycle
  either.

  The suggested frame count is asserted to build, and one fewer to be refused, so it is a threshold
  rather than a plausible number. The message says it is checked from the **initial** state, because
  `max_stable_dt` is state-dependent and a green `--check` is not a guarantee.

### Changed

- **`mesh::Designed` gained a required `domain` field, and that breaks a 0.19.0 struct literal.**
  The one breaking change here, and it is why this release is the same day as the last.

  glTF and USD do not need it — a design is a node beside the panels there — but the HTML report
  draws **one panel per card** and has to know which card a design belongs in. Parsed out of
  `Designed::name` it would have worked today: the site happens to begin with the domain, and a
  `split('/')` is exactly the coupling that holds until somebody changes how a site is spelled.

  A caller constructing one adds `domain: <the panel's name>`. In this workspace that is
  `PlacedMesh::name`, which already carries it beside the site.

- **The HTML report draws the shape the scene designed, over the cells it became.** The exporters
  gained that in 0.19.0 and the report did not, so the one output a reader opens without a 3D tool
  was still showing the staircase alone.

  A **wireframe** and not a solid: the volume view is a raycast into an image buffer, and a shaded
  surface would cover the render it is meant to be compared against. Where the cells reach past the
  outline is the rasterisation, and that is the whole point of having both.

  The outline is projected through the raycast's own camera **inverted** — the ray builder turns a
  screen point into a direction, and a world point needs the other way round. `fwd`, `right` and
  `up` are orthonormal, so `su = b/a` and `sv = c/a` from three dot products.

  `Designed::outline` deduplicates by **position and not by index**: `mesh_surface` gives every face
  its own three vertices so each can carry a flat normal, so two faces sharing an edge share no
  index and an index-keyed dedupe draws every shared edge twice. By the bits rather than by a
  tolerance — the positions are the same `f32` copied from the same source vertex, so there is no
  question of how close counts as the same. The L-bracket: 14 points, 36 edges, 24 faces, and
  `V - E + F = 2`.

  `MAX_OUTLINE_EDGES` **refuses** past 20 000 rather than subsampling, which is the opposite of what
  `MAX_FACES` does and for a stated reason: a field's boundary at a stride is still that boundary at
  a stride, and every second edge of an outline is a picture of a different part. The caption says
  the face count instead.

  `Designed` states its `domain` rather than having it parsed out of the site. glTF and USD do not
  need it — a design is a node beside the panels there — but the report draws one panel per card
  and has to know which card.

  `Placed::apply` is public, so everything on the library side of the wire boundary shares one
  quaternion sandwich. `viewer-core` keeps its own copy, and that is the duplication this workspace
  pays for on purpose.

### Fixed

- **`tools/report-check` had an assertion that could not fail, and CI never ran the new ones.**
  Two holes found by sabotage rather than by reading.

  The outline check is conditional on the run carrying a design, so dropping the outline from the
  report's JSON made it assert **nothing** and pass. A JavaScript harness can check that the viewer
  drew what the data said; it cannot check that the data should have been there, because the page
  is its only input. That expectation is a Rust test now.

  And the `report viewer` job builds eight reports covering the six view kinds, none of which has
  `parts` — so every assertion about a designed outline was skipped in CI. Scene 29 is in that list
  now. A check nothing exercises is not a check, it is a claim.

  A first attempt at the drawing sabotage also "passed" for the wrong reason: emptying the design
  made the viewer throw, which fails loudly and says nothing about the assertion. Keeping the data
  and breaking only the loop is what produced `FAIL bracket-volume: the designed outline was drawn`.

## [0.19.0] — 2026-09-01

### Added

- **The exporters draw the shape somebody designed, beside the cells it became.** Every export of
  a designed part showed its *rasterisation* — a staircase on the solver's grid — while the thing
  it was rasterised from sat on disk beside the scene, reachable by nothing that writes a file.
  `Rasterised` has reported the volume error as a number since designed parts existed, and a
  number in a terminal is not what a designer looks at.

  Measured on scene 29, through OpenUSD:

  ```
  /World/bracket           3640 faces   [0, 0, 0] .. [0.05096, 0.05096, 0.020727]   the cells
  /World/bracket_parts_0_    24 faces   [0, 0, 0] .. [0.05,    0.05,    0.02    ]   the design
  ```

  The overshoot is the grid. It was always there and there was no way to see it.

  **The design is not in the run**, and that is the decision this took. A run is the simulation's
  output; an STL is its input. `mesh::Designed` carries the three ways of putting it in the run
  instead and why each is worse — repeated per frame it is unbounded, a static section is a new
  kind of thing in a format that has never had one, and a path resolved on read makes a run depend
  on a filesystem whose failure is an empty picture. So `mesh::Drawing` is an argument to the
  writer: `gltf_with` and `usda_with` take one instead of a `Surfaces`, which is the third
  positional argument they had gained in two commits and the last chance to reshape them before
  either is published.

  **It is drawn uncoloured, and that is a claim about honesty.** The solver ran on cells; tinting
  a smooth surface from the field would look like a higher-resolution answer than the one
  computed, which is the same mistake as renormalising a colour scale per frame in different
  clothes. `mesh::DESIGNED_GREY` is the grey the editor's viewport already used, chosen to be
  unlike every scale here — so the viewport and the export agree, and the test asserts the design
  carries *no* variation rather than asserting a shade.

  `editor_core::designed` is the reader, with the editor's viewport as its other caller. Its
  `PlacedMesh` now holds the STL's **own** triangles and states its placement, where it used to
  bake the pose into every vertex. Baking is not wrong for a vertex list — a rigid motion loses
  nothing there, unlike a box — but it was the second convention in a workspace that has just
  spent four commits removing one. `Placed::of` is public for the same reason: `capture` and the
  editor must not come to different answers about the same pose.

  `a_pose_moves_a_part_and_does_not_resize_it` asserted the old convention and now asserts the new
  one *and* composes the placement back to the same world box it always had — the same shape the
  `knows_no_physics.rs` change took, and the right kind of breakage.

  Six sabotages, all caught by an assertion: each writer dropping the pass, each writing the
  design at the origin, the colour going per-vertex, and the reader baking the pose in again. The
  fifth needed two goes — the first anchor matched twice, the harness ran the test anyway, and a
  pass with **no sabotage applied** was read as "survived". The check is gated on the edit now.

- **A shipped scene that states a pose**, which is the reason three consumers dropped one. Scene
  30: two blackened copper busbars crossing at a 4 mm clearance, identical but for their current,
  one turned a quarter turn about z and lifted above the other.

  Under the identity a domain's own coordinates and the world's are the same thing. So a reader
  that dropped the placement produced exactly the right answer for **all twenty-nine** shipped
  scenes, and glTF, then USD, then `viewer-core` each did — each found by hand, none by a test.
  This is the file that cannot agree with them.

  The arrangement is chosen so that *not* conducting is the design requirement rather than a
  limitation being worked around: two phases at a clearance must not touch, and `PoseSpec`'s doc
  already warns that placed parts do not interact. A scene of two parts in contact would have been
  the failure that doc exists to prevent.

  **The closed form is the steady balance** — everything generated leaves, and it leaves two ways:

  ```
  P = h A (T − Ta) + ε σ A (T⁴ − Ta⁴)
  ```

  solved for `T` by bisection. Not a second solver: one algebraic equation on a lumped bar where
  the run time-steps sixty-four cells, and the bar is isothermal to **3.4 mK**, so the peak cell
  *is* the mean. It matches to **1.6e-4** and **3.3e-4**, against a budget of 1e-3 that is derived
  rather than chosen: the run stops at 8.0 τ, so e⁻⁸ = 3.3e-4 of the rise is still to come, and
  the discretisation adds 1.5e-4. Bar B's error is the budget's leading term almost exactly.

  **The physics worth a scene** is that four times the heat is *not* four times the rise. Both
  bars are the same bar, so under convection alone the rises would be in the ratio of the watts
  exactly — 4, with no material property left in it. Measured **3.859**, because the hotter bar
  radiates as T⁴ and sheds disproportionately.

  How much radiation carries is asserted by conservation rather than by a threshold: the
  convective half is `h A ΔT` with ΔT the run's own answer, so the rest is radiative. **46.6%** —
  blackening a busbar nearly doubles what it can shed, which is why switchgear busbars are
  blackened. The first version of that check asked for the rise ratio to sit more than 0.1 below
  4, which is a number nothing derived; it tolerated ε down to 0.35 and had 1.4× of room. The
  conservation form fails at ε = 0.2 and reads 4% at copper's own 0.04.

  Five sabotages, all caught by an assertion rather than by an exit code: radiation not applied,
  `capture` dropping the pose again, the turn conjugated, the clearance closed, and the bars left
  bright. **The first run of that script proved nothing** — it backed both `pantometry-scene` and
  `pantometry-thermal` up as `/tmp/lib.rs.orig`, so restoring put one crate's source in the
  other's file, the workspace stopped compiling, and four sabotages "failed" without a test ever
  running. A non-zero exit is not evidence; the script requires `panicked at` now and says so per
  case.

  The scenes README gained a row for scene 29 as well, which never got one, and its opening count
  is corrected: it read "twenty-seven of them ... and one", which is twenty-eight and was never
  the number of files in that directory.

- **The run format states its version, and the version is read before the panels.** A scene has
  carried a `format` since it had consumers; a run never did, and the gap only mattered once
  something wanted to add a shape to it.

  `PanelData` is tagged and `deny_unknown_fields`, so a run carrying a panel kind a build has
  never heard of fails inside serde with a complaint about an unknown variant — which is what a
  truncated file says too. "This run is newer than your pantometry" and "this file is broken" are
  different sentences and a person can only act on one of them.

  So `Run::from_json` reads the `format` key **alone** first, from a reader that ignores every
  other field, and refuses a number it does not understand before parsing anything else. The
  scene format checks its version *after* parsing and has exactly this gap; this does not, and
  `the_version_is_read_before_the_panels` fails if the order is swapped back.

  **Backwards compatible**, because `Run` does not refuse unknown fields: a viewer built before
  the key ignores it. An absent `format` is 1, which is what every run this workspace has ever
  written is — the three recorded runs in `viewer-core/tests/runs` have no such key and still
  read.

  The writer's constant and the reader's are in different workspaces on purpose, since
  `viewer-core` deliberately does not link `pantometry`. `one_run_format_two_crates` is what
  makes them one number, and it checks the **bytes** as well as the constants: a writer that
  agreed on the number and spelled the key differently would pass a constant comparison and
  produce files nothing could version-check.

- **Rings to turn a domain with**, beside the arrows that move it. `Edit > Handles`, or `W` and
  `E`.

  `turn_about_axis` is `drag_along_axis`'s companion and deliberately its shape: a point on the
  ring moves with velocity `axis × (grip − origin)`, so the question is the one the translate
  handle already answers with the tangent in place of the axis. The test is the *rate* — halving
  the drag must quarter the miss — and the whole-radian secant, which is the mistake that made
  the translate handle move things 77% of the way, takes the ratio to 2.88 against 4.

  A turn **composes** rather than replaces. The file states one axis and one angle; a ring asks
  for a rotation about a world axis on top of whatever is there, and rotations do not commute.
  `compose_turns` multiplies quaternions, and it is checked by applying the result to points
  rather than by comparing axes and angles — `(a, 90°)` and `(−a, 270°)` are one rotation spelled
  two ways.

  `ring_points` produces the loop once and the shell both draws and hit-tests from it, because a
  control drawn somewhere other than where it is grabbed works everywhere except under the
  pointer.

- **The undo shortcut's decision is a function.** `edit::shortcut` takes the modifiers, the keys
  and whether anything has focus, and returns a direction. Six cases pinned, including a bare
  letter that must not reach the history and a focused widget keeping its own undo.

  Separated out because verifying `Ctrl+Z` by driving the desktop was tried and abandoned: three
  attempts raised the wrong window twice, missed a menu by a few pixels and clicked an outliner
  row, and the one keystroke that landed could not be read back. What the toolkit does with a
  press is egui's contract; which chord means which direction is this repository's, and that half
  is now checked.

- **A shipped scene whose geometry comes from a file.** Scene 29 names an ASCII STL — an
  L-bracket with a re-entrant corner and a chamfer, so a 2 mm grid has something to round off —
  and cools it to still air. The `parts` path had been in the format for several releases and
  **no scene used it**, so nothing in the twenty-eight-scene suite ran it; verifying the editor's
  new mesh drawing meant building a fixture in a scratch directory.

  ASCII rather than binary, because it is text: it reads in a diff, it greps, and it cannot carry
  anything a reader cannot see.

  The claim is a closed form **independent of the crate**. `Loss::volume_error` compares the
  filled cells against `Mesh::volume`, the divergence theorem over the same triangles the
  rasteriser read — agreeing with it would only say the crate agrees with itself. The shoelace
  area of the seven-point outline times the extrusion is a separate derivation: 1650 mm² × 20 mm
  = 33 000 mm³, so 4125 cells exactly, against 4100 filled. −0.61%, inside the 2% that
  `Loss::CLEAN_VOLUME_ERROR` already calls clean. Plus two bounds no cooling body may cross: the
  peak may only fall, and nothing may end below the air it is losing heat to.

- **Undo in the editor.** It had gained six ways to change a scene before it had any way to
  change one back — a value, a string, a domain added, a domain removed, a placement, and a drag
  on an axis — and two of those destroy something: a removed domain takes its whole block of
  JSON, and a drag rewrites a number nobody is necessarily looking at.

  `Edit > Undo` and `Ctrl+Z`, with redo on `Ctrl+Shift+Z`. Sixty-four states, oldest dropped.

  **The buffer changes two ways and only one of them passes through this crate.** The text pane
  is bound straight to the scene text, so typing never reaches a splice. A history that recorded
  only the splices would hold a snapshot from before the typing, and undoing a deleted domain
  would silently discard whatever had been typed since. So the current text is committed on
  *both sides* of every edit: the first call folds typing into a state of its own and is a no-op
  when there was none, and undo then walks back through the typing rather than over it.

  Committing before a **redo** as well, which the first draft did not. Typing leaves the buffer
  holding a state the history has never seen, and stepping forward from there would overwrite
  it; folding in first makes the typing the newest state, so redo correctly reports there is
  nothing ahead instead of reaching a future by discarding the present.

  egui's own fine-grained undo is left alone for whatever has focus. A focused text field owning
  its undo is what every editor does, and taking it over would mean reimplementing keystroke
  coalescing to no benefit — so `Ctrl+Z` is handled only when nothing has focus, the same guard
  the frame-stepping keys already use. The menu item works either way.

  A load from disk **resets** rather than appends: a file is a different document, and undoing
  from one scene back into another would offer a text belonging to no path.

  Ten tests, every case walked in both directions and three of them past an end, because the
  interesting bugs in a history are at its edges. Verified by breaking it three ways — trimming
  the bound without moving the cursor, letting a no-op commit become a step, and keeping the
  redo tail across a new edit — each caught by the test written for it.

- **An export can be the surface where a field reaches a value.** `pantometry_view::gltf_with`
  and `usda_with` take a `mesh::Surfaces`: `Boundary`, which is what they have always written,
  or `At(level)`. `pantometry run scene.json out.gltf --at 350` on the command line.

  `field_surface` draws the outside of the cells that hold a value, which for a solid block is
  the block — the same shape at every timestep whatever the values are doing. A level is the
  question a field is usually asked, and until now it could be seen in the editor's viewport and
  in nothing that left the workspace.

  **No format change.** The isosurface is derived from the `Field` panel a run already carries —
  counts, extent and values are all there — so this is a parameter rather than a shape. That was
  worth finding out before designing one: the two halves of "let the exporters draw a mesh" have
  very different costs, and only the designed-mesh half needs the run to carry anything new.

  `gltf` and `usda` are unchanged and delegate with `Boundary`, so no existing caller writes a
  different file. Pinned against the box rather than against itself: a distance field's cell
  boundary *is* its box, so the default's bounds must equal the extent while a level's are
  strictly inside. The first version compared the default to an explicit `Boundary` and passed
  when the default was redefined to a level — both sides moved together.

### Fixed

- **The editor's viewport drew a placed run at the origin**, over a wireframe that was in the
  right place. The third of three readers, and the worst of them: `editor-core` applies a scene's
  pose itself, so a placed part's outline was where it belongs and its colours were not — two
  frames in one picture, and nothing on screen saying which anything was in.

  `viewer-core` had parsed `place` since the key existed and read it nowhere. `Placed::is_here`
  was defined and never called, and `Panel::bounds` said **"in world coordinates"** while
  returning the panel's own box — a doc naming the wrong frame, which is worse than none, because
  the number is perfectly usable and silently in the wrong place.

  `Placed::apply` and `Placed::corners_of` are what every caller was about to write. Corners
  rather than a box: under a rotation the box *around* a placed cell is bigger than the cell, and
  not losing that is the whole reason the run format grew a placement. `field_shell` already
  mapped a unit cube through eight corners and never needed them axis-aligned, so a turned field
  is exact. `Panel::world_bounds` is the bounding box, for the one caller that wants one — a
  camera fitting to what is on screen — and is identical to `bounds` to the last bit under the
  identity.

  Eight sites: three in the shaded pass, three in the flat painter, the camera's world, and
  `segments`, which is the standalone viewer's half.

  **Measured with `--drawn-extent`**, a headless subcommand in the pattern `--layout-at` set: it
  loads a run, builds the editor's batches and inverts the framing to say what got drawn, in
  metres. A fixture holds one of each shape placed on a different axis — a field at x = 10, bodies
  at y = 20, a path at z = 30 — so one number says which site failed:

  ```
  baseline   solid x→11.00   solid y→21.25   lines z→30.0
  field unplaced             solid x→ 1.25
  bodies unplaced                            solid y→ 1.25
  path unplaced                                              lines z→ 0.0
  world from local boxes     world → [0 0 0  4.4 3.1 1.0]  — the built-in scene's room
  ```

  A fixture and not a shipped scene because **scene 30 has two fields and no bodies and no
  paths**: run against it, the body sabotage survived and said nothing, and the camera-world one
  was reported as caught on a difference between `-0.0` and `0.0`. A string comparison of two
  boxes is not a comparison of two boxes.

  **Not covered:** the flat painter's three sites, which have no headless hook. Said here rather
  than assumed, because an uncovered site is how this reached three readers.

  An earlier claim in this changelog — that `pantometry view` draws two placed domains on top of
  each other — was wrong. That shell draws **one panel at a time** and refuses fields and point
  clouds outright, so it would decline scene 30 rather than mis-draw it. The overlay is the
  editor's.

- **The USD writer ignored `Panel::place`, so the fix below reached glTF only.** The commit that
  added the placement said the defect was closed. It was closed in one of the two writers: a
  `Mesh` wrote `uniform token[] xformOpOrder = []` whatever the panel said, and a `BasisCurves`
  had no transform at all — a placed domain whose shape is paths was the same defect one prim
  type further along.

  Measured through OpenUSD rather than by reading the file. Two identical 20 mm blocks, one
  turned 60° about z and moved half a metre, exported and then loaded with `pxr`:

  ```
  before   near  [0, 0, 0] .. [0.02, 0.02, 0.02]
           far   [0, 0, 0] .. [0.02, 0.02, 0.02]      the same block, twice
  after    near  [0, 0, 0] .. [0.02, 0.02, 0.02]
           far   [0.482679, 0, 0] .. [0.51, 0.027321, 0.02]
  closed form    [0.482679, 0, 0] .. [0.510000, 0.027321, 0.020000]
  ```

  **USD's quaternion is not glTF's.** A `quatd` is written `(w, x, y, z)`, real part first, where
  a glTF node's `rotation` is `[x, y, z, w]`. That is not remembered, it is measured: `pxr` was
  asked for a 60° turn about z — every component distinct, so the text cannot be read both ways —
  and it wrote `(0.8660254037844387, 0, 0, 0.49999999999999994)`. `xformOpOrder` is outermost
  first, so `[translate, orient]` turns before it moves, which is a glTF node's composition too.

  This crate has no USD dependency and `usd_is_usd.rs` explains why. So the two writers are
  checked against **each other** on one frame instead: same rotation, each asked for the component
  the other format puts somewhere else. A permutation error differs by ~0.5 there, which is four
  orders above the tolerance — and the tolerance is the file's own `{:.6e}`, seven significant
  figures, half an ulp of 5e-8. An earlier 1e-9 was tighter than the format can be and failed a
  correct writer.

  **Still not honoured: the viewer.** `viewer-core` parses `place` into all three panel variants
  and nothing reads it — `Placed::is_here` is defined there and never called — so a window on a
  run draws two placed domains on top of each other exactly as `.usda` did. Named here rather
  than fixed here.

- **A run had no coordinate frame, and two placed domains exported on top of each other.**
  Measured, not argued: two blocks half a metre apart came out of glTF at the same coordinates,
  both `[0,0,0]..[0.03,0.03,0.03]`.

  `capture` held each domain's `Placement` and used half of it. A field was sampled in the
  domain's own coordinates and the pose was dropped — `let _ = pose` — while bodies had it
  multiplied into their positions. One file, two conventions, and nothing in it saying which a
  panel used. `Panel::bounds`' own doc said exactly this and declined to promise anything about
  overlaying them, which is a hazard written down rather than closed.

  Nothing caught it because **no shipped scene states a pose**. Under the identity, local and
  world are the same thing, and all twenty-nine agree.

  `Panel::place` carries a `Placed` now — a translation and a unit quaternion in `[x, y, z, w]`,
  which is a glTF node's `rotation` in a glTF node's order, because that is where it goes. Bodies
  stay in their domain's frame like a field's samples, and the exporter puts the transform on the
  node. Run format **2**; the key is written only when it is not the identity, so every scene here
  produces the bytes it did.

  **A rotated cell is exact again.** Bodies used to pose the cell's two corners and take the min
  and max of the result, which for any rotation is a box *around* the cell — a quarter turn of a
  unit cube gave one √2 across, and nothing written down could recover the cube.

  `knows_no_physics.rs` asserted the old convention and now asserts the new one *and* that
  composing the placement with the local position lands the mote exactly where it always was:
  what changed is where that fact is stored, not the fact.

- **A level nothing reaches said the field was empty.** Reusing the boundary's message meant
  `--at 5000` reported "no cell that holds a value" about a field entirely full of them, sending
  a reader to look for a void instead of correcting a number. It now says which level was asked
  for and what the field actually spans: *"bracket never reaches 5000: the field spans 391.2002
  to 391.2416"*.

- **A `parts` path was resolved against the working directory, not the scene.** Invisible until
  a scene shipped with a part in it, and then a trap: scene 29 ran from `app/pantometry-world`
  and failed from the repository root *and* from its own directory, with an error naming a path
  the user never typed.

  `Beside` resolves a relative `stl` next to the file that names it, which is what every format
  that refers to a sibling does. The CLI and the editor both build one from the scene's own path,
  so a scene checked at a terminal and a scene opened in the editor read the same bytes.
  Measured from four directories, all of which now agree on 4100 cells; before, three of the four
  could not find the file at all.

  An absolute `stl` is used as written, and the error names the path **as the scene spelled it**
  followed by where that was looked for — reporting only the resolved path answers a question the
  reader cannot map back to their file.

## [0.18.0] — 2026-08-28

### Added

- **A designed mesh reaches the screen.** `pantometry_view::mesh::mesh_surface` turns triangles
  into a `Surface`, and the editor draws every `parts` entry of a scene in the place its voxels
  will occupy. Until now the only picture of a designed part was the staircase it rasterised
  into, and a staircase is a picture of the grid as much as of the object.

  It draws **before there is a run**, which is the point: everything else in that viewport is a
  picture of results, and a part is a picture of the scene. Once a run is loaded the field's
  surface wins for that domain, because two surfaces in one place z-fight into a picture of
  neither.

  `pantometry-view` still does not depend on `pantometry-shape`. `mesh_surface` takes the same
  plain arrays the module's other three producers take, and `editor-core` — which links both —
  converts and applies the pose. The **same** pose `Voxels::onto` rasterises against: a mesh
  drawn a millimetre from its own voxels is two pictures that disagree, and nothing on screen
  would say which to believe.

  Flat-shaded, per face, three vertices to a triangle. An STL has no curvature to preserve, only
  the tessellation of one, so smoothing across a machined edge would round a chamfer that is not
  there. A triangle with no finite positive normal is dropped — a zero-area face has no direction
  and a non-finite vertex poisons the cross product to `NaN`, which shades a face black.

  Checked against a cube's `6s²` and `s³`, both computed from the emitted surface rather than
  from the triangles that went in, and against a quarter turn worked out by hand rather than by
  the code under test. Verified by breaking it four ways: negating the normals, doubling the
  positions, removing the degenerate guard, and collapsing every normal to one direction — each
  is caught, and by the test written for it.

  This closes the second half of `ARCHITECTURE.md`'s "a renderer with depth". What a mesh still
  cannot do is travel in a **run**: the report and the glTF and USD writers work from a `Run`,
  which has no shape for a surface, and adding one is a wire-format change to a reader that is
  `deny_unknown_fields` on purpose.

- **An isosurface: the surface where a field reaches a value.** `pantometry_view::mesh::isosurface`.
  `field_surface` draws the outside of the cells that hold a value — a staircase, and a picture of
  the grid as much as of the object. This draws where the *values* reach a number, which is the
  question a field is usually asked: not where the block is, but where it is 100 degrees.

  `ARCHITECTURE.md` has said for several releases that depth buffering "becomes worth doing when
  something here has *surfaces* — a mesh, an isosurface — and nothing does". Something does now,
  and the editor's viewport has had the depth buffer waiting for it.

  **Marching tetrahedra rather than cubes.** Marching cubes needs a 256-entry case table with
  ambiguous entries, where the wrong choice leaves a hole; six tetrahedra per cell have sixteen
  sign patterns and two shapes among them, no ambiguity, and watertightness by construction because
  neighbouring cells cut their shared face the same way. The cost is more triangles, which is the
  right side to be wrong on: a hole in a closed shape reads as a rendering artefact rather than a
  bug, so the cheaper method fails in the way nobody investigates.

  Checked against closed forms, two of them **rates**: a sphere's area converges on 4πR² and its
  enclosed volume on 4/3πR³ at **second order** across three refinements. A plane comes out exactly
  flat, because a linear field makes the edge interpolation exact rather than approximate. The
  enclosed volume is the divergence-theorem sum, which is also the winding test — faces pointing
  both ways cancel toward zero.

  A cell touching a non-finite sample is skipped whole. `Solid3D` returns `NaN` for an emptied cell
  deliberately, and interpolating an edge that ends in one would invent the number the `NaN` exists
  to refuse.

  In the editor as **View → Isosurface at a value**, with the level from the run's own scale rather
  than the frame's — a slider that renormalised per frame would move the surface while the field
  stood still. `field_shell` takes a `level: Option<f64>` rather than gaining a twin, because
  everything after the mesh is the same and that part has been wrong before: a 40 mm cube once
  exported 80 mm across.

### Changed

- **`verify` measures all three of `advance`'s refusals.** It measured one, and said so about
  itself: *"Neither margin is measured here yet, and the report heading says so rather than
  letting 'conservation margins' claim more than it covers."*

  The reason was not laziness, it was that the other two **cannot be seen from outside `advance`**.
  The whole-simulation audit compares the ledger before and after, which a caller can redo. The
  transfer audit reads what is left on the bus, and the per-domain books check snapshots one
  domain's ledger across its own turn — both are gone by the time `advance` returns. So the kernel
  reports them: `Report::transfer` and `Report::books`, each a `Margin` of what the check measured
  and what it was judged against.

  **The two are not in the same units and the report says so.** A transfer margin is an absolute
  amount, because that check is absolute — what is left on a channel *is* the scale, since all of
  it went missing. A books margin is a ratio on the domain's own holdings, which is the entire
  reason that check exists. Ranking either against the other is meaningless; each is read against
  its own tolerance. "Worst" is likewise `residual / tolerance` rather than raw size, or a loose
  quantity would always outrank a tight one that is nearly out.

  Pinned by the case it exists for: a domain that loses `1e-10` of itself against a tolerance of
  `1e-9` **passes**, one digit from refusing, while the whole-simulation audit over the same step
  reports more than eight digits of headroom. The test asserts the *gap* between the two numbers,
  not either alone. Visible on a shipped scene too — scene 22 reads 4.5 digits on the sum and 4.3
  on `buffer/energy`.

  `Exchange::audit_transfers` is untouched; the measurement is a separate `worst_undelivered`,
  because a check that also measures is a check somebody calls for the measurement and gets a
  refusal from.

## [0.17.0] — 2026-08-27

### Changed

- **A domain can be dragged in the viewport.** Three translate handles on the selected domain,
  and the arithmetic behind them — `editor_core::drag_along_axis` — lives in the GUI-free half
  because it is a function of a camera and two vectors, and because its error has an *order* a
  test can measure.

  That mattered immediately. The first version took the handle's screen direction as the secant
  across a whole unit of the axis rather than the derivative, and the object went **77%** of the
  way the pointer asked — at every drag size, converging to a fixed 23.5% relative error instead
  of to zero. A single measurement looks like ordinary discretisation and a tolerance would have
  been written around it; the order creeping toward 1 where a correct inverse gives 2 is what said
  the map was wrong rather than coarse. It measures 2.00 now. The companion test — "a drag along
  the handle moves one unit" — passed the whole time, because its helper computed the direction
  the same wrong way.

  A handle pointing at the camera returns zero rather than turning a one-pixel twitch into a leap
  across the scene, the camera does not turn while a handle is held, and releasing over empty
  space does not clear the selection.

  **The framing holds still while a domain is dragged**, and keeps holding until the view is
  fitted again. Every projection in the viewport is relative to the centre of everything the scene
  contains, so moving the only domain in a scene moved that centre by the same amount and the box
  did not appear to move at all — measured on scene 15: the box goes 3 mm, the centre goes 3 mm,
  apparent motion `0.000000000 m`. It felt like the box resisting the pointer, and it was the
  camera following it. Released on mouse-up the picture would snap back by the whole drag instead,
  which is the same defect arriving late, so the camera now moves when it is asked to.

- **A domain can be turned from the editor.** `PoseSpec` has carried a `turn` — an axis and
  degrees — since it arrived, and nothing in the editor could write one.

  It turned out to be much smaller than planned, because of what a measurement said before any of
  it was built: the inspector's generic walk **already reaches into `poses`**, so a turn is
  editable the moment the key exists. The whole job was creating it, which is `set_pose`'s pattern
  one level deeper. A test pins that, so a walk that stopped reaching in would fail rather than
  making the control silently read-only.

  Checked against a closed form, because the obvious assertions are blind: creating the key,
  reading it back and coexisting with a position all pass if `set_turn` writes somewhere the
  builder never looks — the trap five pose tests fell into. **360 degrees is the identity** for
  any axis, so the bounds return; **90 degrees is not**, so they must move. Neither alone is
  enough. Order is checked both ways too, since writing the position first and the turn first take
  different branches at all four levels.

  An axis of `[0, 0, 0]` is refused before the file is touched, for the reason `PoseSpec::to_pose`
  gives: normalising a zero vector yields a `NaN` rather than an error.

  The handles still only translate; an arc drag is a different inverse and is not here.

- **A domain can be moved from the editor.** The inspector grew a *placement* control, and it is
  the first write in the editor that **creates** a key rather than replacing one: `poses` is a map
  beside `materials`, and **no shipped scene states it** — zero of the twenty-eight — so moving
  anything means writing something the file does not have. Three levels can be missing at once —
  the `at_m` array, the domain's entry, the whole `poses` object — and each is written where it is
  found, appended after the last member at the object's own indentation. A scene that gains a pose
  keeps every other byte, and one that already had it has three numbers replaced rather than a
  second `at_m` appended.

  The drag rate comes from the scene's own extent rather than a constant, because a millimetre a
  pixel is right for a die and wrong for a room. It is still a rate and not a limit: nothing bounds
  where a domain may go, for the same reason nothing clamps the other fields.

- **A domain can be added and removed from the editor.** Before this, creating one meant typing
  an object into the JSON pane from memory — the format defines nineteen kinds and the editor
  offered none of them. **Domain → Add** inserts a starting example, named so it does not collide
  with what is already there; **Domain → Delete** takes the selected one out with the comma that
  held it in place.

  **Where the nineteen examples live is the whole design question**, because Rust cannot enumerate
  an enum's variants and a hand-kept list of that size is wrong between an addition and somebody
  noticing. They live in `pantometry_world::templates`, beside the format, and a test compares
  them against a set maintained somewhere else entirely: the kinds the twenty-eight shipped scenes
  actually use. That works because the two sets are equal today — nineteen variants, nineteen
  distinct kinds — which was measured rather than hoped, and the comparison runs **in both
  directions**, so a template with no scene fires one half and a twentieth domain with a scene
  fires the other. The same shape as `counts_in_prose`.

  They are hand-formatted JSON text rather than serialised values, for the same reason the rest of
  the write path is a byte splice: serialising would insert an alphabetised object into a file
  written `kind` and `name` first.

  Three things the extraction found by trying rather than assuming. Every `structure` in every
  shipped scene layers materials its own scene *declares*, so none of them stands alone and that
  one template is hand-written with the borrowing regions dropped. And two kinds cannot be
  complete on their own at all — a `beam` states what it shines `onto`, a `structure` the block it
  `follows` — which the format catches at two different times: `follows` is refused by the build,
  while a dangling `onto` **builds** and is stopped by the conservation audit at the first step
  with "published but not consumed". Both are pinned by tests, so a third reference arriving is a
  failure rather than a surprise in the menu.

  Adding each of the nineteen and removing it again returns the file byte-for-byte, which is the
  strongest thing either operation can promise and holds only if the insertion takes exactly the
  bytes the removal gives back.

- **The editor's inspector writes.** Counted rather than assumed: across 2101 lines of editor
  shell there were exactly **two** widgets that could change anything — the file path, and one
  multiline text box holding the whole scene as raw JSON. The inspector showed values and could
  not alter one, so every edit was a JSON edit.

  Selecting a row that names a domain now lists every value the scene states about it, and a
  change goes into the text and re-checks on the same frame. Three rows resolve to the same domain
  — its extent before a run, its field and its readings after one — because after a run the reader
  is looking at the last two.

  **What to cover was counted, not guessed.** A census of every value in the twenty-eight shipped
  scenes: 358 numbers, 237 strings, 59 arrays of numbers, 46 arrays of objects, 44 nested objects,
  one array of arrays, and **no booleans at all**. Three things followed. Strings are two fifths
  of the file, so a numbers-only inspector cannot change a material. Nesting is ordinary rather
  than exotic — a `hot_spot`'s `above_k`, a `region`'s `material` — so the walk is a full one and
  not the single level it started as. And there is no flag widget, because the format has nothing
  for one to point at. `material` gets a menu built from the catalogue **plus whatever the scene
  declared**, since a menu that offered only the catalogue could not express a file that names its
  own substance.

  The walk is uncapped, also measured: the widest domain in any shipped scene yields 55 fields and
  the median is 8. A cap would have to drop rows, and a panel that silently omits part of the
  scene it describes is this repository's oldest failure shape.

  **`kind` and `name` are held back, and the second was a defect caught before it shipped.**
  `kind` because the format is `deny_unknown_fields` and changing it makes every key beside it
  unknown at once. `name` because the outliner's selection is keyed by it — a text field
  committing per keystroke would rename the domain on the first character, the row would stop
  existing, and the field being typed into would vanish under the cursor. It is also an
  identifier: five other keys and every key of `poses` refer to a domain by name, so a rename is a
  multi-site edit the inspector cannot make atomically. Both stay jobs for the text box.

  **It is a byte splice, not a re-serialise, and the reason is measured.** `serde_json` without
  `preserve_order`, which this workspace does not enable, backs an object with a `BTreeMap`: a
  `Value` round trip alphabetises every key. The shipped scenes are hand-formatted with `kind` and
  `name` first, so one drag would have reordered and reflowed the whole file. `set_number` and
  `set_text` replace the bytes of the one value and copy the rest verbatim — a string goes through
  `serde_json` on the way in, so a quote or a backslash cannot end it early — and the test asserts
  that literally
  rather than by comparing parsed values — a round trip that produced an equivalent scene would
  pass a value comparison having destroyed the formatting.

  A literal's kind is read off the file, not from a list of key names: `frames` is a `usize` and
  `11.0` is a scene that stops parsing, so a value written without `.`, `e` or `E` stays whole and
  a fraction dragged onto one is refused; a value written with them keeps its decimal point even
  when it lands on an integer.

  Nothing enumerates domain kinds — every value in the domain's object is offered, at whatever
  depth, whatever the object is, which is ARCHITECTURE.md's rule for the inspection half. Nothing clamps either: the
  scene's own check is the authority on what is legal and it runs on every change, measured at
  **0.0–0.4 ms** across the shipped scenes, so a limit invented in the shell could only be a second
  opinion that refuses a value the format takes.

- **`verify` reads the rasterisation report instead of printing it.** `ARCHITECTURE.md` listed
  rasterisation loss among the four errors a passing audit cannot see and said "nothing yet turns
  the measurement into a verdict a person reads". `--check` had printed a `Loss` per designed part
  since parts existed, and a printed number is one somebody has to notice.

  The battery now carries every part's measurement as a row and turns the ones the grid could not
  hold into findings, which are an exit code. The verdict is `Loss::is_clean`'s — the 2% bar moved
  out of the comparison into `Loss::CLEAN_VOLUME_ERROR` so the finding can quote the number that
  was actually applied rather than a second copy of it.

  **The case that had nothing watching it is a part filling zero cells.** `World::build` refuses an
  assembly where *every* part vanished; one of two coming out absent built, ran, conserved, drew
  and reported — a correct answer about an assembly with a piece missing. Measured on a 20 mm block
  at a 2 mm cell: a 12 mm brick takes 150 cells and exactly its own volume, and a 0.4 mm plate
  beside it takes none, because the nearest cell centres are at 13 and 15 mm. That scene is now one
  finding rather than a clean pass.

  The rows print whether or not anything was lost, for the reason this repository keeps
  rediscovering: a section that appears only on failure cannot be told apart from a check that did
  not run. Scenes with no `parts` — which is all twenty-eight of the shipped ones — print no
  section, because a heading over an empty list is a claim to have checked geometry the scene does
  not have.

- **The HTML report shows numbers, not only pictures.** Every view had a colour bar and a cell
  count, which is enough to see *that* something happened and not enough to say *what*. Four
  separate gaps, all now closed:

  - **Axes were in cells.** A room was "61 x 43" and never "4.4 m by 3.1 m", so a hot spot could be
    seen and not located. `PanelData::Field` gained `extent_m`, the box the field was sampled over,
    and every spatial axis is labelled from it — one unit per axis, chosen from the span, so a room
    is metres throughout and a 0.5 mm channel is micrometres throughout.
  - **The scalar chart had no y-axis at all.** Seventy-four pixels were reserved for labels that
    were never drawn, and every series was normalised to its own range, so six of them filled the
    frame and two temperatures forty kelvin apart drew the same line. Series now share an axis when
    they share a unit, with ticks; the legend switches them; and a reading missing from a frame is a
    break in the line rather than a zero, which is what it used to be plotted as.
  - **Nothing could be pointed at.** Hovering reads back the sample under the cursor — its index,
    its position in metres and its value.
  - **The bodies view had no colour bar**, alone among the six. It has one, and every colour bar now
    carries tick marks rather than only its two ends.

  Also: previous/next frame and a speed control, with space, the arrow keys, Home and End; canvases
  sized from the element and the display's pixel ratio rather than from a fixed 1400 baked into the
  markup; a print stylesheet; and only the views that changed are redrawn, so rotating a cube no
  longer re-runs a raycast for a volume nobody is touching.

- **The colour scale is constructed rather than chosen** — `pantometry-view::ramp`, built in CIE LCh
  with lightness linear in the value. The four-stop sRGB gradient it replaces folded back on itself:
  it was brightest at 0.67 of the range and darker above that, so **89 of its 255 steps had a larger
  value looking darker than a smaller one**, and the top third of the scale printed as the same grey
  as the middle. A ten-per-cent change of value also looked 3.07 times bigger in one part of the
  range than another; it is 1.85 now, and the remainder is the shape of the sRGB solid rather than
  the construction. A field whose range straddles zero gets a diverging scale with its neutral at
  **the value zero**, not at the middle of the range, and its two arms mirror in lightness and
  chroma. The table the browser uses is generated by the same code the tests measure.

  The volume renderer had the same off-by-a-datum in its opacity curve and it was invisible: the
  curve is `|2t - 1|`, transparent in the middle, and `t` was the position in the *range*. For a
  field running -100 to +300 that made **+100** the transparent value and drew zero as solid. It
  is transparent at zero now, which is what the comment beside it always said it was for.

- **glTF exports in metres.** It wrote grid indices, with a comment explaining that the extent was
  not in the frame to write, so a 9x9x9 block arrived in Blender nine metres on a side whatever it
  was. The sample positions are also where the samples were, rather than half a cell off: `capture`
  samples corner to corner and the exporter was adding 0.5 of a cell to every axis.

- **`to_json` writes magnitudes.** Its timestamps were `{:.6}`, so a four-nanosecond run wrote two
  hundred frames all at `t = 0.000000` — the same defect `readings_csv` was fixed for, still present
  in the writer beside it. Both writers share one number encoder now, which also drops trailing
  zeros that carry nothing: 8.3% off the largest report this workspace produces and about 1% off a
  typical one. The estimate first written here was "5.21 MB to 1.48 MB" and was a guess; the table
  in `data.rs` is the measurement.

### Added

- **USD.** `pantometry_view::usda` writes a whole run as a `.usda` stage — `pantometry-world
  scene.json out.usda` — so usdview, Omniverse, Houdini and Maya are the viewport.
  `ARCHITECTURE.md` had named this "the next rung, and only worth it if the readers really are
  Omniverse", and that condition is what made it wait.

  It does what glTF structurally cannot. USD has **time samples on any attribute**, so one file
  carries the topology once and the colours per frame, and the timeline scrubs the physics rather
  than a camera move. A body's positions animate too. Verified through Blender's OpenUSD: scene
  23's part is cream at frame 0 and green at frame 20 while its lid warms from navy, both on one
  scale across the run so the two frames are comparable.

  It carries the **numbers**, which no picture reaches. Twelve of the twenty-eight shipped scenes
  have a domain with no field and no bodies, and for several of them the scalar *is* the result —
  so every domain's readings go out as time-sampled custom attributes under `pantometry:`.
  `ARCHITECTURE.md` is right that there is no USD schema for a `Ledger` or a `Violation`, and this
  invents none: a custom attribute is a number with a name, which is what a reading is.

  Still **no dependency**. `.usda` is USD's text serialisation and this writes it by hand, the same
  choice glTF and SVG were written by hand for; `usdcat -o out.usdc out.usda` is one command if a
  pipeline wants the binary crate. `metersPerUnit = 1` is stated rather than defaulted, because
  USD's fallback is centimetres and a file that omits it has every part read a hundred times too
  small.

- **`pantometry_view::mesh`** — the geometry, once, for both exporters. A solid coming out one size
  from the glTF and another from the USD would be worse than either being wrong, because the two
  disagree and nothing in either file says which to believe. `the_gltf_and_the_usd_describe_the_same_solid`
  is that claim.

### Fixed

- **The USD wrote a primvar's interpolation where USD does not read it.** `interpolation` is
  attribute **metadata**, in parentheses; it was written as a separate
  `primvars:displayColor:interpolation` property, which is the spelling a primvar's `:indices`
  uses. With it missing the attribute falls back to `constant`, so a file carrying a colour per
  vertex was drawn with the **first one** over the whole prim.

  Found by rendering the file, where scene 23's hot part and cooled lid came out the same colour
  while the text plainly held two. No check on the document could have said so, which is why this
  one now checks for the metadata form and for the absence of the wrong one.

### Changed

- **glTF exports surfaces, not a point cloud.** A point cloud has no silhouette, takes no light
  and casts no shadow: it is a picture of the *sampling*. A three-dimensional field is now the
  boundary of its present cells — one quad per face whose neighbour is absent or outside — so a
  solid 9x9x9 block writes **486 quads instead of 4,374**, 89% of them being between two cells and
  unseeable, and a void inside a part produces a real interior surface. Bodies are spheres up to
  256 of them. Everything carries normals.

  Rendered from Blender straight off the export: scene 23 is a hot part and a cooled lid with a
  real gap between them, and scene 16's room shows its standing wave's nodal planes as dark bands
  across a solid. Neither was visible before.

  A 2D field is geometry now too. It was refused as "a plane of samples is a graph", and the
  reasoning had expired — it was a graph only because the run recorded no box to give it a size,
  and `extent_m` arrived. A 1D field is still refused.

  The trade is stated rather than hidden: **a surface hides the inside.** A hot spot in the middle
  of a block exports as a block whose faces are all at ambient, because they are. The interior is
  what the report's raycast and slice montage are for.

- **Choices the export makes are reported.** `Exported::notes`, beside `skipped`: a sphere's radius
  is **not in the data** — a body set carries positions and a value, not an extent — so it is a
  quarter of the median nearest-neighbour distance and says so; a subsampled surface says so; and
  the colour scale spanning one frame rather than the run says so, because every other view in the
  crate spans the run.

### Fixed

- **glTF colours were sRGB in a linear attribute.** The specification puts `COLOR_0` and
  `baseColorFactor` in linear space and only textures in sRGB. This wrote `byte / 255` from an sRGB
  ramp, so a mid-grey went out as 0.5 where linear wants 0.216 — every colour this workspace has
  ever exported was decoded about **2.3x too bright in the midtones**, uniformly, and looked like a
  plausible picture the whole time.

- **An exported field was one cell larger than its extent.** Giving every sample a full cell either
  side made a 40 mm cube 80 mm: `capture` samples corner to corner, so the first and last samples
  sit *on* the faces of the box and an end node owns **half** a cell. Cells are clipped to the
  extent now. That arithmetic — a boundary sample owning half a cell — is the third defect of its
  kind here, after `Tube` and `Room` divided a wall's divergence by a whole `dx`.

- **The export held a fifth copy of the four-stop colour ramp**, the one whose lightness folds back
  on itself at 0.67. It is `pantometry-view::ramp` now, like everything else, and a signed field
  gets the diverging scale with its neutral at the value zero.

- **The editor's viewport painted a quarter of every volume in the wrong order.**
  `Camera::project` returns a depth — *distance from the eye, larger is further* — and the native
  shell discarded it, then recovered a stand-in by projecting each point a second time one
  millimetre along world z and taking the reciprocal of how far the two landed apart on screen.
  That is proportional to how far **off the view axis** a point is, not to how far away. On a
  6x6x6 lattice at the camera the app opens at it ordered **6,026 of 23,220 pairs backwards —
  26%** — worst at the centre of the screen, which is where the object is. The browser shell had
  always sorted by the real depth; both call `editor_core::far_to_near` now, so there is one copy
  left to be wrong.

- **The viewer refused every run file the library wrote.** Its reader is `deny_unknown_fields`,
  deliberately, so the moment `pantometry-view` began writing `extent_m` on a field it stopped
  parsing anything. The library's twenty-step gate could not see it — separate workspace, separate
  CI job, and nothing in `crates/` reads that format back. `extent_m` is now `Option` with a
  default, and the committed fixtures are deliberately left in the old format because they are the
  only proof the old direction still works.

- **Two CI jobs were already failing on `main` before any of this** — `native viewer` and
  `gpu accelerator`, both on `cargo fmt --check`, in `runtime/viewer/viewer/src/main.rs:74` and
  `runtime/gpu/src/lib.rs:191`. Both formatting drift; both invisible to the gate, which formats
  only the library's workspace. Fixed, and the general point recorded: a green gate is not a green
  CI, and this tree has five workspaces.

- **The editor printed readings `{:.4}`**, so a cavity holding 3.19e-10 J showed `0.0000` beside a
  field the same run reported at 921 V/m. The third writer in this workspace to erase a magnitude
  that way, after `readings_csv` and `to_json`.

- **A field is drawn where and how big it is.** `Panel::bounds` returned the grid in **cell units**
  for a field — a 9x9x9 block framed as a nine-metre cube, at the origin however far from it the
  part really was — because the run file recorded a grid and not the extent. It records the extent
  now, so both shells prefer the run's own box and fall back to the scene's placed one.

### Changed

- **The editor's viewport is an instrument.** A colour bar with numbers on it, so a colour can be
  read back to a value; a scale bar in model metres, so a wireframe says whether it is a 40 mm die
  or a 4 m room; a probe that names whatever the cursor is over — body index or grid cell, with
  the value and its unit — picking the nearest **to the eye** rather than the nearest on screen,
  so the thing named is the thing visible; and a transport with play, step either way, and the
  space bar, arrow keys, Home and End.

- **One colour scale across the workspace.** Both editor shells carried their own copy of a
  straight blue-to-red line in sRGB, which made four spellings of "cool to warm" here. It covered
  **16.6 L\*** against the library scale's 74, and seventy-five of its 255 steps ran backwards.
  `editor_core::value_colour` is `pantometry::view::ramp`, and a range that straddles zero gets the
  diverging scale with its neutral at the value zero.

### Added

- **`tools/report-check`** — the report's inlined viewer, executed against a stub canvas, with
  assertions on what it drew. It had never been run: every test asserted on the HTML as a string, so
  a syntax error or a view that silently drew nothing would have shipped with every assertion still
  passing. Its README records why the first version of the harness was itself the bug — a
  `vm.runInContext` measurement said the volume renderer took 118 ms a frame and it takes 4, because
  a vm context costs 30x on hot numeric code.

## [0.16.0] — 2026-08-20

### Changed

- **The project is called `pantometry`.** Everything through **0.15.0 was published as `dualis`**, on
  crates.io as seventeen crates and on PyPI as one — and those names are permanent, so this is not a
  rename so much as the same project continuing under a name that works.

  `dualis` collided and did not find. It is a common Latin fragment, and a search for it does not
  arrive here. The new name is a word two books have used and nobody has taken: *Pantometria* (1571)
  and John Dawes's *Pantometry; Or, an Attempt to Systematize Every Branch of Admeasurement* (1830).
  It says **the measurement of everything**, which is both what this library is for — a digital twin
  of physics, chemistry, biology and engineering — and how it answers: by measuring rather than
  predicting. The cell-size chooser rasterises every candidate instead of extrapolating; a view
  factor written down as a known error was measured and turned out exact.

  Checked before it was chosen: free on crates.io and PyPI, and no company or trademark that a
  search surfaces. Six earlier candidates failed exactly there — `clapeyron` is a Julia
  thermodynamics toolkit, `equipoise` is a live pharmaceutical mark, and `conserva`, `virial`,
  `adiabat` and `holonomy` are all company names.

  **Nothing else changed in this release.** The code is 0.15.0's, with 1,971 occurrences of one word
  replaced in 263 files. The `dualis-*` crates are yanked and point here.

## [0.15.0] — 2026-08-19

### Added

- **Geometry a person designed can drive the physics: `pantometry-shape`, the sixteenth crate.** `Mesh`,
  `Triangle`, `Voxels` and `Loss`. An STL in, and out comes `|i, j, k| voxels.contains(i, j, k)` — which
  is *already* the signature `Solid3D::fill`, `Block::fill` and `Waves::fill` take, so **no domain
  changed**. The crate depends on `pantometry-units` and nothing else in the workspace, and nothing in the
  workspace depends on it except the facade.

  Until now the only way to say what shape a thing was, was a closure written by hand. That is a fine
  answer for a test and no answer at all for someone who has a part.

  **The crate's real output is the `Loss` report, not the cells.** A 0.5 mm rib voxelised at 2 mm is not a
  thin rib, it is gone — and the simulation then runs perfectly well, conserves energy, audits clean, and
  answers a question about a different object. There is no symptom. So every rasterisation carries a
  signed volume error, the fraction of the volume in boundary cells, a count of runs one or two cells
  thick, and a count of rows the rasteriser could not decide. A 0.4 mm plate at 2 mm reports a **400%**
  volume error and 400 thin runs rather than quietly becoming a 2 mm plate.

  **The volume error is not a convergent quantity, and that is the finding.** The obvious test — "voxel
  volume converges to mesh volume at first order" — is wrong, and it would have passed. Cells the surface
  bulges out of cancel against cells it cuts into, and what is left is a lattice-point count: a 10 mm
  sphere at 2.5, 2.0 and 1.5 mm gives `+4.9%, +5.8%, −2.3%`, so the first refinement makes it **worse** and
  the second flips its sign. Sliding the mesh a third of a cell off the grid changes none of it. On powers
  of two the sequence looks tidy, which is how three points could have been fitted to any order you liked.

  What does behave is the boundary layer, added as `Loss::boundary_fraction`: a surface area rather than a
  cancellation, so it is `A·dx/V` — first order with a coefficient, measured to halve at 1.897×, 1.892×
  and 1.975× per halving. Its coefficient is asserted as the bracket `[0.5, 1.5]` that counting column
  ends on a convex body *proves*, not as the 0.82 that is measured, because 0.82 is a property of the
  staircase that nothing here derives. On a box the same layer is exactly `1 − (nx−2)(ny−2)(nz−2)/nx·ny·nz`
  and is checked with no tolerance at all.

  It is also the number to put in front of someone choosing a cell size: *forty-three percent of your
  object's volume is in cells that could have gone either way* says what 2 mm on a 20 mm ball means in a
  way a step count does not. It is deliberately left out of `Loss::is_clean`, because no threshold on it
  is right for more than one physics.

  **`Loss::retried_rows` exists because a review asked what would notice if the retry did nothing.**
  Nothing would. `ambiguous_rows` — the case where every perturbation fails — has never been produced by
  any mesh in the suite, so the whole degenerate-ray path was exercised only by inference from a cell
  count, and a retry that silently no-opped would have looked exactly like a mesh that never needed one.
  Counting the rows that *were* retried makes it visible: a cube reports 8 at 1 mm, 20 at 0.5 mm, one per
  diagonal row. That `ambiguous_rows` remains unexercised is now stated in its own doc rather than left to
  be discovered.

  Four things the same review found and that are corrected here rather than shipped: a claim that
  `|volume_error| ≤ boundary_fraction` is an identity — **false**, and the 0.4 mm plate already in the
  suite is the counterexample at `+4.0` against `1.0`; a `[0.5, 1.5]` bracket whose lower half was wrong,
  because a column with one filled cell has one end and not two, so it is `[0.25, 1.5]`; a `1e-6`
  tolerance on the binary STL round-trip that measured **exactly zero**, because 30, 20 and 10 mm are all
  exact in `f32` and the effect it named was never in the configuration; and a comment justifying a
  convergence band with the *box's* boundary fraction while describing the sphere's.

- **A scene can declare a composite, not only a substance.** `Scene::composites`, `CompositeSpec`,
  `CompositePart`, `Palette::with_composites`. `crates/pantometry-world/tests/declared_composites.rs`,
  `scenes/22-wax-in-an-aluminium-matrix.json`.

  **The scene is the pair to read with `21`.** The same wax, now four fifths of a composite whose other
  fifth is aluminium, and `melted` climbs **125.8710 mm³/s against 100.6968** for the pure wax — exactly
  `1/0.8`, a ratio with the density and the latent heat cancelled out of it. Diluting the wax makes a cubic
  millimetre of buffer hold a fifth less latent heat, so the same twenty watts clear it a quarter faster.
  Machine precision on both: `6.0e-15` against the closed form and `1.7e-14` on the energy sum.

  The engineering reason for the metal is the other number: 5 W/m·K against the wax's 0.358, fourteen times
  better at moving heat into the thing that stores it. And the latent heat dilutes by **mass**, so four
  fifths of the volume is 54.67% of the mass and the composite stores 133.4 kJ/kg rather than 195.2 — the
  trap the scene exists to make visible.

  The gap between the two things 0.14.0 added. A scene could bring its own `Substance` and the library
  could mix two of them, and there was no way to say the second thing in a file — so a motor that is
  copper, steel, magnets and air still had to be one hand-computed material.

  `volume_fraction` is spelled out because volume and mass are the trap: volumetric heat capacity is
  volume-additive and specific heat is mass-weighted, and confusing them is worth 46% on a copper and FR-4
  board. Wax filling 80% of a volume is 54.7% of the mass.

  **The conductivity is the caller's and is checked.** No single value exists without the microstructure —
  the same two materials in the same proportions conduct 38 times more one way than the other — so the
  format cannot compute one. What it can do is refuse an impossible one, and the refusal is where the
  value is: a scene file has nowhere else to learn what range is achievable, so the message carries the
  Voigt and Reuss bounds *and* the tighter Hashin–Shtrikman pair. The emissivity is the caller's and cannot
  be checked at all, because a mixture has no surface.

  A composite may not be made of another composite. Not only an ordering problem, though a `BTreeMap` has
  no declaration order: nesting changes what the bounds *mean*, since Hashin–Shtrikman is a two-phase
  result and a three-part mixture has none. The refusal says to flatten it.

  Two things the tests caught, and the first is the one worth reading:

  **A material used only inside a composite was refused as dead weight** — the unused-declaration rule
  firing on the case it exists to protect. `Palette::with_composites` cleared its record of what had been
  asked for, on the reasoning that naming a part is not using the composite. True, and irrelevant: it also
  erased the parts. Three tests failed on one line, and the doc comment above it had described the correct
  behaviour while the code did the opposite.

  **A hand-typed bound was 1.2e-6 too high.** The test asserting that both bounds are accepted used
  `33.68644` against a true `33.6864`, so the assertion caught the test rather than the code. The endpoints
  come from `Mix` now, since what is under test there is the boundary condition and not the bound.

- **Every one of the eleven domains can now be asked a question from a file.** `structure`
  (`pantometry-elastic`), `channel` (`pantometry-fluid`), `cavity` (`pantometry-em`) and `well`
  (`pantometry-quantum`) join the scene format, and each ships a scene written around a closed form the
  domain's own documentation names.

  Measured rather than noticed: those four crates were referenced **zero** times in `pantometry-world`.
  They existed as libraries with their own tests and could not be touched from a scene file, the CLI
  or the browser — which is to say the platform could not ask them anything, while `scenes/README.md`
  had been saying "seven of the library's eleven domains" the whole time.

  Each scene is against the *discrete* answer where one exists, because the difference is a closed
  form and not an error term. A channel's Poiseuille mean is `(gh²/12ν)(1 + 2/n²)`, the `2/n²` coming
  from a no-slip wall imposed by reflecting the first cell — a parabola is not its own linear
  interpolation. A well's third level is `(2ℏ²/m dx²)sin²(nπ/2(N+1))`, below the continuum by
  `θ²/3` — 1.8320e-4 measured against 1.8322e-4 predicted. A cavity's mode comes out **under** the
  continuum frequency because a wave on a Yee grid travels slow, and the bound is that dispersion,
  `(kΔ)²/24`, rather than a tolerance.

- **A temperature field can drive a stress.** `structure` takes `follows`, naming a `block` whose
  temperature becomes its **stress-free strain**, `α(T − T_ref)` element by element.

  The abstraction is an **eigenstrain** and not a temperature: `Block::stress_free_strain` takes a
  dimensionless strain, so swelling, curing shrinkage and a phase change are the same statement and
  `pantometry-elastic` never has to depend on whatever computed it. Four closed forms check it, each
  blind to a different mistake, and all four exact at any mesh because a uniform eigenstrain makes a
  linear displacement field.

  `reference_c` is required rather than defaulted, and that is the whole design in one field: a power
  module is assembled at its solder's reflow temperature and *sits* at room temperature, so it is
  already strained before it is switched on. Scene 25 measures it — the module is most stressed
  **cold**, 0.2314 J of strain energy relaxing to 0.0274 J as it warms back towards where it was
  built.

- **Heat that has a place.** `dissipation` on a block: boxes of cells and their watts. Every other
  source in this format hands watts to the bus, and the bus carries an amount and no location — which
  is right for what the bus is and wrong for a die, a winding, a brake disc or a laser absorber.

  The watts are the box's **total, not a figure per cell**, so the answer does not move when the grid
  does: measured, a power module's junction moves by 2e-5 °C when every grid doubles. Scene 24 checks
  it against the resistance stack, 181.19 °C against 181.21.

- **A grid can have nothing in it, and a clearance radiates across it.** `Solid3D::empty` and
  `"material": "void"`: cells with no capacity, no conduction, no share of the bus, no vote in an
  average and **no temperature**. Two solid cells facing each other across the void exchange the
  parallel-plate series.

  `GapPatch` and `gap_patches` group a clearance into the sheets a view factor is a statement about,
  and `GapPatch::view_factor` gives the exact parallel-rectangle form. What it reports is **how much
  of the answer is a boundary condition**: the block charges the infinite-plate exchange, exact for
  the mirrored sides it has, and a gap open to space would see 0.415 of that for a 32 mm part 16 mm
  under a lid.

- **A block can lose heat.** `cooling` on a face, with an ambient, a film coefficient and an area, so
  a three-dimensional thermal scene can reach a steady state rather than warming for as long as it
  runs.

- **An assembly from CAD.** `parts` fills a block from STL files, one material each, with everything
  outside them void — which is what an assembly in air is. `Voxels::onto` puts two parts on one grid.

- **Where a part's bytes come from is no longer a path.** A `Parts` trait, with `OnDisk` and
  `Uploaded`: `World::build` reads from a disk, `World::build_with` from whatever it is given. The
  same STL voxelises to the same block **cell by cell** from either source, which is what makes a
  scene that runs in a browser the same scene as one that runs in a terminal.

- **Something chooses the cell size, by measuring.** `pantometry_world::fit` rasterises an assembly at
  every candidate cell and reports what each cost — filled cells, volume error, boundary fraction,
  thin runs, features below the cell, undecidable rows. The ladder is a statement about the assembly
  rather than about millimetres: the thinnest dimension of any part gets 1, 2, 4, 8 … cells.

  Predicting would have been wrong in a documented way. `Loss::volume_error` is a lattice-point count
  after the bulges and the cuts cancel, and a sphere at 2.5, 2.0 and 1.5 mm gives `+4.9%, +5.8%,
  −2.3%` — so the rule steers by boundary fraction, which does fall monotonically.

- **The platform's third verb.** `pantometry-world verify` measures what a passing audit does not:
  determinism to the byte, conservation margins, stability margins, and what moves when the coupling
  window halves or every grid refines.

- **An editor, and it runs in a browser.** `runtime/editor` is a fourth workspace: the scene's JSON
  checked as you type beside a wireframe, runs streaming in frame by frame, and the same thing as a
  wasm module with no backend — the whole library already compiled for one. CAD is dropped on the
  window, each file gets a material from the library's own catalogue, and **assemble** writes the
  scene.

- **An eleventh domain.** `pantometry-quantum`: a wavefunction in a well, marched with Visscher's
  staggered scheme, with probability an identity of the update rather than an accuracy claim.

- **A temperature has a colour, and it is computed.** Planck's law through the CIE 1931 observer to
  sRGB, so a glowing body is the colour it would be rather than a ramp somebody picked.

- **A scene can say where things are.** `poses` places a domain in the world, and `environment` makes
  the stage a statement rather than an assumption — gravity had lived in a constructor and no file
  could say otherwise.

### Changed

- **The kernel's exchange guard counts amounts rather than takes.** It refused two consumers of a
  channel by counting `take` calls, which an empty channel also triggers.

- **Three domains now count what was given to them from outside**, and each ledger adds its parts
  separately rather than their sum. `Solid3D::supplied`, `Block::received`, `Channel::driven`: a
  domain holding a conserved quantity that something outside the simulation added to is not a closed
  system, and `Ledger::add` raises an entry's *scale* to the largest thing added to it — so
  pre-summing a near-zero net leaves the first joule of rounding a hundred-percent error.

- **`Cavity::ledger` reports what leapfrog conserves**, `½ε|Eⁿ|² + ½μHⁿ⁻¹ᐟ²·Hⁿ⁺¹ᐟ²`, not the field
  energy. The naive sum swings by `2 sin(ωΔt/2)` about it — 7.4% for a 1.77 GHz mode — and pointing a
  1e-6 audit at it stops every correct run.

- **Seven domains gained `Domain::as_any_mut`.** It is opt-in with a silent default, so a coupling
  that wrote into one of them did nothing and reported nothing. `World::build` now probes a coupling
  and refuses rather than letting a future domain fail the same way.

- **`ScalarField` tells the truth about void.** `at`, `gradient`, `laplacian` and `rate` read the raw
  cell array, and an emptied cell still holds whatever it held when it was emptied — so a clearance
  left this workspace as a piece of the block sitting at its start temperature forever. Measured: a
  glTF carried 252 points for a grid with 120 solid cells.

- **The view layer draws only what is there.** A non-finite sample is no point in a glTF, no square in
  an SVG, and `null` in JSON — which is a token the format has, unlike `NaN`.

- **CSV keeps significant figures.** Both values and timestamps are `{:.9e}`. A fixed format is only
  readable at the scale it was chosen for, and a cavity holding 3.2e-10 J had written its entire
  energy history as a column of zeros.

### Fixed

- **A ray through a face's diagonal was counted twice, and parity approved it.** Found while writing the
  integration test above, and silent in the worst way. `Voxels::of` fills between sorted *pairs* of
  crossings; a ray through the edge two triangles share hits **both**, so an 8 mm cube's crossings come out
  `[0, 0, 8, 8]` — an even count, paired as `(0,0)` and `(8,8)`, filling **nothing**. The row came back
  empty, `ambiguous_rows` stayed at zero, and the loss report said the rasterisation was clean.

  A cube loses a whole diagonal plane that way: **448 cells of 512** at 1 mm, a slab straight through the
  middle, at every resolution. Every box test in the crate passed with the bug in it — a 30×20×10 mm box's
  diagonal is `z = y/2`, which needs an even `y` centre, and a cell-centred grid has none. Pure luck.

  Parity cannot detect this, so a hit on an edge is now a case of its own and the row is retried on a moved
  ray whatever its parity says. The regression test is a cube, at five sizes, asserted on the cell count and
  on the diagonal plane being solid; with the fix reverted it reports 448 against 512.

- **`Mesh::is_closed` reported a watertight mesh as open when a coordinate was negative zero.** The edge
  match is on bit patterns, deliberately — a gap of one bit passes a ray as readily as a gap of a
  millimetre — but `-0.0` and `0.0` are the *same point*, at no distance from each other, and their bit
  patterns differ. It fires on anything symmetric about an axis, where one side's coordinate is a product
  that happened to carry a minus sign. Negative zero is now folded onto zero, which is not a tolerance.

- **A bar and a lump can be made of something other than aluminium from Python.** `add_bar` and
  `add_lump` take a `material`, defaulting to aluminium so existing callers are unchanged, and
  `aluminium_heat_capacity_j_per_k` is now `heat_capacity_j_per_k(material, volume_m3)`.

  Both hardcoded `aluminium_6061`, and `add_lump`'s docstring did not say so — the worse half of a
  hardcoded material, because a caller modelling a copper busbar got an answer for a different metal with
  nothing anywhere saying which. The one-off capacity helper had to generalise with them: a test of a
  copper bar would otherwise have to compare it against aluminium's capacity or hardcode copper's, and one
  of those is wrong while the other is the constant that function exists to stop people copying.

- **Four of the nine catalogue materials could not be named from Python.** `bindings/python/src/lib.rs`,
  `bindings/python/tests/test_pantometry.py`.

  The binding held its own five-arm match — `copper`, `aluminium`, `electrical_steel`, `fr4`, `pla` —
  against a catalogue of nine, so `borosilicate`, `ice`, `stainless_304` and `water` were unreachable and
  nothing said so. The **third** copy of the catalogue's spelling in this workspace, and the third to go
  stale.

  `pantometry-world` had the identical defect and it cost eleven releases. The shape is the same every time: a
  name absent from a lookup is not a wrong answer, it is a substance that never appears, so nothing fails
  and no error is ever raised. `Substance::from_name` exists to be the one place the spelling lives, and
  the binding reads it there now — which also means the refusal lists whatever the catalogue holds rather
  than five names typed out beside it.

- **Five error messages read with a gap in the middle.** A `\` line continuation had been lost from each,
  leaving a run of eighteen spaces mid-sentence in what a caller sees — an unknown finish, a `tracks`
  naming nothing, an unknown tolerance channel, two domains with one name, and an unused material.

  Found by grepping for runs of spaces inside string literals, and the first fix was worse than the bug: a
  regex over the whole file collapsed the indentation of two JSON examples in doc comments, because the
  keys in them are quoted and the pattern could not tell a doc line from a literal. Reverted and redone
  with comment lines excluded outright, which is the check that should have been in the pattern from the
  start.

- **A cooled boundary was first order twice over.** The film was applied as a pass after the
  conduction sweep, which is Lie splitting, and its error carries a coefficient growing as `1/dx`
  against a step falling as `dx²`. Measured before the fix: ratios 1.27, 1.71, 1.87 per grid doubling
  rather than four.

- **A localised source made the steady state depend on the timestep.** It was applied before the
  sweep took its `old` snapshot, so the stencil conducted a share of the generated joules away inside
  the same step — a share equal to `G·dt/C`, about a half at the stability limit.

- **A cavity's mode was seeded at the wrong step.** `release_mode` staggers `H` by half of whatever
  `dt` it is given, and the builder gave it the Courant limit while the scheduler runs at
  `window / ceil(window / limit)`.

- **Two CAD parts could not touch**, because each rasterisation got a grid of its own.

- **A missing file reported what happened and not what it happened to.**

## [0.14.0] — 2026-08-13

### Added

- **Hashin–Shtrikman bounds for stiffness, and a witness that says how tight they are.**
  `Mix::bulk_bounds`, `Mix::shear_hashin_shtrikman`, `Mix::bulk_hashin_shtrikman`.
  `crates/pantometry-elastic/tests/a_checkerboard.rs`, and two tests added to
  `crates/pantometry-core/tests/a_mixture.rs`.

  The tighter elastic pair, assuming the microstructure is statistically isotropic. For aluminium against
  PLA at half and half it takes the shear range from 5.545-fold to 2.821-fold and the bulk range from
  4.57 to 2.37.

  **Checked against Mori–Tanaka**, which is the same equivalence the conductivity pair has with
  Maxwell–Garnett: HS with the matrix as reference *is* the Mori–Tanaka estimate for spherical inclusions,
  a separately derived result written as a different rational function, and the two agree to `2.2e-16`
  across two decades of inclusion fraction. Several plausible wrong versions of the HS expressions are also
  rational functions with the right limits at zero and one, so a dilute check would not have told them
  apart.

  Three things found, and two of them were mine:

  **The textbook prescription for which phase is the reference does not always mean anything.** "Stiffest
  for the upper bound" needs a well-ordered pair, and aluminium against borosilicate is not one — aluminium
  has the larger bulk modulus, 67.5 GPa against 46.5, and the *smaller* shear modulus, 25.9 against 34.0. A
  first version tested which phase had the larger value of the modulus being bounded, and for that pair at
  a tenth aluminium it returned a lower bound of 48.2312 GPa above an upper bound of 48.1922 — **the pair
  inverted**, by 0.08%, which only a sweep over fractions and pairs would see. Both evaluations are
  computed and then ordered now, which is what "interchange which phase is subscripted one" actually
  prescribes.

  **An affine boundary displacement bounds the effective modulus, not the bound.** The docs asserted that a
  kinematic estimate could not fall below HS+, reasoning that an upper estimate cannot cross an upper
  bound. It does not follow — the apparent stiffness bounds the *effective* stiffness, which is itself
  below HS+ — and the bulk modulus is the counterexample that fired the assertion. What is a theorem is
  `HS− ≤ C_effective ≤ C_apparent ≤ Voigt`, and only that is asserted.

  Read correctly, that makes the measurement a **ceiling on the truth**, and the two moduli answer
  differently: for shear the ceiling lands 0.505% above HS+, so the bound is tight to half a percent at
  worst; for bulk it lands 2.83% *below* HS+, so the bound is at least that loose. **The same pair is
  tight for one modulus of one geometry and loose for the other**, which the algebra does not say and a
  caller picking a number inside it would want to know.

  **A checkerboard one cell per phase is aliased to Voigt exactly**, arriving by a different route than the
  thermal one. There the cause was every face being an interface; here it is kinematic — with an affine
  boundary and a microstructure at the element scale a trilinear element has no freedom to relax into, so
  the affine field is the discrete solution and its energy is the volume average. Asserted, as it is in
  `a_composite.rs`.

  Two convergence parameters, separated deliberately: the period count is the fast knob and moves the
  answer 2.7% before plateauing, and cells-per-block is the slow one that moves *where* it plateaus — two
  cells sit 17.5% above HS+ and four get to 0.5%. A sweep changing both at once would report their sum as
  one rate, which is the mistake `a_layered_wave.rs` records making.

- **`Block::fill`, so the static solver takes per-element material too — and gets the same modulus nine
  orders sharper.** `Block::fill`, `Block::materials`, `Block::material_at`.
  `crates/pantometry-elastic/tests/a_layered_block.rs`.

  `Waves` got per-element material first because the closed form that checks a composite's stiffness —
  Backus averaging — is about wave speeds, and there was somewhere to point it. That left one crate with
  per-element material in the wave solver and not the static one, which is an asymmetry with no reason
  behind it.

  Statics turns out to give the **sharper** measurement of the overlapping claim. A traction-driven column
  is an elliptic solve with no time in it — no dispersion relation, no period to fit, no second-order mesh
  error — so a laminate's harmonic constrained modulus comes out to solver tolerance:

  ```text
                                 worst error      debug cost
    a_layered_wave.rs   C33        3.5e-4            67.5 s
    a_layered_block.rs             4.8e-13            0.13 s
  ```

  Nine orders sharper and 520× cheaper. The error there is the conjugate gradient's tolerance accumulating
  over degrees of freedom — `5.9e-15` at eight elements against `4.8e-13` at sixty-four — and not a mesh
  error, which would not care how many iterations it took.

  **The layer thickness does not change the static answer**, and that contrast is a statement neither file
  could make alone. The wave measurement of the *arithmetic* mean is 23% out at eight layers per
  wavelength, because it needs the layers to move together and only a long wave makes them. A static
  uniform stress needs nothing of the kind — equilibrium makes the stress uniform whatever the layers look
  like — so the same column at one, two, four and eight elements per layer gives the same modulus to
  `4.8e-13` with no trend. Backus averaging is a limit for one of those quantities and an identity for the
  other.

  A `fill` on a solved block marks it **unconverged**. The displacement on record solves the previous
  assembly, and left alone it would still read as converged and still return a strain — an answer that is
  present and wrong rather than absent, which is the failure this workspace keeps finding.

  What statics cannot do is the arithmetic end, for the same reason the wave file cannot do `C11`:
  realising it needs the lateral strain zero on average while free locally, and `roller` constrains a face
  rather than a plane through the interior. Statics gets the sharp half; the wave gets both halves less
  sharply. Neither covers it alone.

- **Per-element material in `pantometry-elastic`, and stiffness bounds that something can check.**
  `Waves::fill`, `Waves::materials`, `Waves::material_at`, `Elastic: PartialEq`, `Mix::shear_bounds`,
  `Mix::p_wave_modulus_bounds`. `crates/pantometry-elastic/tests/a_layered_wave.rs`.

  `Mix` shipped without stiffness bounds a release ago, and the reason was stated at the time: `Waves`
  took one material per block, so there was nothing here a bound on stiffness could be *checked* against,
  and a bound nothing can falsify is a comment rather than an API. This closes that.

  The closed form is **Backus averaging** — the exact long-wavelength moduli of a layered elastic medium,
  Backus 1962 — and `Waves::hold` is what makes each one separately measurable, because freezing two
  displacement components leaves a one-dimensional problem with one modulus in it. Aluminium and PLA, half
  and half, a 20-fold contrast in shear modulus:

  ```text
                                            16/32/64 elements        improvement
    C44  harmonic mu   2.4517 GPa      0.558% / 0.143% / 0.036%        15.7x
    C66  arithmetic mu 13.5945 GPa     1.003% / 0.229% / 0.058%        17.2x
    C33  harmonic M    11.1237 GPa     0.533% / 0.143% / 0.035%        15.4x
    M    arithmetic    53.9839 GPa     0.817% / 0.189% / 0.048%        17.1x   (32/64/128)
  ```

  Second order throughout, where a fourfold refinement predicts 16. **Both ends of the shear pair come
  from one block**, 5.5× apart, decided by which way the shear goes — and the *ratio* of the two speeds
  lands within 0.023% of `√5.545`, tighter than either measurement because the mesh error cancels out of
  it.

  **A composite of two isotropic materials is generally anisotropic**, so `Mix` deliberately reports no
  effective `(E, ν)` and `as_substance` leaves the mechanical block absent. A laminate 5.5× stiffer in
  shear one way than the other has no single isotropic pair, and a function returning one would be
  inventing an isotropy the material does not have.

  Two things the measurements corrected:

  I documented the Voigt end of the P-wave modulus pair as a bound that is **not** attained. It is
  attained — holding the lateral strain at zero *pointwise* is exactly the equal-strain condition Voigt
  assumes, and it converges to `⟨M⟩` at second order. What is not measurable here is `C11`, the *free*
  laminate compressed along its layers, which needs the lateral strain zero on average but free locally;
  `Waves::hold` holds a component on every node or none, so this API cannot pose it. Tried anyway and got
  40.87 GPa, which is neither `⟨M⟩` nor `C11`'s 43.77 — a block four elements thick with free faces
  carries plate modes, so that number answers a third question. No claim is made about it.

  And the arithmetic mean is a **long-wavelength** result in a way the harmonic ones are not, which the
  resolution sweep alone cannot see: refining the span shrinks the element *and* stretches the wavelength
  against the layer at the same time, two second-order errors falling together and reported as one. Held
  the element size fixed and varied only the layer thickness: 0.234%, 1.055%, 4.877%, 23.296% at 64, 32,
  16 and 8 layers per wavelength. **At eight layers per wavelength a caller is 23% out, and that is the
  closed form's limit rather than the solver's.**

  The file costs 68 s in debug, and that is written into it, because the first draft cost 5.4× more for
  the same claims.


- **Composites: two substances made into one, with the bounded properties returned as bounds.**
  `pantometry_core::mixture::Mix`, re-exported as `pantometry::prelude::Mix`.
  `crates/pantometry-core/tests/a_mixture.rs`, `crates/pantometry-thermal/tests/a_composite.rs`.

  A motor is copper, steel, magnets and air; a board is FR-4, copper and solder; a buffer is wax in an
  aluminium matrix. Each wants to be **one** `Substance` so a lumped model can hold it, and
  `with_specific_heat` has been saying since 0.4.0 that the bulk `c_p` of such a thing is worth a factor
  of two — while leaving the caller to compute it by hand.

  The properties divide into three kinds and conflating them is the whole failure mode. Density,
  volumetric heat capacity and latent heat are **exact**, from conservation alone. Conductivity is
  **bounded** — Voigt and Reuss for any microstructure, Hashin–Shtrikman if it is isotropic. Emissivity is
  **nothing**: it is a property of the surface and a mixture has no surface, so `as_substance` takes it as
  an argument. A conductivity outside the bounds is **refused**, because no microstructure realises it.

  **One block gives both bounds, to machine precision.** A laminate of alternating aluminium and
  borosilicate measures Reuss `2.213236` across its layers and Voigt `84.057000` along them — same
  material, same fractions, a factor of 38 apart, decided by the direction of the flux. That is why `Mix`
  returns a pair and refuses to pick: a single number would be the upper bound of a 38-fold range
  presented as a measurement. On a half-copper, half-FR-4 board the range is 335-fold.

  Checked against a resolved geometry rather than against itself, and against Maxwell–Garnett for the
  algebra: the Hashin–Shtrikman lower bound with the matrix as host **is** Maxwell–Garnett, to `6.7e-16`
  across three decades of filler fraction.

  Three things found, all of them mine:

  A three-dimensional checkerboard one cell per phase is **not a checkerboard**. Every cell then has all
  six neighbours of the other material, so every face carries the same harmonic mean and the discrete
  operator is that of a *uniform* medium at `harmonic(167, 1.114) = 2.213236` — which at equal fractions
  is Reuss to the last digit. Two drafts reproduced the lower bound and read it as a striking result; the
  microstructure was aliased, not resolved. Asserted now, so it is not rediscovered a third time.

  A coarsely resolved high-contrast composite is under-conductive by enough to **break a bound**: at four
  blocks per axis the measured conductivity goes 3.4926, 5.1396, 6.8213 as the grid is refined, and the
  first of those is below the Hashin–Shtrikman lower bound of 4.3266. Not because a checkerboard violates
  HS, but because two cells cannot represent one. The same pattern in stainless against aluminium does the
  same thing, which is what says it is the discretisation and not the pair.

  A steady-state stopping rule of "the flux stopped changing over N steps" is **not resolution
  independent**: as `dx` shrinks, N steps is a shorter physical time, so the criterion becomes trivially
  true. It reported a 32³ block as converged after 7.5 microseconds of a 5-second process. The rule is an
  absolute residual now — flux in equals flux out, because at steady state nothing is stored — and every
  tolerance in that file traces to a named constant rather than to what made the run pass.

  Elastic and acoustic bounds are deliberately absent: Voigt–Reuss on the bulk and shear moduli is the
  same theorem, but `pantometry-elastic` has no per-cell material, so there is nothing here a stiffness bound
  could be checked against. Yield strength is absent for a different reason — a composite's yield is
  governed by the weaker phase and the interface, so a rule of mixtures for it would be wrong rather than
  imprecise.

- **A scene can declare its own substances, so the catalogue stops being the limit.** `Scene.materials`,
  `Palette`, `Substance::CATALOGUE` and `Substance::from_name`.
  `crates/pantometry/tests/substances_from_a_file.rs`, `crates/pantometry-world/tests/declared_materials.rs`,
  `scenes/21-a-wax-thermal-buffer.json`.

  Nine catalogue entries against hundreds of thousands of materials: enumeration was never going to
  close that gap, and a format that can only *name* a substance can only describe nine kinds of thing.
  A `Substance` is data, so a scene writes one out and uses it exactly as it uses `"ice"`.

  **Checked against Neumann's exact solution, with ice marched beside it through the identical
  harness** — which is the only form in which "as accurate as one we ship" means anything. Gallium and
  n-octadecane, 82× apart in diffusivity, declared as JSON text: worst 0.039% over a 21× range of
  Stefan number, against the catalogue's own ice at 0.035%. Ice sits *inside* the declared range at
  every undercooling, so the error is a function of Stefan number and not of where the substance came
  from.

  Three things found on the way, two of them mistakes of mine:

  `MATERIALS` was a hand-written copy of the catalogue's spelling — eight names beside nine
  constructors — and the missing one was `water`. Unnameable from a scene for **eleven releases**, since
  the format learned to name a material at all. A name absent from a lookup is not a wrong answer, it is
  a substance that never appears, so nothing could have noticed. The lookup lives beside the catalogue
  now and a test checks both directions.

  The first tolerance was `dx/X` — one cell of interface resolution, 8.3%. Measured error: 0.039%,
  **two hundred times smaller**, because the front is read from an integral of a conserved quantity and
  the enthalpy scheme conserves that exactly, so the interface's first-order error never reaches the
  answer. A bound can trace to real physics and still be so loose it checks nothing: that one would have
  passed a conductivity 13% wrong. At 0.06% it catches 0.1%, and that is asserted rather than claimed —
  a test perturbs the conductivity and measures what the front does.

  The teeth test then failed to have teeth. It asked for the same front *depth* from the perturbed
  substance, which recomputes its own target from the wrong conductivity, so the two errors cancelled to
  the last digit and a 1% error reported as 0.0000%. The comment above it had named that exact trap.

  A declaration nothing uses is **refused**, and that is the non-obvious one: `material` on a block is
  optional and defaults to aluminium, so declaring a wax and forgetting to name it gives a block of metal
  that runs, audits and renders while answering about the wrong substance. Same shape as a region
  selecting no cells, one level up.

- **Two-phase conduction, against the two-phase Neumann solution.** `FusionProps::liquid`, `new` and
  `with_liquid`. `crates/pantometry-thermal/tests/two_phase_stefan.rs`.

  The one-phase model is exact while the liquid sits at the melting point, because a face with no
  temperature difference carries no heat whatever its conductivity. It stops being right the moment the
  liquid is warmer, and not by a little: **20 K of superheat slows a freezing front 16%**, from 15.85 mm
  to 13.33 mm at 900 s.

  The closed form was established before any code changed, and its check is the **reduction**: setting
  the superheat to zero must give back the one-phase condition. It does, proportionally — the gap is
  `9.446e-12`, `9.447e-9`, `9.446e-6` at `1e-9`, `1e-6` and `1e-3` kelvin, which is `9.4465e-3` per
  kelvin over three decades. Three decades of exactly proportional error is a stronger statement than
  one equality.

  ```text
    superheat   front      two-phase says     one-phase says
        0 K   15.50 mm   15.85 ( 2.20%)     15.85 ( 2.20%)
        5 K   14.70      15.13 ( 2.89%)     15.85 ( 7.27%)
       20 K   12.83      13.33 ( 3.78%)     15.85 (19.07%)
  ```

  Worst 3.78% against a `dx/X` bound of 6.31%, and at 20 K the one-phase model is **19.1% out** — three
  times the bound. `dx/X` and not its square, because the cell holding the interface has a *mixed*
  conductivity and that is a first-order error.

  **Three things I had wrong, all caught by measurement:**

  - **The kelvin-normalised latent heat breaks with two capacities.** `L/c_p` multiplied by a cell's
    *current* capacity charged `4182/2050` = 2.04× to freeze water, and the front came out **27%
    short** with nothing pointing at the latent heat. The enthalpy map is in joules now, with both
    capacities named.
  - **"A liquid at the melting point cannot influence anything" is false for the scheme.** True of the
    continuum; a fixed grid puts a *wholly* liquid cell where the continuum has an interface, and its
    low conductivity throttles the front — 2.1% at twenty cells, closing 2.99× over a fourfold
    refinement. So `Substance::ice` keeps `liquid: None`: switching it on took the *one-phase* answer
    from 0.43% to 6.9%, sixteen times worse for physics that had not changed, and a default that costs
    a caller a factor of sixteen is the wrong default. Two phases are opt-in.
  - **Melting loosens the stability limit and freezing tightens it.** I checked one direction. A step
    sized on an all-liquid block is refused a few hundred steps later once there is ice in it, water's
    limit being 7.57× ice's. The guard caught it, which is how it was found.

  Cost measured before designing: a `resolve` is 1.6× a `step` at 40 cells and 4.2× at 4096, so a
  two-phase sweep runs 2.6× to 5.2× a one-phase one. Simple and whole-grid; an incremental update
  touching only the mush is the optimisation available if a problem needs it, and none does.

  And `FusionProps` gaining a field broke every struct literal — exactly what `Substance` gaining
  `fusion` did one level up, to tests written the week before. `FusionProps::new` and `with_liquid` are
  the same answer applied one level down.

### Changed

- **A packed bed's conductivity is Maxwell–Eucken now, not an arithmetic mean.** `Puck::bed_conductivity`
  and `Puck::conductivity_at`. `crates/pantometry-porous/tests/the_beds_conductivity.rs`.

  It was `ε k_l + (1−ε) k_s`, which is the **Voigt bound** — exact only when the two phases lie in
  parallel with the flux, which in a bed of spheres nothing does. Nobody chose it; it is what you get if
  you do not think about it. The structural fact about a saturated bed is that the **liquid is the
  continuous phase** and the grains are dispersed in it, which is precisely what Maxwell–Eucken describes
  and is attained by a coated-sphere assemblage. For coffee at 45% porosity the old value was **11.0%
  high**: 0.38625 against 0.34811 W/m·K, and the honest range narrows from 1.674 to 1.184.

  Cross-checked rather than asserted: the puck's Maxwell–Eucken equals `Mix::hashin_shtrikman`'s bound
  built around the liquid to `2.1e-15` over 76 cases, and `Mix`'s form is itself checked against
  Maxwell–Garnett. Two independently written arrangements of one physics; the puck does not call `Mix`.

  **Why it went unexamined for four releases, measured both ways.** Under flow the bed is isothermal, so
  conduction carries no heat and `λ` is multiplied by zero — swinging it over a factor of eight leaves the
  extraction yield identical to `1e-14`, and scene `18` is bit-for-bit unchanged by this commit. It stops
  being free with a gradient: in a 20 °C basket the yield moves 4.9% per unit `ln λ`, so the range the old
  rule left open was worth about 2% in extraction. Both halves are now tests, so the choice cannot go back
  to being invisible.

  What moved: the cold-basket branch of `espresso_shot`, where the cup cools 3.0102 °C instead of 3.1829
  and the yield reads 19.43% instead of 19.44%. Both stay inside the capacity bound the example asserts,
  and no assertion in the workspace changed — every gradient claim here is an inequality or a structural
  statement rather than a recorded value, which is the only reason a model change was safe to make in one
  commit.

  One thing I chased that was not there: a 287% yield, above the 30% soluble ceiling. The domain saturates
  at exactly 30.00000% when driven to it; a percentage in the printout had been multiplied by a hundred
  twice.

### Fixed

- **`CITATION.cff`'s licence field, which failed the 0.13.0 Zenodo deposition.** `license:
  MIT OR Apache-2.0` is a valid SPDX *expression* and `Cargo.toml` is right to use it; CFF's schema
  takes an identifier or a **list** of them and an expression matches neither. Validated against the
  real schema rather than guessed, and now a two-element list.

  0.13.0 therefore has **no DOI** — it went to crates.io and PyPI and Zenodo rejected the metadata,
  reporting it only as a red *Failed* on its own web page. Nothing in the release, the tag or CI knew.
  0.14.0 will be the first version with one; re-depositing 0.13.0 would mean deleting and recreating
  its GitHub release, and it is not worth doing to a version that is already out.

  `crates/pantometry-world/tests/citation_is_valid.rs` is four checks so the next one fails in the gate
  instead of on a web page: the licence is a list, every field a deposition is built from is present,
  the version matches the crate version, and the author block is a list of mappings rather than a bare
  string — that last being the silent half of this class, since a rejection is loud and a record with
  no creator is not.


## [0.13.0] — 2026-08-12

**Elastic waves, any material, and one behaviour change that will bite a deserialiser.**

`pantometry-elastic` had a `density` field its own documentation called unused; it sets the wave speeds,
and `Waves` marches `ρü = ∇·σ` on the element `Block` solves statics with. That closed
`ARCHITECTURE.md`'s depth list — the last entry turned out not to be finite strain but **inertia**.

And `Substance` can describe a material this crate never chose. Four `with_*` builders and a `check`
that says whether a number is *possible*, because enumeration does not reach "every material" and data
does. `Substance::bulk` set all four property blocks to `None` and there was no way to fill any of them
except a struct literal — which 0.12.0 broke when it added `fusion`.

### Breaking

- **`Substance` and its four property blocks now refuse unknown keys.** `deny_unknown_fields`, so a
  JSON material with a mistyped `"thermalz"` is an error instead of a substance whose whole thermal
  block is silently absent. If you deserialise `Substance` from a file that carries extra keys, that
  file stops loading — which is the point: it was previously loading as a *different* material.

- **`Waves::new` takes `(name, counts, cell, material)`**, matching `Block::new`. It shipped in no
  release; noted because the two disagreed for a day on `main`.

### Added

- `pantometry_elastic::Waves` — `ρü = ∇·σ`, central differences, lumped mass. `hold`, `clamp_ends`,
  `release_mode`, `mode_amplitude`, `mode_frequency`, `kinetic_energy`, `strain_energy`,
  `total_energy`, `displacement_at`. Second order in both speeds — 0.161%, 0.040%, 0.010% over 16, 32
  and 64 elements — and the **ratio** holds to `1e-8` because `E` and `ρ` cancel algebraically while
  the mesh error cancels numerically.
- `pantometry_elastic::Axis`, replacing five `usize` axis arguments that disagreed about a fourth one.
- `Elastic::p_wave_speed`, `s_wave_speed`, `speed_ratio`, `from_substance`.
- `Block::as_field` and `Waves::as_field` — the displacement magnitude. Neither half of this crate
  could be drawn before, for as long as the crate has existed.
- `Face::on`, returning an `Axis` where `Face::axis` returns an index.
- `Substance::with_thermal`, `with_mechanical`, `with_acoustic`, `with_fusion`, `check`.

### Fixed

- `bindings/python/Cargo.lock` was stale at 0.10.0 through the whole of 0.12.0. Nothing in that job
  passes `--locked`, so cargo rewrote it silently on every build and the committed copy drifted.

### Documentation

- `RELEASING.md` is new: the release procedure was inside `CLAUDE.md`, which is loaded every session
  for a thing you do once per release. `CLAUDE.md` went from 246 lines to 115.
- `CONTRIBUTING.md` records **six** times the gate reported a pass it had not earned, one of which
  reached `main`, and the measured fact that `set -euo pipefail` does not prevent it in a pasted block.
- No `Co-Authored-By` trailer on commits from `05ed74e` onward.


### Added

- **Any material, not the nine in the catalogue.** `Substance::with_thermal`, `with_mechanical`,
  `with_acoustic`, `with_fusion` and `check`, plus `deny_unknown_fields` on `Substance` and all four
  property blocks. `crates/pantometry/tests/any_material.rs`.

  Enumeration does not reach "every material" and data does. `Substance::bulk` set all four property
  blocks to `None` and there was **no way to fill any of them** except a struct literal naming every
  field — and that path broke earlier in this release when `fusion` was added, at the expense of
  exactly the callers the catalogue is least able to help. A chain of `with_*` is immune to a new
  field appearing.

  `check` is what the library can do for a number it did not choose: not whether it is *right*, but
  whether it is **possible**, reporting every problem at once because a material read off the wrong
  column is usually wrong in several places. Plus one check that is not a bound on a single field — a
  substance stating both a sound speed and elastic constants has three independent numbers describing
  one thing, and the speed must sit within **15%** of the rod or bulk speed those give. Measured
  rather than chosen: every catalogue entry is within 6.2% of whichever it means, and that gap exists
  because a tensile test and an ultrasonic measurement are not the same measurement.

  What it catches: a shear speed in the longitudinal slot (39% out), a modulus from the row below
  (27%), an emissivity written as a percentage, a negative conductivity, an incompressible Poisson
  ratio, and four wrong at once reported as four.

  **The test written to demonstrate the capability found the capability was unsafe.** `Substance` had
  no `deny_unknown_fields`, so a mistyped `"thermalz"` was silently dropped and the material ran as
  one whose conductivity is *unknown* rather than one whose file has a typo. The scene format already
  had that rule written down with its reasoning; `Substance` has it now, on all five types.

  Ti-6Al-4V and paraffin wax run through conduction, elastic waves and a melting plateau against the
  same closed forms a catalogue material does — `dx²/6α` to `1e-12`, `√(2(1−ν)/(1−2ν))` at `ν = 0.342`
  which no catalogue entry has, and `ρLV/P` to within one delivery.

- **`Elastic::from_substance`, and "small strain" spans 130× across the catalogue.** The conversion
  from `pantometry-core`'s material catalogue to `pantometry-elastic`'s material existed **in a test** —
  `two_wave_speeds.rs` built one by hand — which means every consumer wanting to solve an elastic
  problem with a catalogue material wrote the same four lines. It is in the library now, in
  `pantometry-elastic`, because the kernel must not learn that elasticity exists.

  It **drops the yield strength**, which is the one thing to know: `Substance` says where a material
  stops coming back and this type has no yield, so a solve past it returns a displacement that is
  arithmetically correct and physically meaningless with nothing in the answer to say which. So the
  strain where the linear model ends is documented and measured, and the spread is the finding:

  ```text
    ice        0.011%     brittle: 1 MPa against 9.1 GPa
    Cu ETP     0.060%
    N-BK7      0.073%
    Al 6061    0.401%
    PLA        1.429%     a polymer has twenty times a metal's elastic room
  ```

  A first draft asserted every entry was under 1% and failed on PLA at 1.43% — correctly, because a
  polymer is not a metal. The bound was wrong, not the data, so the claim is now the **spread**: 130×,
  meaning a strain bound that would be absurdly conservative for PLA is already past yield for ice.

- **Neither half of `pantometry-elastic` could be drawn, and its two constructors disagreed.** A
  stabilisation pass, continued: `Block::as_field`, `Face::on`, and `Waves::new`'s arguments in
  `Block::new`'s order.

  `Block` offered no field either, and had not since the crate existed — so an elastic run drew nothing
  and nothing said so, for statics as well as waves. Both now offer `|u|`, trilinearly interpolated,
  checked against the displacement the solver reports rather than against a picture.

  `Waves::new` took `(name, material, counts, cell)` where `Block::new` takes
  `(name, counts, cell, material)`. One crate, two orders, and only a caller writing both in one file
  would ever have noticed. Aligned while `Waves` is unreleased.

  `Face::axis` returns an index and predates `Axis`; `Face::on` returns the type. Two spellings of one
  idea inside one crate drift, so a test pins all six faces to both.

- **`pantometry_elastic::Axis`, and `Waves` is drawable.** A stabilisation pass over what landed
  unreleased, which is the moment to change an API rather than after somebody depends on it.

  Five methods on `Waves` took an axis as a `usize` and **disagreed about a fourth one**: four returned
  silently — leaving a body with three free components where the caller believed it had one, and so a
  wave speed that is not the one being asked about — and `mode_frequency` clamped to 2, which is a
  different wrong answer to the same mistake. Three responses to one bad input in one type is not a
  policy, so the input is now unrepresentable.

  And `Waves` offers `as_field` — the displacement **magnitude**, trilinearly interpolated — so the
  analysis layer can draw it. Every other 3D field domain offers one, and a domain that offers none
  gets no picture from a layer whose whole rule is to dispatch on the shape of the data. `|u|` costs
  the sign, which is stated: a full standing wave shows both antinodes bright. `Block` still offers
  nothing here, which is a real gap and older than `Waves`; it is on a published type, so that is a
  separate decision.

  The per-step cost was **measured before deciding** anything about it: 2.47 ms/step at 16³ elements
  and 19.69 ms at 32³, linear in nodes, so the two `Vec` allocations a step makes are under 1% and are
  left alone. A wave crossing a 32 mm body is about 90 steps; a body left ringing for a millisecond is
  8,700 and three minutes.

- **Elastic waves: `pantometry_elastic::Waves`, and `ARCHITECTURE.md`'s depth list is closed.**
  `ρü = ∇·σ` on the same trilinear element `Block` solves statics with, marched with central
  differences and a lumped mass. `Waves::hold`, `clamp_ends`, `release_mode`, `mode_amplitude`,
  `mode_frequency`, `kinetic_energy`, `strain_energy`, `total_energy`, `displacement_at`.
  `crates/pantometry-elastic/tests/two_speeds_marched.rs`.

  A sibling type rather than a mode on `Block`, because `Block` is `Kind::QuasiStatic` and is right
  about itself — a body with velocity has a different lifecycle, a different stability limit and a
  different thing to conserve. What they share is the element, which is the part worth sharing: the
  static tests check that operator against four exact moduli and these check it against two exact
  speeds.

  ```text
    elements   c_p error   c_s error
          16     0.161%      0.161%
          32     0.040%      0.040%
          64     0.010%      0.010%
  ```

  Second order in both, a factor of four per halving to two figures. And the **ratio holds to about
  `1e-8`** at `ν` = 0.2, 0.33 and 0.45 — four orders tighter than either speed alone, because `E` and
  `ρ` cancel out of it algebraically *and* the two modes share a shape so the mesh error cancels out
  of it numerically. What is left being checked is the operator.

  The **leapfrog's own dispersion is removed rather than tolerated**: central differences on one
  eigenmode turn at `Ω = 2·arcsin(ω dt/2)` per step, so the measured `Ω` is inverted through
  `ω = 2·sin(Ω/2)/dt` before anything is compared. Comparing raw periods would add the time error to
  the space error and call the sum an accuracy.

  Energy swings by `2 sin(ωΔt/2)` — measured 0.0348 against 0.0348 predicted — and does **not** drift:
  the mean over the last eighth of a forty-period run matches the first to `2.5e-7`.

  The stability limit is `2/√λ_max(M⁻¹K)` by Gershgorin on the assembled rows, not a borrowed Courant
  number. Measured it is **1.229×** `dx/(c_p√3)`, so the Courant form would have been safe and 23%
  wasteful — the opposite direction from the one assumed before measuring it, which is why it is
  computed.

- **A solid has two wave speeds, and the catalogue's one field means both of them.**
  `Elastic::p_wave_speed`, `s_wave_speed` and `speed_ratio` in `pantometry-elastic`, and
  `crates/pantometry/tests/two_wave_speeds.rs` — the fourth cross-domain test of that shape, after
  `fields_and_rays`, `loss_and_lumps` and `a_slit`.

  `Elastic::density` was documented as "unused by the static solve and carried so a mass can be
  stated once". It sets the wave speeds, which is the only place a mass enters a problem with no
  inertia in it.

  The sharp statement is `c_p/c_s = √(2(1−ν)/(1−2ν))`: both `E` and `ρ` cancel, so it holds at
  **2.2e-16** for all six catalogue entries and a scheme with the wrong stiffness or the wrong mass
  still has to land on it.

  The finding is what `AcousticProps::sound_speed` turns out to mean. A fluid has one speed; a solid
  has two longitudinal ones — bulk `√((λ+2μ)/ρ)` and rod `√(E/ρ)` — and the six entries do not all
  mean the same one:

  ```text
                      rod   stated    bulk    vs bulk   vs rod
    ice              3150     3840    3834      +0.1%   +21.9%     bulk
    Al 6061          5052     6320    6149      +2.8%   +25.1%     bulk
    304 stainless    4912     5790    5623      +3.0%   +17.9%     bulk
    Cu ETP           3614     4760    4483      +6.2%   +31.7%     bulk
    N-BK7            5716     5680    6048      −6.1%    −0.6%     rod
    elec. steel      5113     5100    5853     −12.9%    −0.3%     rod
  ```

  Four bulk, two rod, every one within 6.2% of whichever it is. And read as a bulk wave the stated
  speed implies a modulus: 9.1 against 9.1 GPa for ice, but 72.8 against 68.9 for aluminium and 131.9
  against 117.0 for copper — because **a tensile test and an ultrasonic measurement are not the same
  measurement**, and the dynamic modulus of an annealed metal comes out higher. Ice agrees exactly
  because its elastic constants were back-calculated from velocity, so the one entry that agrees is
  the one where agreement was never independent.

  Two earlier drafts asserted `rod ≤ stated ≤ bulk` with 0.1% and then 2% of slack. Both were false —
  three entries sit above the bulk speed and two below the rod speed — and the second failed on
  aluminium, which is not an edge case. The assertion is now a classification: each entry is within 7%
  of one of the two, and that the split is *mixed* is itself asserted, so a future edit cannot make
  the point disappear by making them uniform.


## [0.12.0] — 2026-08-12

**Two of `ARCHITECTURE.md`'s three depth entries, and both the same shape.** Single-material and
single-phase are answered; small-strain is the one left. Neither cost a new crate, neither touched the
kernel or either layer above it, and each turned on a **closed form** rather than on a feature — the
harmonic mean making a layered wall's resistance exact, and Neumann's solution putting a freezing
front in a place rather than at a rate.

Also diffraction, which measures where scalar optics stops being true rather than assuming it is.

Twenty-two new public items across four crates. `Solid3D::max_stable_dt` reports a different number
for a thin or filled block than 0.11.0 did — read the entry below before pinning a step by hand.

### Added

- **Latent heat, against Neumann's exact solution of Stefan's problem.** `LatentHeat` in
  `pantometry-units`, `FusionProps` and `Substance::fusion`, `Substance::latent_energy` and
  `Substance::ice` in `pantometry-core`, `Solid3D::set_melted_fraction`, `melted_fraction_at` and
  `melted_volume` in `pantometry-thermal`, `"ice"` in the scene format's materials.
  `crates/pantometry-thermal/tests/a_freezing_front.rs` and `20-melting-a-block-of-ice`.

  The second of `ARCHITECTURE.md`'s three depth entries — every domain was single-phase — and it left
  one: small-strain. No new crate, and nothing outside `pantometry-thermal` but the material data.

  **Why this problem:** it has an exact solution, which almost nothing with a moving boundary does. A
  semi-infinite liquid at its melting point, surface dropped to `T_s`, freezes to `X(t) = 2λ√(αt)`
  with `λ e^{λ²} erf(λ) = St/√π`. A *position*, not a rate and not a limit. Measured for ice under 20 K
  of undercooling:

  ```text
    dx = 1.0 mm    t = 100 s    5.311 mm against 5.288 mm     0.431%
                   t = 400 s   10.570 mm against 10.567 mm    0.032%
                   t = 900 s   15.843 mm against 15.848 mm    0.033%
    dx = 0.5 mm    t = 100 s    5.285 mm against 5.283 mm     0.032%
                   t = 900 s   15.847 mm against 15.848 mm    0.012%
  ```

  And nine times the time gives **2.9984×** the depth against `√9 = 3` — the `√t` law with no closed
  form in it at all, which a scheme with the right coefficient and the wrong power would fail.

  **The scheme is enthalpy**, bookkept as a temperature and a melted fraction: a cell's state is one
  monotone number — `T − T_m` below, `φ·ℓ` inside, `ℓ + T − T_m` above — so the sweep adds energy and
  inverts it. Nothing tracks the front, nothing iterates, and energy is conserved because energy *is*
  the state. Not **apparent heat capacity**, which smears `L` over a temperature interval and lets a
  step big enough to cross the interval skip the latent heat silently, running the front fast. Here an
  overshoot's remainder lands on the far side because the inverse map says where that much energy
  goes — verified at machine precision with ten times the latent heat in one delivery.

  A cell fed constant power holds at its melting point for exactly `ρLV/P`: measured 305.9 s against
  305.9 for a cubic millimetre of ice at a milliwatt. For ice `L/c_p` is **163 K**, so a scheme that
  dropped the latent heat would not be slightly wrong.

  Four things the tests found, all of them mine:

  - **The profile's convergence order is not clean** between adjacent resolutions — 0.78 then 1.84,
    and averaging five instants did not settle it. A fixed grid makes the front advance in a
    staircase, so the field at any moment depends on where the front sits between two cell centres and
    that phase is not a smooth function of `dx`. Over the full fourfold refinement it is 6.2×, between
    first and second order. Two earlier drafts asserted second order on no evidence and then first
    order on evidence that did not earn it.
  - **"Every cell behind the front" is the wrong metric.** It measures a different set of points at
    each resolution — the nearest one to the front moves inward as `dx` shrinks, into exactly the
    region a cell-centred profile represents worst. That alone read as order 0.66. Fixed depths.
  - **The front position is more accurate than the field it comes from** (0.03% against 0.8%), because
    `melted_volume` is an integral of a conserved quantity where a temperature is a point sample.
  - **A checkpoint of temperatures alone is not a checkpoint.** 0 °C is ice, water or any mixture, so
    `checkpoint`/`restore` carry the phase — a live path, since `Schedule::Iterative` restores every
    iteration and the conservation audit restores on a violation. Without it a half-frozen column
    comes back as a fully liquid one at the same temperatures, and the ledger balances, because the
    ledger reads the state that was corrupted.

- **A block can be made of more than one material.** `Solid3D::fill`, `substance_at`, `substances`,
  `face_conductance`, `heat_capacity` and `stability_ratio` in `pantometry-thermal`, and `regions` plus
  `material` on the scene format's `block`. `crates/pantometry-thermal/tests/a_layered_wall.rs` and
  `19-a-coating-stops-the-heat`.

  The first of `ARCHITECTURE.md`'s three *depth* entries — single-phase, single-material,
  small-strain — and it cost no new crate and nothing outside `pantometry-thermal`.

  The number it turns on is the conductivity **on a face**, which is the harmonic mean
  `2k_Lk_R/(k_L+k_R)` and not the arithmetic one. For aluminium against borosilicate those are
  2.21 and 84.1 W/m/K, a factor of **38**, and the arithmetic mean short-circuits the interface.
  The harmonic mean earns an *equality*: with the material interface on a cell face the discrete
  chain of face resistances is exactly `Σ Lᵢ/(kᵢA)` at every resolution — measured at exactly zero
  for twelve cells and `2.2e-16` for twenty-four, forty-eight and ninety-six. A layered wall's
  resistance has no discretisation error at all.

  The arithmetic mean's is **first order**, which is the dangerous kind — 8.852%, 4.233%, 2.072%,
  1.025% over those four. It vanishes on refinement, so a single-resolution check would have read it
  as a discretisation error and given it a tolerance. Reaching 0.1% would take about **984** cells;
  harmonic is there at twelve.

  And that rate is itself an equality, which asserting `≈2` had obscured. The wrong resistance is one
  face's worth against `(n−1)/2` cells of the two layers, so the relative error is exactly
  `2(1/H − 1/A)/(1/k_a + 1/k_g) / (n−1)`. Multiplying it by `n−1` gives **0.973670** at all four
  resolutions and matches that closed form to `1e-12` — so the coefficient is pinned, not just the
  order. The first draft asserted the consecutive ratios were 2 within 0.1 and passed at 0.091,
  because the exact ratios are `23/11`, `47/23` and `95/47`: a tolerance absorbing a systematic
  offset it could have predicted.

- **Diffraction, and where scalar theory stops being true.** `single_slit_intensity` and `slit_zero`
  in `pantometry-optics`, and `Cavity::obstruct` in `pantometry-em` for a perfect conductor with a hole in
  it. `crates/pantometry/tests/a_slit.rs` holds one against the other.

  `sinc²(π a sinθ/λ)` rests on **Kirchhoff's** boundary condition — the field in the opening is the
  incident field — and Maxwell's equations do not say that. The metal carries currents and the
  opening's field is perturbed over a strip of order `λ` at each edge, so the fraction of the
  aperture that is wrong goes as `λ/a`:

  ```text
    a = 12λ    0.0057     largest absolute difference in normalised intensity
    a =  6λ    0.0125
    a =  3λ    0.0311
    a =  1λ    0.2772     wrong by more than a quarter of the pattern it predicts
  ```

  A factor of **48**, monotone. And the first dark fringe lands at `sinθ = λ/a` to 0.66% at twelve
  wavelengths — the sharpest number in the pattern, and the one with no special constant in it,
  where a circular aperture's `1.22 λ/D` has the first zero of `J₁` in it.

  Below a wavelength a slit has **no** zeros at all, which is what `slit_zero` returning `None` is
  saying: a 0.75λ opening puts 79% of its axial intensity at 30° off axis where a 12λ one puts
  0.05%. That is the difference between an aperture and an antenna.

- **The far field out of a small box.** Not by making the box `a²/λ` deep — hundreds of wavelengths
  for a wide slit — but by recording the aperture plane one wavelength behind the screen and Fourier
  transforming it, which **is** the Fraunhofer limit exactly.

  Magnetic walls make the incident plane wave possible next to a screen, and their mirror images
  would make this a grating — except that the transform integrates over one box width and the
  neighbouring slits lie outside it.

- **Three things about the cost of asking.** A one-cell screen leaks, because the electric edges
  either side of it are shared with the vacuum beyond; two cells hold. A `pattern(width, samples)`
  that marched per call marched the same box eight times where five would do, at 120 s in a debug
  build — the mode four of CI's jobs run; separating the march from the sampling took it to 64 s.
  And the numbers in this entry are measured rather than remembered: the first draft of
  `single_slit_intensity`'s own tests claimed a sidelobe of `4/(2.25π²)` where it is `1/(2.25π²)`,
  and a grazing intensity of 0.8 where it is 0.2545.

### Changed

- **`Solid3D::max_stable_dt` sums the face conductances instead of dividing by a diffusivity**, so
  it is `minᵢ Cᵢ/(dx·Σ_f k_f)`. Two answers moved, and a caller stepping by hand will notice.

  It is now **shape-aware**: `dx²/(2α)`, `dx²/(4α)` or `dx²/(6α)` according to how many axes have
  more than one cell. This domain charged every shape the three-dimensional rate, so a bar-shaped
  block — the shape its own closed-form tests use to check the axes against a 1D answer — paid three
  times the steps it needed. At each of the three limits the sharpest representable mode amplifies
  by exactly `−1.000000000000000`, which is what makes the limit an equality rather than a caution.

  And a filled block is limited by its worst **cell**, which is usually far *looser* than
  `dx²/(6·α_max)`: one aluminium cell inside borosilicate is stable at exactly `k/k_face` = 167/2.21
  = **75.45×** aluminium's own limit, because heat cannot reach a cell faster than its worst face
  delivers it. It can go the other way by at most 2×, and reaching that needs neighbours that conduct
  better *and* store more — impossible for real solids, whose volumetric heat capacity spans one
  order of magnitude where conductivity spans four. Measured at 1.9945× with a doctored specific
  heat, which is the only way to get there.

- **`Solid3D`'s ledger and its placeless heat are weighted per cell.** `Σ CᵢTᵢ` is what a
  two-capacity block conserves, not the mean temperature, and heat off the plain channel spreads to
  a uniform *rise* — in proportion to capacity, not in equal joules, which would have warmed the
  poorer-storing material more and so claimed a location the bus never carried.

- **`ScalarField::rate` on a block is the conductance form**, `(1/Cᵢ)Σ_f G_f(T_f − Tᵢ)`, read off the
  same faces the sweep uses rather than reconstructed as `α·∇²T` from a diffusivity a filled block
  has no single value of. `laplacian` stays `∇²T`, which is what the trait asks for.

- **The scene format's material names live in one table**, `pantometry_world::MATERIALS`, with
  `stainless_304`, `borosilicate` and `ice` added. It was a `match` inside the network builder; the
  moment a second domain wanted a material that would have been two lists that agree until they do not.

- **`Solid3D`'s ledger is in enthalpy**, so a melting front is on the books. A cell holding at its
  melting point while it absorbs 306 mJ per cubic millimetre has taken that heat in and its
  temperature says nothing about it; an audit reading temperature alone would call it a leak of exactly
  that. Every door heat comes in through — the sweep, `deposit`, the plain channel — goes through one
  place, so a phase change cannot be forgotten at one of them.

- **A `Solid3D` that can melt reports a `melted` reading** in mm³, and one that cannot does not. A
  column of zeros in every report tells a reader nothing; the condition is a property of the block's
  materials fixed at construction, so it is not a mode that can surprise anybody mid-run. Without it a
  phase change was the one thing the domain does that a report could not see.

## [0.11.0] — 2026-08-11

**Four new domains, and the list `ARCHITECTURE.md` opened with is closed.** Electromagnetism,
elasticity, fluids and flow through a packed bed are crates on the kernel now, and none of them
needed the kernel or either layer above it to change — rule 4 held ten times. Ten domains, sixteen
crates, 546 tests.

The other half of the release is two tests that make the *relationships* checkable:
`fields_and_rays.rs` puts a Yee grid against Fresnel's algebra, and `loss_and_lumps.rs` puts a
field's decay in a conductor against a `Winding` that answers the same number at every frequency.
Neither closes a gap by making one domain do another's job.

### Added

- **A waveguide, and its whole dispersion relation from one march.** A box with conducting side
  walls and open ends *is* a rectangular guide, so this needed nothing but a launcher that takes a
  transverse profile and a driven source.

  ```text
    f_c = c/2a                    exact, and a property of the cross-section alone
    β = √(k² − k_c²)              0.06% to 2.2% from 1.4 f_c to 3.1 f_c
    α = √(k_c² − k²)              0.4% to 5.5% below cutoff
    v_p v_g = c²                  1.0015, with v_p = 1.14c and v_g = 0.88c
  ```

  A pulse carries a band and the mode's `sin(πx/a)` profile is orthogonal to every other, so one
  march gives `β` at every frequency in it. `λ_g/λ` runs from **1.43** near cutoff to 1.05 at the
  top — a check at one high frequency would pass for a solver with no cutoff at all.

- **Five things about measuring a guided wave, four of them mistakes.**

  **A driven source is not a convenience below cutoff, it is the only way to ask.** An evanescent
  field is a near field: it does not travel, and with nothing driving it it decays to zero, so there
  is no steady state and no spatial profile to fit. The first attempt used an initial-value pulse and
  returned decay constants of 0.285, −0.891 and −0.102 against a closed form of 2.7 — a *negative*
  decay being the measurement saying it was fitting noise. `Cavity::impress` adds to `Ey` on a plane;
  it must be called **after** the step, because the electric update rewrites every interior face and
  erases a source added before it.

  **Mur's boundary is tuned to `c` and a guided mode is not.** Its phase velocity is
  `c/√(1−(f_c/f)²)` — 2c at `1.25 f_c` — so the absorbing faces are worst matched exactly where the
  measurement is most interesting. `β` came out 15% low there against 2.3% at `1.4 f_c`. Truncating
  the march before each reflection returns recovers the band above `1.4 f_c`; below it the near end's
  reflection and the direct pulse overlap and no window separates them. That is a limitation of the
  boundary rather than of the scheme, and which one it is is now written down.

  **A band that does not cover the frequency asked about returns a transform of rounding.** The first
  version read `β` correctly at its carrier, to 0.84%, and 87% wrong two octaves away.

  **A transform window that is not a whole number of cycles leaks between bins**, and with the
  amplitude falling by `e⁻²` across the stations that leakage was the difference between 13% off and
  0.4%.

  **A longer lever on a noisier point is not a better fit.** Widening the stations to hold a fixed
  number of decay lengths made the evanescent measurement worse — 12.4% where the tight span gives
  5.5% — because the far stations sit among the boundary's leakage.

- **What a lumped resistance cannot say, measured.** `Winding` states a resistance; `Conductor`
  solves for one from a shape but is quasi-static, so neither has a frequency in it. A field does.

  `crates/pantometry/tests/loss_and_lumps.rs` measures the amplitude decay inside a conducting
  half-space against the closed form that holds for **any** conductivity, not just a good one:

  ```text
    p = σ/(ωε)                          the loss tangent
    α = ω√(εμ/2)·√(√(1+p²) − 1)        exactly
  ```

  A pulse carries a band, so one march gives the decay at every frequency in it: a discrete
  transform at each depth separates them, and three frequencies come out 0.04%, 0.89% and 1.15%
  from their own `α`. Measuring the *dependence* rather than one number is the point — a single
  driven sinusoid would hide whether the frequency enters correctly.

  Then the handbook `δ = √(2/(ωμσ))` as what it actually is: the `p ≫ 1` limit, off by **55.4%** at
  `p = 1`, 5.1% at 10, 0.50% at 100 and 0.050% at 1000. The 55% is the part worth knowing — a lossy
  dielectric is not a conductor, and using the skin depth on one is not a small error.

  Beside it, copper: 9.22 mm at 50 Hz, exactly a tenfold thinner per hundredfold in frequency, and
  a `Winding` reporting 0.1724 Ω at every one of them.

- **`pantometry-em` carries structures now, and a Yee grid and Fresnel's algebra agree on a
  reflectance.**

  Per-cell media, so a slab or a coating goes into the box, and a *magnetic* boundary so a plane
  wave can exist in one at all: a wave along `z` polarised along `y` has `Ey` tangential to the `x`
  faces, and a conductor there turns it into a waveguide mode with a cutoff and the wrong phase
  velocity.

  `crates/pantometry/tests/fields_and_rays.rs` is the check the two crates existed to make possible.
  `pantometry-optics` gets a reflectance from algebra on two indices; `pantometry-em` marches Maxwell and
  knows nothing about either. They share no code and do not depend on each other:

  ```text
    20 cells per wavelength    −10.543%
    40                          −2.556%     ratio 4.13
    80                          −0.634%     ratio 4.03

    n = 1.5   3.981% against 4.000%
    n = 2     11.041% against 11.111%
    n = 3.5   30.519% against 30.864%
  ```

  The claim is the **rate**, not a tolerance: a single resolution with a loose bound passes for a
  scheme converging at first order, and first order is what an interface placed half a cell wrong
  would give.

  Beside it, two things the algebra alone cannot say. A **quarter-wave coating** of the geometric
  mean index takes a 14.4% reflection to 0.23% — a factor of **61**, and it is interference between
  two surfaces rather than a property of either. And a slab **delays** a pulse by `(n−1)d/c` to
  2.2%, measured on the energy centroid.

- **Three things the cross-domain test found about how to measure.**

  **A monitor plane inside the pulse measures nothing.** The first geometry put a pulse of
  half-width 1.6λ at 4λ and the monitor at 5λ, so the "incident" window opened with the pulse
  already on it: the reflectance came out 34.7% where Fresnel says 14.8%.

  **The peak of an oscillating pulse is a peak of its carrier.** It moves in steps of half a period
  — 0.83 fs here, against delays of a few femtoseconds — so the slab delay read 9% high for that
  reason and no other. The energy centroid is immune to it.

  **A one-dimensional problem does not want a two-dimensional cross-section.** The testbed was
  4 × 4 cells across, every one holding the same field by construction, and the file took **486 s**
  in a debug build — the mode CI's OS matrix runs. Two by one is the smallest grid with an interior,
  gives the same numbers to the last digit, and takes 36 s.

- **An open boundary for `pantometry-em`, and the source that showed why the first measurement of it
  was meaningless.**

  Mur's first-order condition: exact for a wave arriving along the normal, worse away from it as
  `(1−cos θ)/(1+cos θ)`. Not a perfectly matched layer, which would be right to four digits and is
  more machinery than this crate has earned — so the figure is **measured** rather than claimed. A
  line source in an open 24³ box is down to **0.149%** of its energy after two and a half crossings,
  against a conducting box's 101%.

  It unblocks the field formulation of optics: diffraction cannot be computed in a box whose walls
  send everything back.

  **The first source left 34% behind, and the boundary was not at fault.** It was a Gaussian blob
  of `Ey`, which has `∂Ey/∂y ≠ 0` and therefore `∇·D ≠ 0` — a *charge distribution*, whose near
  field does not propagate. No absorbing boundary can remove what was never going anywhere, and a
  measurement using such a source measures the electrostatics instead of the boundary. The source
  is a blob of `Hz` uniform along `z` now, which is divergence-free on the discrete grid exactly.

  Opening a face does not cost the divergence identity, and the test says so: `∂(∇·B)/∂t = −∇·(∇×E)`
  holds for *any* `E`, so changing `E` at a boundary cannot touch it. A measurement that said
  otherwise would mean the boundary was writing to `H`.

  `Cavity::mur_coefficient` is `(cΔt − Δ)/(cΔt + Δ)` — a function of the Courant number alone, which
  a test checks by getting the same number out of vacuum at 2 mm and glass at 6 mm.

- **A tenth domain: `pantometry-fluid`, and the list `ARCHITECTURE.md` opened with is closed.**
  Incompressible Navier–Stokes by projection on a staggered grid — velocities on cell faces,
  pressure at centres, which is the arrangement Yee uses for electromagnetism and for the same
  reason.

  It was the hardest of the three to make checkable, and the crate is built around admitting why:
  fluids has few exact solutions, its schemes trade stability against numerical diffusion, and
  **"it looks like a fluid" is the easiest wrong answer in computational physics to accept.** So
  each test is chosen to be blind to a different mistake, and the docs say which.

  ```text
    Poiseuille      exact to 6e-15 — against the *discrete* parabola
    Couette         exact to 1.2e-14 — and blind to advection entirely
    Taylor–Green    0.272% → 0.065%, ratio 4.17 — the only one that sees the nonlinear term
    uniform flow    unchanged to 0.000e0 after 500 steps
    momentum        drifts 1.1e-15 in a periodic box, because the advection is in flux form
  ```

  Two limits on the step — viscous `dx²/6ν` and Courant `dx/|u|` — and one on the **mesh**: the
  cell Reynolds number `|u|dx/ν ≤ 2`, which no amount of shortening the step fixes. Refused rather
  than run, because past it central differences produce a sawtooth that reads as turbulence.

- **The discrete Poiseuille profile is not the continuum one, and the difference is a closed form.**

  The first test compared to `gh²/12ν` and measured `19/18` of it at six cells. That is
  `1 + 2/36` to the digit — and the scheme was right.

  Any quadratic satisfies the interior equation exactly, so the discretisation contributes nothing
  there. The **wall** does: a no-slip condition imposed by reflecting the first cell makes the
  *linear interpolation* between them vanish at the wall, and a parabola is not its own linear
  interpolation. The discrete answer is `u_j = (g/2ν)[(h²+Δ²)/4 − (y_j−h/2)²]`, whose mean is
  `(gh²/12ν)(1 + 2/n²)`.

  The test now makes two statements instead of one: the discrete form holds to machine precision at
  any mesh, and the gap to the continuum shrinks by exactly `(16/6)² = 7.111×` from six cells to
  sixteen. Comparing to the continuum alone and calling the difference a tolerance would have hidden
  both.

- **A periodic seam that was averaged instead of copied cost 4.7% of a decay rate.**

  With `Walls::None` the `y = 0` and `y = h` velocity faces **are** the same face. `apply_walls`
  set both to their mean — but the update writes the low one and leaves the high one stale, so the
  mean moved it only half as far as the physics did. A half-step lag dressed as a boundary
  condition.

  It showed up as a Taylor–Green decay 4.68% below the discrete Laplacian's own eigenvalue.
  Diagnosed by running the same vortex at an amplitude eight orders smaller: the rate was identical
  to six digits, which cleared the advection and left diffusion holding the bag. After the fix the
  measurement sits **1.02% and 0.26% above** that eigenvalue at the two meshes — which is forward
  Euler's `σΔt/2` to the digit — and 0.27% and 0.07% below the continuum, because the spatial error
  is larger and has the other sign. Two errors, each with a formula, which is what makes the
  comparison mean anything.

- **A ninth domain: `pantometry-em`.** Maxwell's equations in three dimensions, on the grid Yee built
  for them, and the reason that grid is the grid.

  `E` on cell edges and `H` on cell faces, each half a cell from the others, so every curl is a
  difference of quantities that already sit where the curl belongs. The consequence is the point:
  the **discrete** divergence of the **discrete** curl is identically zero — every term appears
  twice with opposite signs — so `∇·B = 0` is a property of the update rather than something the
  scheme approximately respects.

  Demonstrated directly, from both sides. A released mode has `8.6e-17`; two thousand steps later
  it is `1.1e-12`, which is `f64` rounding accumulating linearly and twelve orders below the field.
  And injecting `1.000e-3` A/m of divergence by hand, then running 500 steps of live
  electromagnetism, leaves `1.000e-3` A/m — **1.000000000×**. The update cannot touch it in either
  direction.

  Beside that: cavity resonances against `(c/2)√((m/a)²+(n/b)²+(p/d)²)`, converging at 0.536% →
  0.134% → 0.033% for refinement ratios of **4.002 and 4.001**; the Courant limit `dx/(c√3)`,
  refused above and stable at; and equipartition to 0.998.

- **Four things the electromagnetic tests found, three of them in the tests.**

  **`H` has to be the discrete curl of `E`, not a sampling of the continuum one.** The identity is
  a property of the operators, and a continuum field that happens to be divergence-free does not
  inherit it: setting `H` analytically measured `3.3e-4` where the identity says zero. Applying the
  scheme's own half-step update to a zero field gives both the right initialisation *and* exact
  divergence-freedom, and deletes the hand-derived curl.

  **Joule heating is `∫σE²dt`, not the change in field energy.** The latter is the loss plus
  whatever flowed in from the magnetic field on the same step, and it reported dissipating
  thirty-six times the energy the cavity ever held.

  **Eight half-cycles cannot see second order.** A peak located to the nearest step carries an
  error of `dt`, and eight of them divides it by eight — leaving 0.6% of measurement noise on a
  0.55% signal. The first version reported the *same* frequency at 8, 16 and 32 cells, to seven
  digits, and it took a probe to establish that the scheme was fine and the ruler was not. Two
  hundred half-cycles with the peak interpolated through a parabola resolves it.

  **A 12.8% energy swing is not a defect.** `½εE² + ½μH²` is not what leapfrog conserves — `E` and
  `H` are half a step apart — and the naive sum oscillates by exactly `2 sin(ωΔt/2)`. The first
  version asserted "under 5%", measured 12.779%, and called it a failure where the formula says
  12.8167%. The test checks the closed form now, and that halving the step halves the swing:
  ratio 2.000. An arbitrary bound cannot tell a correct oscillation from an incorrect one.

  A fifth, smaller: a *normalised* divergence divides by `max |H|`, which oscillates through the
  cycle, so comparing it at two instants compares two denominators. The injection test read 1.099
  and then 0.5387 for a quantity that had not moved.

- **An eighth domain: `pantometry-elastic`.** What a shape does under load — `∇·σ = 0` with
  `σ = λ tr(ε)I + 2μ ε` — solved rather than stated, so a stiffness is a property of a geometry the
  way `Conductor` made a resistance one.

  The third elliptic domain here and the first whose unknown is a **vector** at every node. That is
  the only thing that is new, and it is why the operator is assembled from the strain energy rather
  than differenced from the Navier–Cauchy equation: a Hessian is symmetric by construction, which
  conjugate gradients needs, and `∇(∇·u)` differenced on a collocated grid admits a checkerboard
  displacement that costs no energy.

  **Four moduli come out exactly**, at any mesh, because trilinear elements reproduce a linear
  displacement field exactly:

  ```text
    uniaxial stress    E = 68.90 GPa      and nu = 0.330000000 with it
    confined           M = 102.09 GPa     = lambda + 2mu, 1.48x E
    hydrostatic        K = 67.55 GPa
    simple shear       G = 25.90 GPa      where 2G would be the classic factor-of-two error
  ```

  Four combinations of two constants. A solver with `λ` and `μ` transposed reproduces none; one
  whose shear rows carried `2μ` passes the first three and fails the fourth. Checking one modulus
  is checking that a stiffness exists.

  Beside them: the **patch test** on a general linear field with a rotation mixed in — the interior
  lands on it to `1.9e-20 m` — **Clapeyron's `2U = Σf·u`** to `4e-16`, and a null space that is
  exactly the six rigid motions and nothing else, each costing `1e-32 J` against a reference
  strain's `9e-6`.

  **Bending only converges**, from the stiff side, and that is stated rather than tuned away: a
  fully integrated element develops shear where it should flex, so a cantilever runs 0.866, 0.933,
  0.960, 0.981 of `PL³/3EI` as the mesh refines. Reduced integration would cure it and buy hourglass
  modes, which is trading a stiffness error for a *singularity*.

- **Three things the elasticity tests found about their own design.**

  **A shear rig is not simple shear.** Clamp the base, drag the top, read `τ/γ` — and a cube gives
  **0.40 of `G`**, because the sides are free and the block bends. `prescribe_boundary` puts the
  analytic field on the whole surface instead, which is the patch test and removes the rig from the
  answer. Displacement control needed the solver to start its iterate at the prescribed values
  rather than at zero, and the relative residual to be measured against the initial one when there
  are no loads at all.

  **A reaction has to be read by component.** A node on an edge belongs to two faces, so
  `reaction(YLow)` sums the `y`-low share of the `x`-low rollers' 30 N and reports 2.5 N of a force
  it is not carrying. `normal_reaction` is the unambiguous one; the magnitude is not, and the test
  now asserts both — that the normal is zero *and* that the magnitude is not, so the distinction
  cannot quietly stop mattering.

  **Two domains wanted the name `Body`.** `pantometry-mechanics` had it for a rigid one. The prelude is
  where that shows, and only because both are exported there. The elastic one is `Block`.

- **The native viewer opens the espresso run, and two things `runtime/viewer`'s README already
  claimed turned out not to be true of its renderer.**

  `portafilter_flow flow.json` writes the wire format, and `runtime/viewer` reads it — that viewer
  does not link `pantometry` at all, so a run that opens there is the format demonstrating it carries
  enough to draw a run. It did, unchanged.

  What it drew was wrong in two ways, both found by pointing it at a run shaped unlike the
  fixtures: tall, thin, and spanning two orders of magnitude in value.

  **The colour scale was per frame.** `Run::scale_of` computes the range over the whole run, is
  documented as the thing that stops a decay looking like a steady state, is tested — and
  **nothing called it**. `segments` measured the values it had in hand. Water leaving a screen
  clean and arriving at a spout at 83 kg/m³ therefore rendered mid-ramp the whole way down. The
  new test checks the *consumption* rather than the accessor, which is the difference that let
  this survive.

  **The camera was set up for a cube.** `Framing` normalises by the longest side, so a tall thin
  run fills one axis and a fraction of the others: the portafilter came out at 15% of the frame
  height and 0.29% of its pixels. Backing in does not fix it — at this field of view a unit
  subject fills the frame at a distance of 0.35, inside its own bounding box. Distance sets how
  strong the perspective is and that was never the problem. `Camera::fit` sets the **focal
  length**, the projection is linear in it, and one pass is exact. 0.29% to 2.52%, with a cube, a
  plate and a portafilter all landing their furthest corner at 0.850 of the half-frame.

  `--snapshot` takes `--frame N` now, because a run that fills up over its length has nothing in
  it at `t = 0` and a snapshot of that still counts as "the renderer works".

- **The stream out of the spout is drawn as a stream.** A parcel falls the 75 mm from basket to cup
  in 0.12 s against a 0.7 s frame, so sampling parcels there caught roughly none: the first version
  drew the basket beautifully and had nothing between the spout and the cup. What is there is a
  continuous jet, coloured by the concentration leaving the basket — the domain's number, not the
  parcels'.

- **`portafilter_flow`**, which is the picture the domain was built to be able to draw: a shower
  screen, a basket, a body and a spout, with parcels of water leaving the screen, working down
  through the grounds and coming out the bottom darker than they went in.

  Nothing places a streamline. A parcel is advected by the **pore** velocity of the solved Darcy
  field — `u/ε`, not `u`, a factor of 2.2 — and its colour is what it has picked up from the cells
  it crossed. Two checks make that a readout rather than a drawing: the mean transit comes out at
  **0.981×** the closed form `εL/u`, and the load a parcel arrives with agrees with the domain's
  own outlet TDS to **0.9996×**. A Lagrangian tracer and an Eulerian reading, agreeing.

  The hardware and the water share one panel, so they occlude each other properly. What is solved
  and what is only drawn is stated in a table rather than left to look the same.

  With the wall gap the ring's water is through in 4.8 s against the core's 14.4 — **3.01× faster**
  — and arrives carrying 18.0 kg/m³ against the even bed's 82.9. The same water, a fifth of the
  coffee.

- **Two defects the picture found that no test had.**

  **`repack` did not re-solve.** The constructor does, because a quasi-static field read before its
  solve is a field of zeros wearing the shape of an answer. Leaving the solve merely *stale* is
  worse: the field left behind is the previous answer — smooth, bounded, the right order of
  magnitude, and wrong by 42%. Widening the ring to 0.60 left `flow_rate` reporting the even bed's
  flow to the last digit, 1.0000 where the closed form says 1.4188. Every existing test stepped the
  puck before reading it, and a step re-solves. Every mutator ends solved now, and
  `a_mutator_leaves_the_flow_solved` checks all three against their closed forms.

  **An out-of-range probe returned zero, which made the outlet an attractor.** RK2's midpoint lands
  further along than the parcel is, so a parcel one cell from the outlet probes past it. Zero there
  parks it a hair inside the last cell — measured, 19.979 mm of a 20 mm bed — for the rest of the
  run. Nothing about that looks wrong: the stream is steady, the colours are right, the flow rate is
  right, and the water never reaches the cup. Visible only because the example counts arrivals and
  the count was zero.

- **`Substance::stainless_304`, and a basket made of what a basket is made of.**

  `Basket::espresso` used aluminium because the catalogue had it. Steel conducts a tenth as well —
  16.2 W/m/K against 167 — while holding **more** heat per unit volume, 4.0 MJ/m³/K against 2.4:
  a better reservoir and a worse spreader, which is why a group head is brass and a basket is not.

  The wall sets the explicit solver's step, and the limit goes as the diffusivity, so steel's is
  seventeen times larger. `pantometry-porous`'s suite went from **152 s to 5.8 s** in release and 52 s
  to 6.7 s in debug — the mode CI actually runs — and the model got more accurate in the same
  edit. It had made the debug matrix jobs ten times slower than the whole workspace had been.

  Reaching for the metal already in the catalogue cost an order of magnitude in run time *and*
  understated the thermal mass by 65%.

- **A seventh domain: `pantometry-porous`.** Flow through a packed bed, the heat it carries, and the
  dissolution that rides on both. An espresso puck, and also a filter, a catalyst bed, a leaching
  heap and an aquifer.

  Darcy's law is elliptic — `∇·((k/μ)∇p) = 0` — which is the same operator `Conductor` solves for
  electric potential with mobility in place of conductivity. Reproduced exactly: a uniform bed
  gives `Q = kAΔp/(μL)` to 3×10⁻¹⁵, the way a uniform block gives `ρL/A`.

  On top of it, two things that are marched rather than solved: heat advected by the liquid and
  conducted through the bed, and solute dissolving out of the particles and carried away.

  **Grind is two lengths, not one.** Extraction is diffusion out of a particle, so it uses the
  sieve diameter: `K = 4π²D/d²`. Flow is not — a real grind is a coarse mode plus a tail of fines
  that lodge in the gaps and carry the pressure drop, and Kozeny–Carman from the sieve diameter
  says an espresso puck passes twenty litres a second. So the permeability uses a *hydraulic*
  diameter, 1/160 of the sieve one, and that ratio is one of exactly two fitted numbers.

  Those two `d²` point opposite ways, which is why espresso is hard: finer is less permeable
  *and* extracts faster. Neither is coded — they fall out of the two being derived separately.

  Seventeen tests, every one against a closed form, an exact limit or a conservation law.

- **The two things the first version of that crate got wrong**, because each was invisible until
  something was drawn or weighed.

  **Extraction had no equilibrium term**, so `dm/dt = −K·m` depended on nothing outside the
  particle: every cell at the same temperature extracted at the same rate however much liquid
  passed it. A bed with a channel through it then extracted *perfectly evenly*, which is the
  opposite of what a channel does. The driving force is a difference now — `dm/dt = −K(m −
  m₀C/C_sat)` — and a cell the flow rushes past is kept dilute and keeps going while a cell the
  flow avoids fills up its own pore liquid and stalls at `1/(1+β) = 0.577`. That ceiling is the
  entire reason a channelled shot under-extracts.

  **The dose was 52 g in an 18 g basket**, because `solid_density` was coffee's *skeletal* density
  of 1400 kg/m³. A ground particle is itself porous; its apparent density is 600, and with an
  inter-particle porosity of 0.45 that is a 330 kg/m³ puck, which is what 17.6 g in 20 mm of a
  58 mm basket weighs. Three densities, and the middle one is the one this model wants.

- **The statistic that looked like a channel detector and was not.** The obvious measure is the
  spread of the per-cell extraction. An evenly packed bed already sits at 0.105 — water that
  entered clean is loaded by the time it leaves — and a wall gap that halves the yield takes it
  only to 0.128. The signal is a fifth of the noise it rides on.

  `Puck::radial_contrast` divides the axial gradient out: the ring and the core span the same
  depths. It is 1.0000 on an even puck and 1.20 with the gap.

  Nor does the peak extraction rise, which is the story everybody tells. It **falls**, 0.936 to
  0.836, because the channelled basket reached the same weight in 15 s instead of 25. "The channel
  over-extracts" is about the ring relative to the core, and at equal weight the absolute numbers
  go the other way. Reported rather than asserted, so the test does not encode a plausible story.

- **"Finer grind, higher yield" is false at equal time.** Measured: at 25 s a 175 µm bed reaches
  12.9% against a 350 µm bed's 20.2%, its liquid sitting at 10.7% TDS and going nowhere, because
  four times less water crossed it. Pulled to the same *weight* the statement holds and both
  mechanisms show at once — 4.000× the time, exactly the permeability ratio, and 24.1% against
  11.0%.

- **`espresso_shot`**, which is what all of that is for. A basket from the pump to the cup:
  17.6 g in, 38.0 g out, 24.6 s, 19.6% yield, 8.3% TDS. Grind, temperature and pressure each swept
  on their own against the exponent each carries — time to a fixed weight goes as `1/d²` to 3×10⁻⁴
  and as `1/Δp` to 10⁻¹³. Then a gap at the wall, and a portafilter left on the counter.

  `espresso_shot shot.html` draws a vertical cut through three baskets in nine fields as the shot
  runs. A cut and not a volume: nine volumes on fifty frames is 45 MB to animate three baskets
  that are not moving, and the question was about a cross-section.

- **`pantometry_scene::sample_field` is public.** One field per domain is not enough and the limit is
  the trait's rather than the physics'. A bed under flow has a temperature, a pressure, a speed,
  an extraction state and a concentration on the same grid, all true at once, and `as_field` can
  nominate one. A caller that wants the other four builds them here.

  Which one `Puck` nominates is a decision with the same shape. It is **extraction**, not
  temperature: a bed is isothermal unless somebody deliberately cooled the basket, so a
  temperature panel is a flat rectangle on every ordinary run — a picture that renders, looks
  fine, and carries nothing.

- **`Basket` replaces ten positional arguments.** Two of them were `Length` and two were
  temperatures, so transposing a pair compiled and ran; `clippy::too_many_arguments` had been
  silenced to allow it. `Basket::espresso()` is a conventional double basket and `..` covers the
  nine things you are not asking about.

  It also made the geometry statable. The basket radius was "the largest circle that fits", which
  puts the metal in the *corners* of the grid — the right heat capacity in the wrong shape, and a
  cut through the axis crosses none of it. That was invisible until the first cross-section came
  out flat.

- **`scenes/18-an-espresso-shot`**, so all seven domains have scenes, and `DynamicViscosity`,
  `MassFlow` and `Concentration` in `pantometry-units`.

- **The scene format carries a version**, and `--check` validates a file without running it.

  `Scene::format`, where **absence means 1** — what every scene written before the field existed
  is, so old files are readable by construction rather than as a special case. A version this
  build cannot read is refused rather than half-run.

  `deny_unknown_fields` already catches a key that was *added*. It cannot catch a key whose
  **meaning changed** — same name, same type, different semantics — and that is the gap a version
  number closes. The promise attached is narrow and stated: within one version, a file that loads
  today loads tomorrow.

  `pantometry-world --check scene.json` parses and builds without running, and reports a parse failure
  as `file:line:column` with the keys that were expected. CI runs it over every shipped scene, so
  it is not the one entry point nothing exercises.

  This was the editor's blocker: an editor writes files, and the format had changed under its own
  users once already, silently, when `mode` became `release`.

- **`runtime/gpu`**: `Solid3D`'s seven-point stencil as a WGSL compute shader, implementing
  `Domain` so it drops into a `Simulation`. Its own workspace, and the rule is that **the CPU
  domain is the reference and this is a cache of it**.

  Measured, 400 steps: 3.5× at 16³, 24× at 32³, 85× at 48³, **191× at 64³**. The GPU column is
  flat at ~0.055 s — it is bound by dispatch overhead, not the stencil, so what the table really
  shows is the CPU's `n³` growing away from a constant.

  WGSL has no `f64`, so this is a *different computation*, not a faster one. It conserves to
  `5.0e-11` where the CPU holds `9.1e-15`, which is below `Simulation`'s default `1e-9` audit — a
  scene using it must loosen `conservation_tolerance_for(ENERGY, ..)`, and `GpuSolid` declines
  `books_balance` for the same reason.

  Reductions stay on the CPU. A mean summed with atomics depends on which workgroup finished
  first, and addition is not associative, so `ledger` reads the grid back and sums in index order.

### Fixed

- **Single precision was never the problem; spending it on an offset was.** The GPU buffer first
  held absolute kelvin and diverged from the reference by `1.4e-3` after two hundred steps — a
  thousand times what accumulation predicts.

  The update is `centre + F·(sum − 6·centre)`. Near 293 K that `sum` is about 1759, where `f32`'s
  resolution is `1.2e-4`, and the difference being extracted is of order `1e-3` K. Subtracting two
  numbers that agree to five digits keeps **less than one digit** of the answer, every step.

  The buffer holds `T − T₀` now, where the same numbers sit near 1 K and the subtraction keeps
  about four digits. Divergence `1.449e-3` → `8.7e-7`, conservation drift `7.4e-7` → `1.2e-10`.
  The stencil is linear, so subtracting a constant commutes with it exactly and the fix cost
  nothing.

- **`pantometry_view::gltf`**: a frame as glTF 2.0, so Blender, three.js, Omniverse, a USD pipeline
  or macOS Quick Look can open a result. **No new dependency** — glTF is JSON with the binary
  base64'd into a `data:` URI, and this crate already writes JSON by hand, so `pantometry-view` still
  depends on exactly one thing and it is `pantometry-scene`.

  Paths become `LINES`, points and 3D fields become `POINTS`, and a 1D or 2D field becomes
  **nothing, reported**: a line of samples is a graph, not something to put in a 3D viewer, and
  `Exported::skipped` says which panel and why. Eight of the seventeen shipped scenes export
  geometry; the other nine are refused with a reason rather than written as an empty file.

  One frame, not the run. glTF animates node transforms and morph targets — a thing moving or a
  mesh deforming between fixed vertex counts — and a field whose values change is neither.
  Encoding a run as an animation would mean choosing a lie about what is moving.

  Checked against what a loader enforces rather than against whether it looks plausible: buffer
  length against the decoded blob, accessor counts against the geometry, indices read back out of
  the bytes and range-checked, four-byte alignment on every view, `min`/`max` on every `POSITION`,
  and the hand-written base64 against the RFC 4648 vectors.

- **`runtime/viewer`**: a native wgpu window for a run — rotate, zoom, scrub — in **its own cargo
  workspace**, excluded from the library's. Measured: the library resolves 12 external crates, the
  python bindings 15, and this 86.

  It does **not depend on `pantometry`**. It reads the JSON a run wrote and nothing else, so "the wire
  format carries enough to draw a run" is demonstrated rather than claimed.

  `viewer-core` holds everything a renderer gets wrong the same way twice — one colour scale
  across the run, one framing across the run, and a projection that clamps a point behind the eye
  instead of turning it into a streak — with seven tests against real run files.
  `--snapshot out.ppm` renders one frame headlessly, because a window nobody can photograph only
  proves the program did not panic.

- **`PanelData::Paths`** in `pantometry-scene`, and a **`layout`** view in `pantometry-view`: runs of
  connected points in space — a ray through a lens train, a trajectory, a field line — drawn
  rotatably and depth-sorted beside the field and body views.

  The third shape, and it took an optical bench to need it. A field is defined everywhere and a
  body is somewhere; a **path** is a thing that went from one place to another, and drawing a
  traced ray as a scatter of its vertices loses the one property that makes it a ray.

- **`optical_bench`**, an example that draws the *instrument* rather than a graph. A doublet, a
  fold mirror turning the axis through 90°, three field angles, and an image plane — prescribed,
  traced, refocused, and then **bent** until the spot lands inside the Airy disc:

  ```text
    effective focal length   97.8 mm   measured from a traced ray, not assumed
    RMS spot, as prescribed  209.2 um  47x the diffraction limit
    RMS spot, bent           0.57 um   0.13x — diffraction-limited
  ```

  `cargo run --release --example optical_bench bench.html` gives a layout you rotate in a browser.

- **`busbar_rating`**, an example shaped like an engineer's working day rather than a
  demonstration. A bolted busbar joint, geometry to production yield:

  ```text
    contact resistance   3.440 uohm, of which 37% is the joint itself
    thermal path         0.0140 W/K, from the network's solved balance
    continuous rating    445.8 A  (2.23 A/mm2, 105 C limit, 40 C ambient)
    thermal runaway      1018.8 A — 2.29x margin
    yield at nominal     36.6% of 20 000 units
    derate for 99.9%     384.6 A, 86% of nominal
  ```

  Every step has a closed form behind it: `rho L/A` for the bar, Maxwell's `rho/2a` for the
  constriction as a limit the solve is shown converging on, `dT = I²R20/(g − I²R20·alpha)` for the
  electro-thermal fixed point, and that expression's pole for the runaway current.

  The finding it is built to make is the last two lines. **A rating computed from nominal values
  is a coin toss in production** — 36.6% here — and the derating that fixes it is what the Monte
  Carlo is for.

- **Two 3D examples**, `heat_in_three_dimensions` and `room_in_three_dimensions`, run by CI like
  the rest.

  The first is built on the closed form for an instantaneous point source: the peak falls as
  `t^(-d/2)`, and **that exponent is the dimensionality**. A bar gives `-1/2`, a plate `-1`, a
  block `-3/2`. Fitted at `-1.514` inside the window where the source is still a point and the
  block still looks infinite, against `-1.362` before it and `-0.245` after — so the window is
  demonstrated rather than asserted.

  The second is `room_modes` with a ceiling: the floor-to-ceiling mode at 71 Hz that a floor plan
  does not have *at all*, and a mode count checked against Weyl's three-term estimate.

- **`Solid3D`, `Reading`, `Bodies` and `Tolerances` in the prelude.** Each was reachable only
  through `pantometry::core` or its own crate, and each is a type a consumer meets while building a
  frame, setting an audit or writing a domain. Found by examples reaching for them.

- **A `volume` view** in `pantometry-view`: a 3D field is raycast — trilinear sampling, front-to-back
  compositing, rotatable with the same camera the bodies view uses — **beside** the slice montage
  rather than instead of it. A render shows shape and a reader cannot get a number back out of it;
  a montage is quantitative and unreadable as a shape. `ARCHITECTURE.md` gap 6, and the answer
  turned out not to be a depth buffer.

  The opacity transfer function is chosen from the run's own range: transparent in the middle for
  a signed field, or a standing wave renders as a solid block; transparent at the low end for a
  one-sided one, or a block at ambient does the same for the opposite reason.

  When a feature occupies less than 3% of the frame the caption **says so**, with the figure. A
  single hot cell in a block of 729 is a small bright dot and everything else is transparent,
  which is correct and reads exactly like a broken renderer.

## [0.10.0] — 2026-08-10

Three domains gained a third dimension and electricity gained a field, so
`ARCHITECTURE.md`'s six gaps are down to one. **Breaking**, in two places, both in
`pantometry-scene`: `Extent::new` takes an `nz`, and `Panel::grid` returns a triple. Everything
else is additive.

### Added

- **`Solid3D`** in `pantometry-thermal`: conduction through a block in three dimensions. A seven-point
  stencil on cubic cells, insulated faces by mirroring, and the explicit limit `dx²/6α` — a third
  of `Bar1D`'s, because the limit tightens with every axis.

  The first 3D field domain, and the thing a bar cannot do: heat spreading *sideways* out of a
  hot spot is the whole job of a spreader plate, and a one-dimensional model has nowhere for it
  to go but along.

  Checked against the **exact eigenvalue of its own discrete operator**, not against a second
  implementation. A separable cosine mode on a cell-centred grid with mirrored faces decays by
  precisely the same factor every step, so the test is an equality at machine precision. Then the
  continuum: that discrete rate approaches `α·π²/L²` at second order, checked as a *rate*, since
  a first-order scheme also converges.

- **`Domain::books_balance`** in `pantometry-core`: an opt-in claim that a domain's ledger changes by
  exactly what it took from the bus minus what it published. A domain that makes it is checked
  **on its own scale** every step, rather than inside the sum of every ledger.

  The failure that closes, demonstrated in `per_domain_books.rs`: a domain holding a microjoule
  beside one holding a kilojoule loses a fifth of itself, and the total moves by `2e-10`. No
  tolerance catches that — tightening to `1e-12` refuses the run for floating-point noise long
  before it can see a leak of that shape, because the problem is the scale and not the number.

  Opt-in because not every honest ledger is an exact one. `LumpedMass` loses heat to an
  environment that is not on the bus; that is a boundary being modelled, and a check that accused
  it would be the wrong check. Every other domain in the workspace takes the claim and passes.

  `Exchange::traffic` and `Exchange::total_published` are what make it attributable: the scheduler
  visits domains one at a time, so the bus traffic between the snapshot before and the snapshot
  after belongs to exactly one domain.

- **`Tolerances`** in `pantometry-core`, and `Simulation::conservation_tolerance_for`: a relative
  tolerance **per conserved quantity** rather than one for the whole simulation.

  The failure it closes, demonstrated in `per_quantity_tolerances.rs` in both directions: a
  Barnes-Hut tree gives up exact momentum by construction, so at `1e-9` a correct run is refused;
  loosen to `1e-6` and a real energy leak passes. A quantity's achievable accuracy is a property
  of the scheme carrying it, and different quantities in one simulation are carried by different
  schemes.

  `audit_with` is the per-quantity form; `audit` keeps its signature and delegates, so no existing
  caller changes. `Violation` now carries the tolerance that *actually applied* rather than the
  default. A `BTreeMap` inside, because a violation's message must not depend on the order a
  builder was called in.

  A scene can set them too, as `tolerance_for`, and a channel name the kernel does not have is
  **refused** rather than ignored — the same failure `aluminum` for `aluminium` produced in this
  format once already, where one character turned off the check the library exists for.

- **`Conductor`** in `pantometry-electrical`: current as a field. `∇·(σ∇φ) = 0` solved by conjugate
  gradients on a block with two electrodes, `J = −σ∇φ` read off it, and the dissipation as
  `∫σ|∇φ|²dV`. **Nobody states a resistance** — it comes out of the shape.

  For a uniform block it comes out as `ρL/A` to machine precision, which is what makes it
  checkable; for a notched one it comes out as whatever the notch gives, which is the point.
  Series adds resistances and parallel adds conductances, and neither is coded — both fall out of
  the same solve, on materials four orders of magnitude apart so a face conductivity that used an
  arithmetic mean instead of a harmonic one would show.

  `V·I` equals `∫σ|∇φ|²dV` to machine precision — Tellegen's theorem, and the sharpest single
  statement that the discretisation is self-consistent, since the two are different sums over
  different things.

  The first **elliptic** domain here, and the first whose failure mode is a solver rather than a
  stability limit: a solve stopped at its cap returns a field that is smooth, bounded and shaped
  exactly like an answer. `step` refuses one that did not converge, `residual` is a *reading*, and
  `with_solver` exists so the refusal path can be provoked from a test.

- **`Resistivity`, `Conductivity`, `ElectricField`, `CurrentDensity`** in `pantometry-units`, with
  `product!` declarations for `J = σE` and `E = ρJ` — those lines compiling is the check that
  (S/m)·(V/m) is A/m².

- **`DomainSpec::Conductor`** and scene 17 — a copper busbar with a notch. Seventeen scenes.

- **`Hall`** in `pantometry-acoustic`: the wave equation in three dimensions. A staggered grid with
  pressure on nodes and velocity on the faces between them, rigid surfaces, and `dx/(c√3)`.

  It is not a more accurate `Room`. A floor plan **does not have** the vertical and oblique modes
  — not less accurately, at all — and a 2.4 m ceiling puts the first one at 71 Hz. The mode count
  also grows as `f³` rather than `f²`, which is why a real room's resonances merge into a hiss
  where a two-dimensional model keeps them separable much further up.

  Checked against the rigid-wall mode frequencies, which are exact, and against a second-order
  convergence rate measured **across three doublings** — 13 to 97 nodes falls 54.6×, against 64×
  for second order and 8× for first. One doubling was tried first and proved nothing: the
  per-doubling ratios bounce between 2.9 and 5.6, because "worst departure over a run" is a
  maximum and therefore noisy. `Room`'s own convergence test reaches the same conclusion by the
  same route.

  It carries the leapfrog startup fix from its first line rather than inheriting the `O(h)` defect
  `Room` and `Tube` shipped with, and a mutation confirms three separate tests would catch it.

- **`DomainSpec::Hall`** and scene 16 — the same 4.4 × 3.1 m room with a 2.4 m ceiling, released
  in its oblique (1,1,1) mode. Sixteen scenes.

- **`Extent` and `PanelData::Field` gained a third axis** in `pantometry-scene`. Breaking:
  `Extent::new` takes `nz`, `Panel::grid` returns a triple, `PanelData::Field` carries `nz`.
  `Extent::volume`, `Extent::count`, `Extent::dimensions` and `Panel::slice` are new.

- **`DomainSpec::Block`** in `pantometry-world`, and scene 15 — a hot spot in a 9×9×9 aluminium block.
  Fifteen scenes now, all run by CI through the real binary.

- **A `slices` view** in `pantometry-view`: a 3D field is drawn as every z-slice at once rather than
  one plane behind a slider, because a viewer who never touches the slider would see a picture of
  a solid that was really a picture of one plane. The filmstrip has no room for a montage, so it
  draws the middle slice and labels it `z-slice 5/9`.

### Removed

- **`Solid3D::as_bodies`.** It existed as cover for the capture gap below — a way to get a block's
  cells out as a point cloud when a field would have come back as one slice. With `Extent` now
  three-dimensional the cover is unnecessary, and it was never free: a domain that is two shapes
  at once makes the picture depend on whether somebody remembered to set an extent, which is a
  mode nothing announces. It is a field, and only a field.

### Fixed

- **`Solid3D::max_stable_dt` was documented as a limit and read as a recommendation.** At exactly
  `dx²/6α` the sharpest mode the grid can hold has an amplification factor of `-1`: marginally
  stable, so it flips sign every step and never decays. A point source excites it as hard as
  anything can, and the peak comes out **1.96×** the closed form — from a scheme that never
  diverges and whose conservation audit is exact to the last bit. At half the limit it is 1.005×.

  Nothing is wrong with the limit; it is a stability limit and stability is all it claims. The
  documentation now says so, and `heat_in_three_dimensions` runs at half of it with the numbers
  for all three cases in its header.

- **A quasi-static domain reported an answer before it had one.** `Conductor::new` left the
  potential at zeros, so the first captured frame reported a resistance 24× below the floor
  `ρL/A` puts under it — beside a residual of `inf` that nothing was reading. A quasi-static
  domain has no state before its solve, so it solves at construction now.

- **A three-dimensional field was captured as its `z = 0` face**, silently. `Extent::samples` was
  a pair and the sampler built its position as `(u, v, 0)`. For six domains that was exactly
  right; for the seventh it produced a 9×9 plane of a 9×9×9 block — a perfectly plausible picture
  of a block, two thirds of the samples missing, nothing anywhere to say so.

  `FRICTION.md` 23, and the lesson is not about `Extent`. A layer's assumptions are only visible
  from below: `pantometry-scene` names no domain and succeeded at that, and could not have discovered
  that it assumed flatness, because everything it had ever been handed was flat.

- **A test that was green and measured nothing.** Replacing `Solid3D`'s stability constant `1/6`
  with the one-dimensional `1/2` left nine of ten closed-form tests passing — none of them ever
  excited a mode sharp enough to care. `the_limit_is_where_the_sharpest_mode_stops_growing`
  releases the fastest-alternating mode the grid can hold, steps it sixty times at exactly the
  reported limit, and measures. It also asserts the limit is *marginal*, since a scheme that
  damped that mode comfortably at its own limit would be leaving stability unused.

## [0.9.0] — 2026-08-10

### Added

- **`pantometry-view`**, the eleventh crate and the top of the three layers: a filmstrip as SVG, a
  self-contained HTML report, a CSV of every domain's scalars, and the frames as JSON. No
  dependencies — SVG and HTML are text, so a `format!` and a file write is the whole renderer.

  **The view is chosen by the shape of the data**: scalars over time become a chart, a 1D field a
  profile, a 2D field a heatmap, points in space a rotatable scene. Its tests are driven by frames
  written out by hand rather than by a simulation, which is the only way to tell "a heatmap
  because the data is a 2D grid" apart from "a heatmap because that domain was a room".

  Every view holds one scale for the whole run. A picture that renormalises per frame makes a
  decay look like a steady state, and it is what you get if you do not think about it.

- **`pantometry-scene`**, the tenth crate and the middle of the three layers `ARCHITECTURE.md`
  describes. `Placement`, `Extent`, `Frame`, `Panel`, `PanelData`, `capture`, `settle_framing`
  — where a domain sits, and what one instant of a run looks like.

  It **names no domain**, and `knows_no_physics.rs` demonstrates that rather than asserting it:
  the test defines a physics inside the test file — a field, two bodies and a reading, in a crate
  `pantometry-scene` cannot possibly know about — places it and captures all three shapes. If a
  physics invented in a test comes back whole, a real one costs one crate and nothing else moves.

  Both crates were modules in `pantometry-world`, which is `publish = false`. A consumer who could
  state a simulation and run it could reach neither the shape of the answer nor any view of it —
  the largest gap between what was built and what was usable.

- **`Domain::readings`** and **`Reading`**: the named scalars a domain has when it has no picture.
  Eight of the fourteen shipped scenes contain a domain that draws nothing — a heater, a lamp, a
  winding, a thermal network — and for several the scalar *is* the result.

- **`Domain::as_bodies`** and the **`Bodies`** trait: count, position, a value to colour by, and a
  *real* wall or `None`. The counterpart to `as_field`, which covered only the domains that are
  continua. `FRICTION.md` finding 11, recorded and unfixed for months, and paid the moment the
  layers were separated: a scene layer that must name three physics to find out where anything
  *is* needs editing every time a fourth arrives.

  The trait draws a line the old code could not. A periodic cell is a boundary condition and the
  domain reports it; an orbit's box is a property of the picture, and nothing physical sits at
  its edge, so a view measures that one over the whole run instead of being told.

- **`Simulation::domains`**: enumerate what is in a simulation. There was no way at all — a
  caller could ask for a domain *by name*, which is no use to a layer that must visit every one.

- **`ScalarField::unit`**: two characters for a legend, and the fifth place a layer had been
  matching on domain types to get something a domain already knew.

- **`Pose`** in `pantometry-core`: a rigid motion, rotation and translation, no scale and no shear.
  An isometry preserves every distance and angle exactly, which is the only class of placement a
  physics can be moved by without its physics changing. The first test is
  `placing_something_cannot_change_a_distance`.

### Fixed

- **A second domain could empty a channel another had already emptied, and the audit could not
  see it.** `Exchange` counted takes but nothing compared the count across a turn, so two
  consumers of one channel each reported a consistent ledger while the amount was delivered
  twice. `Simulation::sweep` now compares takes per channel across the turn and raises a
  `Violation` naming the channel.

- **`Bar1D`'s field was labelled `"C"` and returns kelvin.** The application converted before
  drawing, so the offset and the label were applied in the same expression and nothing could
  disagree with anything. The field now says `"K"`, which is what the cells hold; the celsius a
  picture wants is a view's conversion, and `FRICTION.md` 22 records that the library gives it
  nowhere to live yet.

## [0.8.0] — 2026-08-09

### Added

- **`Ensemble`** in `pantometry-core`: many independent samples, run in parallel, with an answer that
  does not depend on how many threads produced it. The other axis of parallelism — `TreeNBody`
  splits one evaluation across cores, this splits many evaluations, which is the shape a Monte
  Carlo study and a parameter sweep both have. Measured 8.74× on sixteen threads.

  It is bit-for-bit across thread counts because two existing decisions meet: `Rng::for_index` is
  stateless and index-addressed, so sample `i` draws the same numbers wherever it runs, and
  results land in a slot chosen by index rather than being appended. The failure this avoids is
  the usual one and is nearly undetectable — a Monte Carlo drawing from a shared generator gives
  a different answer per core, and the difference looks exactly like statistical noise.

  `estimate` folds in **fixed-size blocks**, so a study is bounded by its block count and not its
  sample count: ten million samples in kilobytes, tested. The block size is fixed rather than
  derived from the thread count on purpose — a per-thread split combines a different number of
  partial sums on four cores than on sixteen, and floating-point addition is not associative.
  `Ensemble::blocks` is public so a caller building a histogram or a quantile has the same
  discipline available.

  Welford within a block, Chan's merge between them, rather than `sum(x²) − n·mean²` — which
  subtracts two large nearly-equal numbers and loses every digit exactly when a Monte Carlo has
  converged. Pinned against a case with an exact answer.

- **`Estimate`**, with `mean`, `standard_error`, `samples`, and `standard_deviation()` kept
  distinct from the error on the mean, because confusing those is the usual way to misreport a
  Monte Carlo result.

### Changed

- **`Fluid` is about 1.95× faster**, in changes that move no result: the Lennard-Jones potential's
  loop invariants hoisted out of the pair loop, the force quotient taken once instead of twice,
  and the periodic wrap skipped where it provably does not apply. Verified on four platforms
  against a pinned digest. Two techniques that should have helped did not — a cell-ordered copy
  of the positions is *slower*, because the counting sort is stable in index order and the reads
  were never a random gather.

- **`detector_snr`** runs on `Ensemble`, which also fixes a variance it was computing as
  `sum(k²)/N − mean²` — at a mean of 900 that subtracts 1.6e11 from itself to reach 900.

- **`where_the_time_goes`**, a new dependency-free example, because this workspace had never
  measured itself and every claim about which loop mattered was a guess. It takes the best of
  five trials: consecutive runs vary by 8%, which is wider than several differences that were
  nearly reported as wins.

## [0.7.0] — 2026-08-09

### Added

- **`ThermalNetwork::path_conductance(node, at)`** — the conductance of the whole heat path from
  a node to ambient, as the slope of its own solved balance. Exact and operating-point
  independent when nothing radiates; the local slope when something does, which is the right
  answer because everything asking for this is asking a derivative question.

- **`Volume::cm3`/`mm3`/`m3`/`litres` and `Area::cm2`/`mm2`/`m2`.** Building a three-node network
  was six lines of `Volume::from_si(x * 1e-6)`, because the constructors take dimensioned types
  and the numbers a person has are cubic centimetres.

### Fixed

- **A number this workspace was quoting was wrong, and the fix is an API rather than an edit.**
  `runaway_current`'s documentation said a motor's threshold falls from 4.95 A to 4.11 A once the
  joints are counted. The 4.11 is a *convection-only* path: the real one is 0.220 W/K rather than
  0.203, because the housing also radiates at its operating temperature, and the true threshold
  is **4.28 A**. The hand-assembled formula understated the margin by 4%.

  Found by a sizing tool written against the published 0.6.0 — a consumer deliberately unlike
  `pantometry-world`: no scenes, no rendering, no fields, asking for settled answers rather than
  stepping. It had to assemble that conductance out of numbers the network already held, which
  is `FRICTION.md` 20 and is what `path_conductance` now answers.

  It also cross-checked something worth keeping: the tool's fixed-point iteration on
  `steady_state` lands the winding at 99.0 °C, and scene 13's marching with a between-frames
  feedback lands at 99.02 °C. Two unrelated routes to the same coupled answer.

## [0.6.0] — 2026-08-09

### Added

- **`Domain::as_any_mut` and `Simulation::domain_as_mut`.** A caller could read a domain and not
  write one, which made a whole class of coupling closable from *nowhere*: not inside the step
  loop by design, and not outside it by omission.

  The case is a copper winding whose resistance rises 0.393%/K. Its temperature lives in
  `pantometry-thermal`, its resistance in `pantometry-electrical`, and neither can see the other's state
  — correctly, since domains meeting only on `Exchange` is the property the crate split defends.
  The caller between frames can see both.

  This does **not** weaken that rule: it is about what happens inside `step`, and this runs in
  code holding `&mut Simulation` that could drop the domain and rebuild it. It is also
  deliberately *not* a state channel on the bus, which a true in-loop coupling would need and
  which stays undecided.

  `as_any_mut` defaults to `None`, so a domain that forgets it is silently unwritable — the
  opt-in hazard of `FRICTION.md` findings 7 and 12, handled in the same change rather than
  rediscovered: all twelve implementors got the counterpart beside `as_any`, including the
  `Box<dyn Domain>` forwarding impl. `FRICTION.md` 18.

- **`Winding::dissipation_at(T)` and `resistance_at(T)`** — pure functions rather than `Domain`
  methods, so `P(T)` composes for whoever holds both sides. `dissipation()` is now
  `dissipation_at(its own temperature)` rather than a second copy of the arithmetic, checked on
  `to_bits()` at four temperatures.

- **`Winding::runaway_current(g)`**: `√(g/(R₂₀α))`, where `dP/dT` overtakes `dQ_out/dT`. The test
  measures the slope from two dissipations a kelvin apart and asserts the inequality *flips*
  across it, rather than reproducing the formula from itself. `None` for a voltage drive, which
  cannot run away because `V²/R` falls as it warms.

  `g` is the **whole** path to ambient. A winding reaching air through 0.9 and 2.4 W/K of joints
  and then 0.294 W/K of convection has a series conductance of 0.203, and the threshold falls
  from 4.95 A to 4.11 A — 17% of margin a lumped model reports as present when it is not.

- **Scene 13**, the feedback closed between frames and measured. The amplification is `1/(1−g)`,
  checked as a ratio against the same scene without it: 1.281 measured. Convection alone predicts
  1.310; including the housing's linearised radiative conductance at its operating point gives
  1.280, so the 2.2% is radiation stiffening the heat path rather than error.

## [0.5.0] — 2026-08-09

### Added

- **`pantometry-electrical`**, the sixth domain and the tenth crate. `Winding` computes
  `R = ρ(T)·L/A` and publishes `I²R` onto the channel `pantometry-thermal` already takes from, with
  neither crate naming the other.

  It closes a real gap rather than adding a sixth for its own sake. Every other producer of heat
  here answers a question about *something else* that happens to warm a thing — light landing on
  a mirror, a dashpot damping a bounce. A winding is the case where getting hot is the entire
  subject, and until now the workspace's own examples stood a stated number of watts in its
  place. A stated number cannot be wrong, which is another way of saying it is not a model.

  Two mistakes it made, both caught by something other than its author. It first declared its own
  `HEAT = "heat"` channel, which reads correctly and is a *different* channel from
  `quantity::ENERGY` — so it published joules nothing consumed, and the audit named it on the
  first step with the amount. And an infinite reserve turned out not to *fail* the audit but to
  **disable** it: `inf` before, `inf` after, `inf` equals itself, and a winding pouring joules
  into a plate runs green at any tolerance. `Winding::step` refuses that itself, because the
  audit structurally cannot.

  The electro-thermal feedback that causes thermal runaway is deliberately **not** expressible:
  a domain would have to read another's temperature inside the step loop, and `Exchange` carries
  amounts, not state. Resistance is evaluated at a temperature the caller states.

- **`Resistance`** in `pantometry-units`, with `product!(Resistance, Current => Voltage)` — that line
  compiling is the check that ohms times amperes are volts. Plus `Current::a`, `Voltage::v` and
  `Resistance::ohm`/`milliohm` constructors.

- **Scene 12**, `12-winding-heats-a-motor`: the same motor as scene 11 with the watts computed
  from a length of wire rather than stated. The two settle within a fifth of a kelvin of each
  other, so the guess was good — and this is the scene that would have caught it if it had not
  been. Twelve scenes now cover all six domains.

- **Two tests that keep prose from aging.** `documented_version.rs` reads every `pantometry = "x.y"`
  in the documentation and every "the tree is X.Y.Z" in `.claude/agents/`, and compares them
  against `CARGO_PKG_VERSION`. Both fail if they find *nothing*, because a check that stopped
  matching would pass forever. The `invariant-guard` line they now cover had been a release
  behind twice running, in the one file whose subject is checking things.

## [0.4.0] — 2026-08-09

### Added

- **`ThermalNetwork::steady_state(power)`** — where a network settles, solved rather than marched
  to. The question a designer actually asks is *will the winding survive*, and stepping to it is
  slow and approximate: reaching a part in a thousand takes about seven time constants of
  explicit Euler, each accumulating its own error. This solves the same balance the step loop
  converges to, exactly, and returns a `SteadyState` read with the same `Node` handles.

  Not the implicit stepping this workspace declines to have: `Domain::step` is unchanged, the
  kernel is untouched, no schedule learns anything, and the network is not modified by being
  asked. Newton, because radiation makes the balance `T⁴` and a single solve would answer the
  linearised problem — the mistake `LumpedMass::equilibrium_rise` exists to correct on one body.
  The Jacobian's radiative part is the `linearised_loss_conductance` the step limit already uses.

  Refuses a network where no node loses heat to an environment: it warms without limit, so there
  is no steady state, and a plausible finite number would be the worst possible answer.

  Exposed in Python as `steady_state(name, watts)`, where it is worth more than in Rust — a
  marching loop crosses the binding once per step.

### Fixed

- **The Newton bound was set from the wrong measurement and refused a kilowatt.** It was eight,
  twice the worst case the first tests exercised — but none of them loads a node hard enough for
  the `T⁴` term to dominate. At ambient the radiative slope `4εσAT³` is tiny against what the
  balance needs, so the first solve overshoots enormously and Newton walks down at the `3/4`
  ratio a quartic gives: twelve iterations at a kilowatt, sixty-six at a terawatt. Now 100, from
  counting them, with the table in the method's documentation.

  Found by instrumenting the iteration count, which also showed the loop had no way to *report*
  exhausting itself — it returned the last iterate, a plausible temperature for a balance that
  was never struck. It returns a `Violation` now, and the test that would have caught the
  original bound checks a radiation-dominated solve against a root found by **bisection**, which
  shares no arithmetic with Newton.

## [0.3.0] — 2026-08-09

### Added

- **`ThermalNetwork`**, the third domain in `pantometry-thermal`: n lumped bodies joined by
  conductances, as *one* domain. Winding, stator, housing — and the **drop across each joint**,
  which is the number that decides whether a motor survives and the one a `LumpedMass` cannot
  give, because it reports the whole assembly as a single temperature. It also expresses a
  *contact* resistance between different materials, which `Bar1D`'s uniform grid cannot.

  One domain rather than a `conducting_to(peer)` on `LumpedMass`, because a conductance carries
  `UA(T₁ − T₂)` and needs **both** temperatures — and domains here meet on an `Exchange` that
  carries amounts rather than state, so neither side can compute the flux alone. Adding that
  method would have broken the property the crate split exists to hold. A network is a single
  coupled system of ODEs, which is what a thermal network physically is.

  Nodes are `Node` handles rather than names, and that is the load-bearing decision. A link
  contributes `+q` to one node and `−q` to another **in the same sum**, so they cancel
  identically and the conservation audit is blind to links *by construction*: a sign error, a
  transposed index or a link dropped altogether passes at machine precision, and the winding
  simply runs at a plausible wrong temperature forever. A handle can only come from a
  constructor, so a dangling link is not representable. `node_named` and `handles()` are the
  bridge for callers building from a file — the JSON scene format and the Python binding both
  resolve names once, at construction, and raise before any stepping happens.

  Seven tests, six of them closed forms, every one per-node or against a formula computed in the
  test file rather than on a total. `n = 1` reduces to `LumpedMass` bit for bit over a 4000-step
  trajectory including `max_stable_dt`, so the new domain inherits every check the old one
  already passes; `linearised_loss_conductance` is shared between them rather than written twice,
  because two copies is the obvious way for that to stop being true.

  Closes #2.

- **`Conductance`** and **`HeatCapacity::j_per_k`** in `pantometry-units`. `Conductance × Temperature
  = Power` and `Conductance × Time = HeatCapacity` are declared, and the declarations compiling
  is itself the check that `UA·ΔT` is watts and `C/UA` is a time.

- **`ThermalNetwork` in the Python bindings**: `add_network(name, nodes=[…], links=[…],
  absorbing=…)`, with `node_temperatures`, `node_temperature` and `heat_flow_w` to read it back.
  A node given `ambient_k` without `area_m2` — or the reverse — is refused rather than quietly
  becoming an interior node that looks like it is cooling and is not. `temperature()` refuses a
  network rather than averaging it, and names the calls that answer.

- **The Python bindings are on PyPI**, as `pantometry`. `pip install pantometry` gets an abi3 wheel for
  Linux x86_64/aarch64, macOS x86_64/aarch64 or Windows x64, plus an sdist to fall back on.
  Built by `.github/workflows/release-python.yml` on a tag, because a wheel built on one
  machine is a wheel for one platform — uploading the Windows one alone would have made
  `pip install` fail on Linux and macOS in a shape that reads as an unsupported platform rather
  than a botched release. Trusted publishing, so no token lives in the repository.

- **A `network` domain in the scene format**, and scene 11: 12 W into a copper winding, out
  through electrical steel and an aluminium housing. The first scene with **nothing to draw** —
  `as_field` declines, because nodes have capacities rather than positions and a conductance is
  not a distance — so the scene test's "produced a panel" guard now takes an explicit list, and
  being on it costs a named check rather than buying a pass.

- **Python bindings**, in `bindings/python`, as their own cargo workspace. `pip install` the
  wheel and `import pantometry`: a `Simulation`, the library's heater, bar and lumped-mass domains,
  and the conservation audit as a `pantometry.Violation` carrying `quantity`, `site`, `before`,
  `after`, `scale` and `tolerance` — addressable rather than a sentence to parse. A refused step
  does not move the clock.

  SI floats at the boundary with the unit in the parameter name, because the dimensional types
  are a compile-time thing Python cannot have, and a runtime wrapper would cost per operation to
  catch an error a Python caller does not make. What crosses instead is the audit.

  A domain cannot be written in Python yet, and the reasons are in its README. Enough to *run and
  audit* coupled physics, not enough to *extend* it.

  Separate workspace because pyo3 brings about fifteen crates and links libpython, and the
  library's twelve external dependencies, its `deny.toml` allow-list and its WebAssembly jobs are
  promises that should not have to accommodate a Python extension. Verified rather than assumed:
  the library workspace still resolves to exactly twelve external crates. An abi3 wheel, so one
  build serves 3.10 upward; CI builds it, installs it and runs its ten tests, each against a
  number computed in the test file rather than read off the simulation.

### Fixed

- **An `O(h)` bias in `ThermalNetwork`'s steady state that the conservation audit could not
  see.** Heat arriving on the bus was added to the absorbing node's temperature *before* the flux
  snapshot was taken, so that node drove its link from an already-raised value. Explicit Euler
  otherwise reaches a steady state exactly — the fixed point is where the right-hand side
  vanishes, with no step-size dependence — so the joint next to the source sat `K·h/C` low:
  predicted 0.0031006, measured 0.0031005, while the far joint and the environment drop were
  exact to six figures. Every total stayed right throughout, because the excess simply landed in
  the neighbour. Found by the series-resistance closed form, not by the audit. The arriving heat
  is a term of the same right-hand side as the fluxes and now joins the same sum.

- **A NaN check in the Python bindings written as `!(x > 0.0)`**, which rejects NaN by the
  negation being true rather than by saying so. The nested binding workspace is excluded from the
  root one, so the `lint` CI job had never reached it and it had gone unlinted since it was
  written. `cargo fmt --check` and `clippy -D warnings` now run in the bindings CI job.

### Changed

- `FRICTION.md`'s header and footer disagreed on how many findings were fixed, and both
  disagreed with the file. Counted: twelve of seventeen. `AGENTS.md` gained a heat-model
  selection table and quotes a CI-run example function rather than a hand-written snippet, with
  the number it prints pinned by an assertion — prose stating a figure that nothing checks is
  how a document goes stale.

## [0.2.0] — 2026-08-08

Breaking, and almost nothing broke: `&str: Into<String>` meant not one existing call site
changed. Everything in it came from the workspace acquiring its first consumer, and then from
two subagents built out of what that consumer taught.

The headline is not the ergonomics. It is that a first-order accuracy defect in the kernel's own
scheduler — in the schedule chosen *for* accuracy — was found by an application comparing a
coupled run against the closed form of its own recursion, having survived every test the library
had while the conservation audit reported clean to 1e-12.

### Fixed — the kernel scheduler

- **`Schedule::Multirate` does not refine a coupled quantity.** `sweep` steps one domain to
  completion before the next, so a quasi-static publisher puts a whole outer step's joules on the
  bus at once and a subcycling consumer takes all of them on its *first* substep. Refining the
  substep therefore does not move the answer: the error is first order in the **outer** step and
  independent of the substep. Measured on a lumped plate under a lamp — 26.2% low at a 300 s
  outer step, 13.8% at 150 s, 7.1% at 75 s, whatever the substep count — and at 300 s the
  schedule chosen *for* accuracy is worse than `Staggered`, with the errors on opposite sides.

  Every one of those runs passes the conservation audit at ~1e-12. The total that crossed is
  exactly right; only its distribution in time is wrong, and a `Ledger` has no representation for
  *when*. The time-domain twin of the reason `audit_transfers` became a per-face check.

  **Fixed** by `Exchange::take_share(channel, dt)`, new in the kernel: `advance` tells the bus
  what interval the sweep covers and a subcycling consumer asks for its substep's share instead
  of the lot. `Bar1D` and `LumpedMass` use it. Error at a 300 s outer step went from 1.89 K to
  0.304 K — from the worse of the two schedules to fourteen times better than the alternative.

  The share is apportioned against the time *remaining* rather than the whole interval, which is
  what leaves the channel exactly empty: `A·dt/T` with both reduced keeps `A/T`, so the last
  substep takes the remainder. Against the whole interval, `n` shares strand `O(n·ε·A)` where
  `audit_transfers` uses an absolute tolerance.

  Built rather than deferred because the recommendation to wait for a second consumer did not
  survive checking: `04-heater-and-bar` already pairs a quasi-static heater with a bar that
  subcycles hard, so a shipped scene had the defect.

### Fixed — what the two new subagents found

- **`Violation::at`'s cases printed an ungrammatical sentence.** They carry a *message* in
  `quantity` and `Display` had no branch for them, so the first error a consumer ever saw read
  "substance has no heat capacity is not conserved at plate: inf". A third branch.
- **`Report` could not be named without a module path**, though it is what `advance` returns.
  Added to `pantometry-core`'s root and the prelude.
- **`Substance` was in the prelude and unbuildable from it**: `bulk` leaves `thermal: None`,
  which `LumpedMass` refuses to step, and the three types needed to supply one were not exported.
  `ThermalProps`, `MechanicalProps`, `AcousticProps`, `ThermalConductivity` and
  `ThermalExpansion` now are.
- **The scene format discarded unknown keys.** `serde` does that by default, which is right for
  a wire protocol and wrong for a saved document: `main.rs`'s own built-in scene kept the
  pre-`release` spelling for two commits and nothing failed, because the keys were dropped and
  the field fell back to its `Default`. Editing them was a no-op that reported success.
  `deny_unknown_fields` on `Scene`, `DomainSpec`, `Release`, `ScheduleSpec` and `Boundary`.
- **The round-trip test could not have caught that.** Both sides of its byte assertion were
  serialiser output, so the hand-written spelling never entered any comparison. It now parses the
  text a person would type and requires each stated value to survive.
- **An unrecognised `finish` produced a silent zero-watt lamp.** The early return also skipped
  `with_reserve`, so the reserve stayed infinite, so `Light::ledger` reported nothing, so the
  audit had nothing to compare — and the scene ran green at `conservation_tolerance(0.0)`, the
  strictest setting expressible, with the lamp doing nothing. One character, `aluminium` against
  `aluminum`. `DomainSpec::build` is fallible now and names the finishes it knows.
- **Two domains could share a name.** `Simulation::domain` takes the first match, so the second
  was never sampled and the first was drawn twice under the second's label and geometry — a
  500 °C bar reported as 20 °C, twice. `World::build` refuses it.
- **A scene whose every domain lacked a field wrote a zero-byte SVG and exited 0**, and "0 KiB"
  could not distinguish that from a legitimate 937-byte strip. The report is a row per *domain*
  now, naming the ones with nothing to draw, with the run-wide extremum beside the final value —
  because a ball that bounced half a metre and one that never moved both end at zero. An empty
  picture is refused rather than written.
- **One colour scale spanned panels of different units**, so a 1 Pa room beside a 7546 m/s orbit
  rendered as an empty bordered square while the numbers beside it looked fine. One extent per
  panel now, still shared across frames.

### Changed — breaking

The first consumer went looking for the API's shape and found the same decision three times:
`&'static str` and `impl Domain` are free when every name is a literal in a test and every
domain type is known when you compile. They are not free for an application. All three are
reversed, and the cost was one order of magnitude smaller than the argument for keeping them —
**no existing call site changed**, because `&str: Into<String>`.

- **`Domain::name` returns `&str`.** Domains store a `String`; constructors take
  `impl Into<String>`, so `Bar1D::new("bar", ..)` and a name read out of a scene file both
  work. `Interface` followed, and `Exchange`'s spatial channel map is keyed by an owned
  interface name.
- **`Report::substeps` is `Vec<(String, u32)>`.** The only breakage across 349 tests: five
  comparisons against string literals. A report outlives the borrow it would otherwise hold.
- **`Simulation::with_boxed(Box<dyn Domain>)`**, plus `impl Domain for Box<dyn Domain>`. The
  simulation always stored boxes; now a caller who built one at run time can hand it over.
- **`Domain::as_field() -> Option<&dyn ScalarField>`** and `Simulation::field(name)`, both
  opt-in and `None` by default, in the style of `as_any`. `ScalarField` was written as the
  interface a visualiser reads a simulation through and was unreachable from `&dyn Domain`, so
  the workspace's own visualiser downcast to concrete types instead — precisely what the
  interface existed to avoid. It no longer names a single domain type.
- **`Room` is in the prelude**, which it should always have been; `Tube` already was.
- **`NBody`, `TreeNBody`, `RigidBody` and `Rolling` answer `as_any`.** All four returned the
  default `None`, so `Simulation::domain_as` could not reach any of them and a renderer got
  nothing back. The failure was a picture with no bodies in it — not an error, not a violation,
  just an empty frame, which is the least debuggable outcome there is. Optics, thermal,
  acoustic and molecular had all opted in because tests inside the workspace had reached for
  them; mechanics had not, because none had.
- **`Bar1D::exposing` takes `impl Into<String>`.** The sweep behind `Domain::name` had matched
  on the parameter being called `name` and missed this one, so a boundary name — the thing two
  domains have to agree on, and therefore exactly the kind of thing that comes from a file —
  was still compile-time only. The kernel's own worked examples were also still teaching
  `fn name(&self) -> &'static str`.

### Added

- **`pantometry-world`** — the first consumer, and not published. Scenes described as JSON, built
  into a `Simulation`, run, and drawn as an SVG filmstrip with no dependency. It exists to use
  the SDK from outside rather than to be a good application, and it reports what that was like
  in `crates/pantometry-world/FRICTION.md`: twelve findings, seven fixed. Five are recorded rather
  than actioned, and the reasons differ: sharing the examples' SVG plotting would mean
  committing to a public drawing API in a workspace whose scope excludes rendering; validating
  a scene's schedule against its domains at build time is something an application can already
  do with `Domain::max_stable_dt`; and a duplicated face count needs no format change because
  `Exchange::publish_on` already refuses the mismatch by name. The report also says what it
  does *not* cover.

  **Ten scenes ship**, in `crates/pantometry-world/scenes/`, covering all five domains,
  and CI runs every one through the real binary as well as the test harness. Three acoustic —
  two room modes and a clap that reflects off all four walls; two thermal, which are the same
  heat told over a plain channel and over a shared boundary and are the argument for
  `Interface` in one picture; two mechanical — four satellites and a bouncing ball whose
  dashpot heat a thermal lump takes; and two molecular, the same 108 atoms at `T* = 0.15` and
  `T* = 1.4`, which side by side are melting. Each has one number asserted, chosen to be a
  property of the physics rather than of the file, and a scene that ships without a claim
  fails the test rather than passing quietly.

  A room can be released as a Gaussian pulse now and not only as a mode, and a panel can hold
  bodies as well as a sampled field — an orbit is a countable number of things at places, and
  rasterising one would invent a continuum it does not have.

  **Bodies are drawn in three dimensions.** The physics always had them: `NBody`,
  `ContactSystem` and `Fluid` all carry `DVec3`, and flattening to a plane was the renderer's
  simplification. They are projected axonometrically now, sorted back to front, with radius
  growing toward the viewer and colour mixed toward the plate for distance — all three are
  needed or the picture is flat however true the coordinates are. Periodic cells get a
  wireframe, because they are a real boundary. The orbit scene tilts its satellites out of one
  plane, which is what makes the third axis carry anything.

  **Optics has a scene**, and it is the fifth `Domain` written outside the library: a
  blackbody lamp on an aluminium mirror whose reflectance falls off in the blue, so the
  colour temperature decides how much of a hundred watts becomes heat. Asserted as the
  difference between 2800 K and 6500 K rather than as one number, because a flat reflectance
  would make `Spectrum` and `SurfaceOptics::absorptance` an expensive way to multiply by a
  constant.

  The scene format couples both ways. On a plain channel, a heater defined **in the
  application** publishes joules and a bar takes them. Over a shared boundary, a beam
  publishes a Gaussian `Flux` onto an `Interface` the bar exposes, and the bar ends up hotter
  in the middle than at the ends by a ratio the scene never states. Both are audited at 1e-9,
  and the totals are checked against `Q / ρVc_p` computed outside the library.

  That closes the gap the report itself had flagged — until then no `publish`, `take`,
  `Exchange`, `Interface` or `Flux` appeared anywhere in the consumer, so the part of the API
  this workspace exists for had only ever been driven from inside. Writing both domains from
  outside needed nothing beyond `pantometry::prelude`.
  Excluded from the wasm, determinism and 1.78 jobs, which are promises the *library* makes to
  the people who depend on it.

### Found and fixed

- **A first-order startup error in `Room`, and in `Tube`.** A mode released from rest follows
  `|cos(2πft)|`, but the gap converged at first order against grid resolution where the
  scheme's interior is second. `released_from` left the velocity at `t = 0`; a staggered
  leapfrog carries it at `t = −h/2`, so the first velocity update travelled a whole step where
  it was owed half. `O(h)`, permanent, and `h` follows `dx` through the CFL condition.

  Fixed: the first velocity update takes half a step. Second order now, and the worst
  departure over 20 ms fell from 0.0528 to 0.00238 at 31 cells. Found by the workspace's own
  application checking itself against a closed form, because nothing inside the library was
  checking a rate.

  Two things came out of the fix. **A test had turned the bug into the specification** — one
  step from rest was asserted to move the pressure by `h²c²∇²p`, where Taylor gives `½h²c²∇²p`
  since `ṗ(0) = 0`; the test was missing the half because the scheme was, and `Tube` had the
  same pair. And **the old startup conserved energy exactly while the correct one does not**:
  with `v = 0` read as the half-step value, `Σ∇·(p∇p) = 0` at a rigid wall makes the first
  step's energy change cancel to the last bit. Starting correctly breaks that by `O(h²)` at
  the first step only — 0.42% at 31 cells, quartering on refinement, and 1e-15 thereafter. The
  old code had bought exact bookkeeping by making the scheme first order, which is this
  workspace's own documented trap appearing a second time in the same crate. Energy is now
  reported against the released state as its datum, with the difference available from
  `Room::startup_adjustment` and bounded at 25% so a real first-step bug cannot hide there.

- **The API is comfortable only when the set of domains is known at compile time.**
  `Simulation::with` takes `impl Domain` and there is no `impl Domain for Box<dyn Domain>`;
  domain names are `&'static str`, so a name from a file has to be leaked; and a renderer
  cannot get a `&dyn ScalarField` from a `&dyn Domain`, so it downcasts and knows every domain
  by name — which is what `ScalarField` existed to avoid. Three symptoms of one position. It
  may be the right position, but nobody chose it.

## [0.1.0] — 2026-08-07

First release. All eight crates published to crates.io together.

Everything below was in this release. Two notes on the publish itself, since they cost time
and are not obvious from the outside:

- crates.io requires a **verified** email address, not merely a registered one, and reports
  its absence as a `400` at the first upload rather than at login.
- New crates are rate limited to a burst of five, then roughly one every ten minutes. A
  workspace of eight publishes five, stops, and has to be resumed — so `cargo publish
  --workspace` is not atomic and a partial publish is the normal outcome, not a fault.

### Added

- **`pantometry-units`** — dimensional analysis with the SI exponents in the type, so `Length +
  Time` does not compile. Const-generic `Qty` and `QVec3`, macro-generated products, and
  unit-bearing constructors as the only place a factor of a thousand may appear.
- **`pantometry-core`** — the kernel. Conservation as an audit (`Ledger`, `audit`, `Violation`),
  fixed-step integrators including `velocity_verlet`, multi-domain scheduling with
  quasi-static, multirate and iterative coupling, deterministic sampling through
  `Rng::for_index`, an accurate discrete Fourier transform, closed-form rigid motion, scalar
  and vector fields, and shared boundaries (`Interface`, `Flux`).
- **`pantometry-optics`** — spectral radiometry, surface optics with Fresnel and coatings,
  Sellmeier and Cauchy dispersion, ray geometry, Airy diffraction and the ideal MTF, Zernike
  wavefronts and aberrated PSFs, angular-spectrum propagation, partial coherence, and a
  detector with the four noises that come with counting photons.
- **`pantometry-thermal`** — lumped masses and explicit one-dimensional conduction, with radiative
  and convective loss.
- **`pantometry-mechanics`** — exact N-body, Barnes-Hut with a quadrupole, penalty contact with
  Coulomb friction, rolling, and rigid-body rotation with Euler's equations.
- **`pantometry-acoustic`** — the linear wave equation on a staggered grid, in a tube and in a
  room, with impedance boundaries.
- **`pantometry-molecular`** — Lennard-Jones fluids in periodic boxes, cell lists, a Langevin
  thermostat, virial pressure, and radial distribution functions.
- **Five examples**, each of which asserts its numbers and is run by CI. Give one a path and
  it writes an SVG; the plotting has no dependency. Two further examples are checks rather
  than showcases: `agents_quickstart` and `readme_check`.
- **`AGENTS.md` and `CLAUDE.md`** — the API on one page for a consumer, and the gate and
  conventions for a contributor. Written after an AI agent looked for pantometry on `PATH`, as a
  Python package, and in a consuming repository, found it in none of the three, and used
  MuJoCo instead. That is a distribution failure and not a documentation one, but the
  quickstart it now lands on is `examples/agents_quickstart.rs`, which CI runs, so it cannot
  drift from the library the way a hand-written snippet does.
- **CI** across Linux, macOS, Windows, two WebAssembly targets and Rust 1.78, with formatting,
  clippy at `-D warnings`, rustdoc at `-D warnings`, licence and advisory checks, and the
  examples.

### Fixed

- **The acoustic wall weighting.** A pressure sample sits *on* a wall, so it owns half a cell,
  and both `Tube` and `Room` divided its divergence by the whole `dx`. Every mode read low —
  1.4% on an 89-cell room — and the scheme converged at first order despite a second-order
  interior. Found by the *rate* of convergence rather than by its size.
- **`Tube`'s absorbing ends** consequently needed their own step limit, `Z·dx/2ρc²`. At the
  full CFL limit the corrected boundary inverts a wave instead of absorbing it: stable, silent
  and wrong.
- **`Bar1D`'s enthalpy reference.** Measured from absolute zero, the bar in these tests holds
  1.42 kJ, so a millijoule arriving is a change in the seventh significant figure and the
  audit's relative check asked for precision the arithmetic had thrown away — a floor that
  grew from 1.6e-12 J at 41 cells to 7.3e-12 J at 161. Measured from the initial temperature,
  the number being summed *is* the change.
- **A statistical test that passed on one seed.** The ideal-gas check asserted that doubling
  the density doubles the departure from `PV = Nk_BT`; across four seeds the ratio came out
  1.35, 1.61, 2.34 and 2.92, averaging to 2.06. It now averages over seeds, because the fix for
  a noisy statistical test is more samples and not a wider tolerance.
- **Licence texts** now ship inside each crate. The packages declared `MIT OR Apache-2.0` and
  contained neither, which `cargo package` does not warn about.
- **Five tolerances that were loose rather than wrong.** A tail-correction ratio checked
  against the `rc⁻³` power law inside 0.02, where the exact ratio differs from that law by
  0.0188 — a known discrepancy filling 94% of the budget. It is now checked against the closed
  form to 1e-12, with the power law asserted separately as the limit it actually is. The
  Langevin settling test averaged one seed and allowed 5%; the seed-to-seed spread is 0.96%,
  so it now averages four seeds and allows 2% — more samples, and a *tighter* tolerance. Two
  gradients expected to be zero were checked with a relative comparison against zero. A
  crystal's negative pressure was asserted as `pressure.max(0.0)`, which an identically zero
  virial satisfies; it is now a band that excludes zero. And the bouncing ball's energy
  handover compared the heat that crossed against `start_energy - 0.0f64.max(0.0)` — the final
  energy was never computed, so the equality its comment promised was a one-sided bound.
  `ContactSystem` now answers `as_any` so the equality can be checked; publishing 5% extra
  heat fails it.

### Documented

- Every public item, with `#![deny(missing_docs)]` in all eight crates so it stays that way.
  244 items had no doc comment; the concentration was in `pantometry-units`, which is the API
  nobody can avoid.
- A `compile_fail` doctest proving `Length + Time` does not build — the workspace's reason for
  existing, previously asserted only in prose.

[Unreleased]: https://github.com/YounghyeonPark/pantometry/compare/v0.22.0...HEAD
[0.22.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.22.0
[0.21.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.21.0
[0.20.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.20.0
[0.19.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.19.0
[0.18.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.18.0
[0.17.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.17.0
[0.16.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.16.0
[0.15.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.15.0
[0.14.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.14.0
[0.13.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.13.0
[0.12.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.12.0
[0.11.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.11.0
[0.10.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.10.0
[0.9.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.9.0
[0.8.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.8.0
[0.7.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.7.0
[0.6.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.6.0
[0.5.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.5.0
[0.4.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.4.0
[0.3.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.3.0
[0.2.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.2.0
[0.1.0]: https://github.com/YounghyeonPark/pantometry/releases/tag/v0.1.0
