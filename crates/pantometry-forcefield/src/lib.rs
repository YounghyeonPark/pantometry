//! pantometry-forcefield: a small molecule, read from the chemical dictionary and typed for a
//! force field.
//!
//! The atomistic counterpart to `pantometry-protein`'s coarse network: one body per atom,
//! hydrogens included, with the bonds the chemistry has rather than the ones a cutoff implies.
//! The aim, over the steps that follow this one, is molecular mechanics on drug-sized molecules —
//! an energy, its minimum, and which conformations are stable — with every term checkable.
//!
//! **This step computes UFF's whole energy, minimises it, moves the molecule in time, and
//! computes the free energy of decoupling a ligand from its surroundings.** What is here:
//!
//! - [`Component::from_ccd`] reads one entry of the wwPDB Chemical Component Dictionary, strictly:
//!   every malformed or unsupported input is a [`CcdError`] naming what it refused. See [`ccd`].
//! - [`uff::assign`] gives every atom its Universal Force Field type from its element, its bond
//!   orders and the dictionary's aromatic flags, and [`UffType::parameters`] gives each type's
//!   Table I row in SI. See [`uff`] for the citation, the table and the rules.
//! - [`ForceField`] holds all six of UFF's terms — harmonic bond stretch, angle bend, torsion,
//!   inversion, Lennard-Jones van der Waals, and Coulomb electrostatics with 1-2 and 1-3
//!   exclusions — and gives the energy, term by term, and the analytic force on every atom. See
//!   [`energy`] and [`angular`] for the equations, the page each is on, and the places this
//!   departs from the paper as printed: the sign of `r_EN`, the sign of the linear angle term,
//!   and a torsion switched off continuously rather than abruptly near a straight angle.
//! - [`ForceField::minimise`] relaxes a geometry: L-BFGS with an Armijo line search, so the
//!   energy never rises, to a stated force. [`DihedralRestraint`] holds a dihedral while the
//!   rest relaxes, which is how the paper's torsion barriers are compared. See [`minimise`].
//! - [`qeq`] computes charge-equilibration charges (Rappé and Goddard 1991), the ones UFF
//!   prescribes, at a stated geometry: the exact Slater Coulomb integrals, hydrogen's
//!   charge-dependent orbital iterated as the paper does, and each misprint and ambiguity in the
//!   paper settled by reproducing its tables. [`ForceField::with_qeq_charges`] puts them in the
//!   electrostatic term, fixed at that geometry; zero charges stay the default.
//! - [`solvation`] is implicit water: the generalized Born electrostatic solvation energy of
//!   Onufriev, Bashford and Case (2004), OBC II, with ε = 80 and optional salt, and its analytic
//!   force through the Born radii at fixed charges. [`ForceField::with_generalized_born`] adds it
//!   to the energy and the force; vacuum stays the default. A Shrake–Rupley surface gives the
//!   paper's nonpolar term beside it, as a number, not a force. The evaluation is the direct
//!   sum's to the bit, computing each pair integral once and, given the atoms that do not move
//!   ([`GeneralizedBorn::with_frozen`]), their descreening once, and calling the kernel's own
//!   `exp` and `ln` ([`pantometry_core::math`]): on 3b's 987-atom complex a step of dynamics costs
//!   18.2 ms, against the direct sum's 87 with the platform's, and 11.4 times vacuum's.
//! - [`System::from_pdb`] reads a protein chain and its ligand out of a PDB entry and makes them one
//!   [`Component`] the rest of this crate takes unchanged: hydrogens placed by superposing each
//!   residue's dictionary template on the crystal, bonds with orders (peptide bonds and disulfides
//!   included), formal charges at pH 7 by stated rules, and every atom marked [`Part::Protein`] or
//!   [`Part::Ligand`]. Strict, as the dictionary reader is: a missing heavy atom, a gap in a chain,
//!   a residue that is not the sequence's or an unnamed `HETATM` residue is a [`PdbError`] that
//!   names it. See [`pdb`].
//! - [`Binding::new`] cuts a pocket of whole residues out of a [`System`] around its ligand, gives
//!   the pocket and the ligand QEq charges separately, and computes `ΔE_bind = E(complex) −
//!   E(pocket) − E(ligand)` at fixed geometry — in vacuum exactly the protein–ligand cross terms,
//!   which is tested — with OBC II's polar desolvation and the buried area beside it, and
//!   minimises the ligand in the rigid pocket with every protein atom frozen
//!   ([`Minimiser::with_frozen`]). [`RigidMotion`] moves the ligand.
//!   [`Binding::relaxing_hydrogens`] relaxes each of the three systems' own hydrogens on its frozen
//!   heavy atoms — the complex, then the pocket alone and the ligand alone from where the complex
//!   left them — and every quantity is then taken with each system at its own positions, so the
//!   vacuum ΔE_bind is the cross terms plus a hydrogen reorganisation that is tested to be exactly
//!   that. See [`binding`].
//! - **The congener series** (`tests/the_congener_series.rs`): the nine T4 lysozyme L99A entries
//!   of Morton and Matthews (1995) — benzene, benzofuran, indene, isobutylbenzene, indole,
//!   n-butylbenzene, p-xylene, o-xylene and ethylbenzene — are built by the same pipeline with no
//!   change to it. Their binding energies are set beside the experimental ΔG° that Mobley et al.
//!   (2007) tabulate. The comparison is reported, not asserted.
//! - [`MolecularDynamics`] integrates Newton's equations for the atoms, with the elements'
//!   standard atomic weights as masses: **BAOAB Langevin dynamics** (Leimkuhler and Matthews
//!   2013) under a heat bath, velocity Verlet without one, frozen atoms held to the bit, every
//!   random kick keyed by `(seed, step, atom)` so that, at a fixed step, a run is the same however
//!   it is cut into calls — a `Simulation` whose frames are not whole multiples of the step takes
//!   other substeps, and is another run — and the
//!   bath's work booked so kinetic + potential − work is conserved to the integrator's error. See
//!   [`dynamics`].
//! - [`Complex`] puts a [`Binding`] in motion: a **mobile zone** of whole residues within a cutoff
//!   of the ligand, every atom of them free, and the rest of the binding's pocket a **frozen
//!   buffer** held to the bit, in vacuum or OBC II over every atom. A mobile residue bonded across
//!   the binding's cut is refused by name, because its backbone end would be free in vacuum.
//!   [`Complex::run`] takes a [`Frame`] of observables every so many steps — the ligand's RMSD from
//!   its crystal pose without superposition, the same blind to which atom is on which site, its
//!   turn about its crystal plane's normal, its centroid's displacement, its distance to the cavity centre, its heavy-atom contacts, the
//!   interaction energy, the mobile atoms' temperature and the books — into a [`Record`] that also
//!   keeps every atom's mean-square fluctuation. [`Estimate::of`] gives a mean and its standard
//!   error from the autocorrelation time, and [`complex::mean_square_displacement`] turns a
//!   B-factor into `3B/(8π²)`. **Benzene stays bound in T4 lysozyme L99A**, its centroid within
//!   1.16 Å of the crystal's in every run. **Its ring does not sit where the crystal puts it**: in
//!   vacuum it sits 20–30° turned from the crystal orientation, between the crystal's sites, and
//!   turns, which the crystal's B-factors exclude. That is a finding about the model. See
//!   [`complex`].
//! - [`Decoupling`] is the **alchemical Hamiltonian** of one group of atoms — a ligand — and the
//!   rest, with three couplings in a [`Lambda`]: the group–environment Coulomb pairs scaled
//!   linearly, its van der Waals pairs through Beutler's soft core (α = 0.5, p = 1, as Mobley,
//!   Chodera and Dill 2006 used; the form read in the GROMACS manual, the paper not opened), and a
//!   [`Boresch`] restraint switched on linearly. The group's own interactions are kept at every λ:
//!   decoupling, not annihilation. Energy, forces and `∂U/∂λ` analytic. **Under OBC II** the
//!   group's charges are `λ_e q` in every GB term and the descreening across the partition is
//!   scaled by `λ_v`, so the decoupled end state is the environment alone in solvent and the
//!   group alone in vacuum, exactly — the end state openmmtools' alchemical GB builds, reached by
//!   another path. See [`alchemy`].
//! - [`Boresch`] holds a ligand's six relative coordinates to three receptor atoms, chosen by a
//!   stated geometric rule, and gives the analytic free energy of releasing it to 1 M — Boresch et
//!   al. 2003's closed form, read secondarily in Clark et al. 2023 (eq 7), with its two
//!   approximations undone in closed form beside it. See [`boresch`].
//! - [`Windows`] runs BAOAB at each state of a λ schedule — deterministic per window, the same
//!   however it is cut into calls, resumable — and estimates the free energy by **thermodynamic
//!   integration** (trapezoid, or Simpson on a uniform line) and **Bennett's acceptance ratio**,
//!   solved to the last bit by Newton's method inside a bisection bracket. **Its uncertainty is
//!   the delta method's, window by window**, because summing the intervals' variances — the
//!   common practice — leaves out the covariance of two intervals that share a window: measured,
//!   that understated σ̂ 1.7-fold on a decoupled particle's MD windows (with their samples thinned
//!   as well) and 1.26-fold on exact samples along a five-state chain. See [`free_energy`].
//! - [`Molecule`] is the kernel [`Domain`]: the atoms as [`Bodies`] with their names and bonds,
//!   so the scene layer draws a ball-and-stick molecule without knowing what a molecule is, the
//!   energy terms and the force as readings, and **one minimiser iteration per step** by default,
//!   so a run is the molecule relaxing, frame by frame — or, once [`Molecule::thermalised`], **a
//!   time step**, so a run is the molecule moving thermally.
//!
//! # What it is checked against
//!
//! Facts about the molecule and the paper that this code did not produce. Aspirin (`AIN`) parses
//! to the atom count and composition its own `_chem_comp.formula` states; every hydrogen has one
//! bond; every heavy atom's bond orders sum to its valence, with an aromatic bond counted as one
//! and a half; the types of aspirin's thirteen heavy atoms are the ones a chemist reads off the
//! structure — its acid and ester oxygens resonant, as the paper types a conjugated ether oxygen
//! (see [`uff::assign`]); and the natural angles of `C_3`, `C_R` and `C_1` are the tetrahedral, trigonal and
//! linear angles of geometry. Every refusal has a test that feeds it the input it refuses.
//!
//! The energy terms against closed forms: a bond's energy is zero at its natural length and
//! `½ k Δr²` away from it; `C_3–C_3` is 1.514 Å; the paper's own worked Si–O correction is
//! reproduced; the van der Waals minimum is `−D_IJ` at `x_IJ` and crosses zero at `x_IJ / 2^(1/6)`;
//! two unit charges at 1 Å are 332.0637 kcal/mol; the exclusions leave exactly the 1-4 pair of a
//! four-atom chain. Every angle form is zero with zero slope at θ₀ and has curvature `K` there,
//! with `K` worked by hand for methane; ethane's nine torsions sum to a barrier of exactly
//! `V_C_3` = 2.119 kcal/mol; eq 17 by hand for ethylene and benzene; the inversion is zero when
//! planar and matches a hand value when not. The forces against central finite differences of
//! the energy, to a tolerance derived from the step; the energy against translation and rotation;
//! and the net force against zero.
//!
//! QEq against closed forms — every `ns` self-repulsion as an exact fraction, Roothaan's 1s–1s
//! `J(R)`, `J → 1/R`, very unequal exponents against independent references, the paper's
//! two-atom eq 18, both ends of the range clamping, the total charge, rigid-motion invariance —
//! and against its own paper: Table II's eighteen alkali halides, Table III's four hydrogens and
//! 36 Table IV charges at stated experimental geometries, each to its printed precision plus its
//! measured sensitivity to the geometry. See [`qeq`].
//!
//! Generalized Born against closed forms: eq 5's pair integral against its own quadrature, in
//! every branch; a lone ion's radius `ρ − 0.09 Å` and Born energy; two ions far apart as two Born
//! terms and a screened Coulomb, to a bound derived from the descreening; two at one point as one
//! of their summed charge; a charge inside a dielectric sphere, where HCT gives the Coulomb-field
//! radius exactly and Kirkwood's series shows it 5.73% long; eq 6 with eq 8's constants; the
//! paper's 30 Å bound; and the force against central differences. See [`solvation`]. And the
//! evaluation against the direct sum it replaced (`tests/the_generalized_born_is_its_direct_sum.rs`),
//! with no tolerance: the energy, every radius and every force the same bits on aspirin, a cluster
//! with a coincident pair and 181L's pocket, under six frozen masks, with the frozen atoms where
//! they were kept and moved by 0.1 Å and by one ulp, and a complex's run the direct evaluation's
//! run.
//!
//! A protein against its own entry: T4 lysozyme L99A with benzene (PDB 181L) has every heavy atom
//! its templates name at the file's coordinates bit for bit, its sequence is `SEQRES`'s with the
//! two residues `REMARK 465` lists as unlocated, each residue's formula after protonation is the
//! textbook residue's adjusted by the rule, every heavy atom's bond orders meet its valence with its
//! formal charge, and the total, +9, is `Arg + Lys − Asp − Glu` counted from `SEQRES` and the
//! salt-bridged His31 the rule makes +1. Every
//! placed hydrogen keeps its template's bond length exactly and its angles within a bound its own
//! fit earns, no hydrogen is within 1.5 Å of a heavy atom, every superposition is no worse than one
//! built by hand from three atoms, and every refusal is fed the input it refuses.
//!
//! The congeners against their own entries and 181L: each entry's `SEQRES` is 181L's and its
//! `SEQADV` the same three conflicts (C54T, C97A, L99A). Its `HETATM` residues are waters, two
//! chlorides, one HED and the ligand only. Every one builds with His31 +1 by the rule, a total
//! of +9 and 2604 protein atoms. Each ligand's dictionary entry is held to its molecule by its
//! bond graph, written by hand, which tells o-xylene from p-xylene where a formula cannot.
//!
//! The minimiser against geometry whose minimum is known exactly: a diatomic relaxes to eq 2's
//! natural length, water to its two natural lengths and θ₀ (no non-bonded pair is left in
//! either), and methane to the regular tetrahedron — each to within the displacement its final
//! force allows, `|F| / λ_min` with λ_min the measured smallest Hessian eigenvalue. Then against
//! the paper's own minimised numbers: the relaxed torsion barriers of its Table II and the
//! structures of its Figures 4–7, compared row by row, and asserted only where the comparison is
//! earned.
//!
//! The dynamics against closed forms (`tests/the_dynamics.rs`): velocity Verlet's largest energy
//! departure on a harmonic well is its shadow's `E₀ (ωh)²/4` to a part in 10⁶ at three steps; one
//! C–H bond, released along its axis, follows Verlet's exact discrete cosine at `√(k/μ)` to 5e-13
//! of its amplitude; BAOAB samples a harmonic well's positions at `k_BT/k` at ωh = 0.5, 1 and 1.9,
//! and its momenta at `k_BT (1 − (ωh)²/4)` at the step and `k_BT` after `O`, each within four
//! standard errors from the series' own autocorrelation time; a bath at 0 K is the paper's
//! Appendix worked by hand; the kicks are uncorrelated between atoms and components; thermalised
//! velocities are χ² on `3N − 6` degrees of freedom, in mean and variance; aspirin under a bath has
//! `⟨½mv²⟩ = (3/2) k_BT` per element, and `⟨KE⟩ = (3N_free/2) k_BT` with and without frozen atoms;
//! NVE aspirin's energy error falls as `h²` (and per step as `h³`); each BAOAB step on a harmonic
//! well moves the books by exactly the closed form `(a²/2)(P₁² − P₂² + Q₂² − Q²)`, to rounding;
//! every atomic weight is CIAAW's, and the masses sum to every dictionary entry's formula weight
//! within CIAAW's stated uncertainties; at a fixed step a run is the same bits however it is chunked; frozen atoms
//! do not move; and a moved start gives the moved trajectory to 8e-15 Å.
//!
//! The complex in motion against what must hold (`tests/the_complex_in_motion.rs`): the buffer
//! keeps its bits and each atom its element's mass; the degrees of freedom are `3 N_mobile`; the
//! potential's force on every mobile atom is the whole complex's bit for bit, in vacuum and in
//! OBC II; a run is its seed however it is cut into calls, and rebuilt from the file; the NVE energy
//! error of the mobile zone falls as `h²`, within a bound earned over sixteen starts; every frame is
//! recomputed from its positions, and every observable checked on a hand-built input — a translation's
//! RMSD is its length, a ring turned by `θ` has RMSD `2a sin(θ/2)`, site RMSD zero at 60° and an
//! in-plane turn of `θ`, a uniformly turning ring has site RMSD² `2a²(1 − 3/π)`, contacts count
//! strictly below the cutoff, `3B/(8π²)` on a typed B, and the autocorrelation time of an AR(1)
//! series is `(1 + φ)/(2(1 − φ))`, on one series and as a mean over 32; the frozen-only terms are
//! gone from the potential, to a traced allowance; every atom bonded across the cut is refused; a
//! cutoff's zone is exactly the pocket `Binding::new` cuts at the same radius, at debug speed on the
//! 3 Å binding and in release on 10 Å; and, in release, equipartition holds per element over the
//! mobile atoms.
//!
//! The free-energy engine against closed forms (`tests/the_free_energy_against_closed_forms.rs`,
//! `tests/benzene_decoupled_from_its_pocket.rs`): BAR's estimate and both its variances
//! calibrated over 1600 seeds of exact samples of two harmonic wells,
//! `ΔF = (3/2) k_BT ln(k₁/k₀)`, in bands of ±0.14 on the z-scores' variance; the delta method's
//! also on correlated chains, where Shirts's formula, which assumes independence, gives under a
//! third of the variance the estimates have, and along a chain of five states, where summing the
//! intervals' variances understates σ̂ by 1.26; TI's and BAR's σ̂ through [`Windows`] itself,
//! over 600 small campaigns on a schedule with a corner; TI and BAR from MD windows on the same
//! wells, each quadrature's bias computed from the exact integrand; Boresch's closed form against
//! a quadrature of `e^(−βU)` with `U` the restraint's own energy, and the excluded tails by
//! quadrature, to 4e-13 `k_BT`, and that integral's Jacobian against Cartesian grids to 1e-13;
//! the soft core equal to UFF's pair to the bit at
//! λ = 1, zero at λ = 0, its closed form at `r = 0`, and its derivatives against finite
//! differences; one Lennard-Jones particle decoupled from a fixed atom against its configurational
//! integral by quadrature, and coupled again from the other end; `∂U/∂λ_e` the coupled cross
//! Coulomb energy to the bit in 181L's pocket; and a campaign the same bits however it is cut up.
//! Under generalized Born (`tests/benzene_bound_in_generalized_born.rs`): fully coupled, the
//! decoupled solvation term is the force field's to the bit; with the charges off it carries none
//! of the group's, to the bit; fully decoupled it is the pocket alone in solvent plus benzene
//! alone in vacuum, the solvation term to the bit and the whole to its rounding; the frozen
//! buffer's kept terms change no bit; a Born ion decouples to `(q²/2R)(1 − 1/ε)` by TI and BAR
//! through [`Windows`]; and two ions follow eqs 2–8 written out at every `(λ_e, λ_v)`.
//!
//! # Determinism, and where it stops
//!
//! No clock, no randomness, no hash order and no threads, so a run repeats bit for bit on one
//! machine, at any optimisation level. **Not across platforms.** IEEE-754 defines `+ − × ÷` and
//! `sqrt` exactly, but this crate also calls the platform's `cos`, `sin`, `ln`, `asin` and, under a
//! restraint, `atan2` — on every torsion evaluation as well as when a [`ForceField`] is built.
//! Those are not correctly rounded and are not required to be the same function on two
//! machines, so no digest of a minimum is pinned; see [`minimise`]. [`qeq`]'s integrals call
//! `exp`, with the same consequence for a charge. [`solvation`]'s `exp` and `ln`, one of each in
//! every pair, are [`pantometry_core::math`]'s and the same everywhere; its `tanh`, one per atom,
//! is the platform's. **The rest of the crate's `exp` and `ln` are still the platform's**:
//! [`qeq`]'s integrals, the BAR solver in [`free_energy`], the Langevin coefficient in
//! [`dynamics`], the Boresch correction and UFF's bond-order term. None is where the time goes —
//! each runs once per system, per integrator or per sample — and switching them would make no
//! result of a run the same across platforms, since a minimum, a trajectory and a free energy all
//! pass through UFF's torsion `cos` as well, and QEq through `powi`, whose rounding Rust does not
//! specify. It would change their bits, and is left to be decided with those.
//! [`pdb`]'s superposition is arithmetic and `sqrt`, but turning a rotor calls `sin` and `cos`.
//!
//! # What is deliberately not in it
//!
//! - **Dynamics without constraints, and with nothing but a Langevin bath.** [`dynamics`] is BAOAB
//!   and velocity Verlet at a 0.5 fs step, with every X–H bond free: no SHAKE or RATTLE, no
//!   barostat, no multiple time steps, no periodic box. The measurement that chose the step is in
//!   [`dynamics`]; constraints would buy a step about four times longer and are not needed yet.
//! - **Free energies by TI and BAR, not MBAR, and no nonpolar term in them.** MBAR would need
//!   every sample's energy at every state; see [`free_energy`] for why BAR between neighbours is
//!   enough here. The solvent-accessible surface has no useful gradient, so it is in neither leg
//!   of a cycle and is estimated at the end states apart ([`alchemy`]). **Benzene's absolute
//!   binding free energy to T4 lysozyme L99A under OBC II comes out at +6.1 ± 0.3 kcal/mol, +4.6
//!   with that nonpolar estimate, against experiment's −5.19**: the model does not bind benzene,
//!   and `tests/benzene_bound_in_generalized_born.rs` says which of its parts are candidates. It
//!   is reported, not asserted.
//! - **No conformer search.** The minimiser finds the minimum downhill from where it starts; the
//!   torsion scans hold a dihedral to find a barrier, and nothing looks for the global minimum.
//! - **No reading of the paper but this crate's in use**: [`Variant`] keeps one other — `r_EN`
//!   added as eq 2 prints it, or dropped — only as the evidence that settled the sign, which the
//!   paper's own minimised structures did (see [`energy`]). The other question 1d left open, the
//!   group-6 sp³–sp² torsion's minimum, was settled by typing a conjugated O or S resonant as the
//!   paper does, after which the row it asked about is reached only by an oxonium oxygen.
//! - **No ring perception, no bond-order guessing, no protonation of a dictionary entry.** The
//!   dictionary states aromaticity, bond orders and explicit hydrogens, and a second opinion
//!   computed here could only disagree with it. A file without bonds is refused, not guessed at.
//!   The one place protonation is decided is [`pdb`], for the twenty amino acids at pH 7, by a
//!   table of stated rules applied to the dictionary's templates — no pKa is computed, a histidine
//!   in a salt bridge is +1 by a stated geometric rule and any other is neutral, and no
//!   hydrogen-bond network
//!   is optimised.
//! - **No elements beyond H, C, N, O, F, P, S, Cl, Br and I**, and no hypervalent sulfur types —
//!   [`uff::assign`] documents the sulfone gap.
//! - **No charges unless asked for, no solvent unless asked for.**
//!   Electrostatics is computed
//!   from charges a caller supplies, and they default to zero, because UFF's valence parameters
//!   were fitted without them (p. 10031 of the UFF paper). QEq's are an option
//!   ([`ForceField::with_qeq_charges`]), **fixed at the geometry they were computed for**: the
//!   force treats them as constants, as the QEq paper uses them, so a minimisation with them does
//!   not re-equilibrate as it moves. Not here: charges that follow the geometry, and QEq's
//!   polarisation extensions. A hypervalent sulfur is refused by name ([`Unsupported`]) rather
//!   than given the wrong radius.
//! - **Implicit solvent is generalized Born only, its electrostatic part only in the force
//!   field, and opt-in** ([`ForceField::with_generalized_born`]). The nonpolar `0.005 × SASA` term
//!   is computed ([`solvation::surface_area`], [`solvation::nonpolar_energy`]) but not added to
//!   the energy a minimiser sees, because a point-counted surface has no gradient worth the name.
//!   Not here: Poisson–Boltzmann, explicit water, GB with a solute dielectric other than 1, a
//!   and a cutoff.
//!   The radii and scale factors are AMBER's as OpenMM holds them, not read in primary, and Br and
//!   I are outside the set OBC was fitted with; see [`solvation`]. **With QEq charges, OBC II
//!   over-solvates** small molecules against experiment — mean −2.1, RMS 3.4 kcal/mol, ethers by
//!   about −4 to −6 — so absolute solvation and binding energies built on it are not
//!   quantitative.
//! - **A binding energy at fixed geometry, not a binding free energy** ([`binding`]): no entropy,
//!   no protein flexibility (the pocket is rigid at the crystal's coordinates), no polarisation or
//!   charge transfer between the partners (QEq on each separately), and a solvated estimate whose
//!   polar part converges with the pocket's size only slowly, because a cut pocket solvates its
//!   own cut surface — +7.3 kcal/mol at 6 Å, +14.9 at 20 Å and +15.0 for the whole protein, for
//!   benzene in T4 lysozyme L99A. GB solvates the apo cavity as water; [`ApoCavity::Empty`] is the
//!   reference for a cavity that is empty, as L99A's is. No comparison with an experimental affinity
//!   is asserted. **On the congener series it does not rank the ligands.** Across nine ligands
//!   whose measured ΔG° spans 2.1 kcal/mol, the minimised vacuum ΔE correlates with experiment at
//!   r = +0.46 at 6 Å, whose 95% interval is [−0.29, +0.86]. The solvated total correlates at
//!   r = −0.77, interval [−0.95, −0.22]: the wrong sign. The polar term grows with the ligand and
//!   with the cutoff, so the ligands that bind best are the ones it penalises most. **At the
//!   crystal pose the placed hydrogens decide the vacuum ΔE** for three of the nine: indene's
//!   H12 is 1.41 Å from Val111's HG13, which puts van der Waals at +224 kcal/mol. Neither
//!   placement sees the other partner's hydrogens, and [`Binding::relaxing_hydrogens`] is the
//!   fix: relaxed, no ligand hydrogen is within 1.94 Å of a protein hydrogen in any of the nine.
//!   It does not make the series rank: relaxed and minimised, the vacuum ΔE correlates at
//!   r = +0.48 at 6 Å, interval [−0.27, +0.87], and the solvated total at r = −0.81.
//!   **Neither does converging the polar term or emptying the cavity**: at 20 Å, where the
//!   polar term has converged, the solvated total correlates at r = −0.93, and at −0.87 with the
//!   cavity empty. The polar term is largest for the ligands that bind best, and its largest part
//!   for the two largest is the screened cross terms, which the empty cavity leaves in.
//!
//! Nothing here opens a file: every crate in this workspace compiles to `wasm32`, so a caller
//! reads the text and passes it in.

#![deny(missing_docs)]

pub mod alchemy;
pub mod angular;
pub mod binding;
pub mod boresch;
pub mod ccd;
pub mod complex;
pub mod dynamics;
pub mod energy;
pub mod free_energy;
pub mod minimise;
pub mod pdb;
pub mod qeq;
pub mod solvation;
pub mod uff;

pub use alchemy::{
    Alchemical, AlchemyError, AtLambda, Coupling, CrossPair, Decoupling, Lambda, SoftCore,
};
pub use angular::{Bend, Inversion, Torsion};
pub use binding::{
    ApoCavity, Binding, BindingError, Desolvation, HydrogenRelaxation, Interaction, RigidMotion,
};
pub use boresch::Boresch;
pub use ccd::{Atom, Bond, BondOrder, CcdError, Component, Coordinates, Element};
pub use complex::{Complex, ComplexError, Estimate, Frame, Record, Solvent};
pub use dynamics::{Bath, MolecularDynamics, Potential};
pub use energy::{Energy, Evaluation, ForceField, Unsupported, Variant};
pub use free_energy::{Bar, FreeEnergy, Protocol, Quadrature, Sample, Window, Windows};
pub use minimise::{DihedralRestraint, Minimiser, Progress, Status};
pub use pdb::{Histidine, Part, PdbError, Placement, Residue, Selection, System};
pub use qeq::{Charges, Qeq, QeqError};
pub use solvation::{DecoupledSolvation, GeneralizedBorn, Rescaling};
pub use uff::{Parameters, TableI, UffType};

use pantometry_core::integrator::substeps_for;
use pantometry_core::{Bodies, Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{LengthVec, Qty, Time};

/// A molecule as a domain: atoms at places, typed, bonded, and drawable.
///
/// # How it moves
///
/// **By default downhill, one minimiser iteration per step; in time once it is given dynamics.**
///
/// [`Molecule::thermalised`] (or [`Molecule::with_dynamics`]) puts it in **molecular dynamics**:
/// [`Domain::step`] then advances `dt` of time by BAOAB Langevin dynamics, or velocity Verlet with
/// [`Bath::Isolated`] (see [`dynamics`]), in steps of at most [`Molecule::time_step`] — 0.5 fs
/// by default, and [`Domain::max_stable_dt`] says so, so the kernel cuts a frame into steps of it.
/// The ledger is then `energy` = kinetic + potential − thermostat work, which moves only by the
/// integrator's error: `O(h²)`, which no [`Simulation`](pantometry_core::Simulation) default of
/// 1e-9 can hold, so a simulation running a molecule in motion sets
/// [`Molecule::ENERGY_TOLERANCE`] for energy, or refuses on the first frame — on purpose.
///
/// Minimising: [`Domain::step`] takes one L-BFGS step of
/// [`Minimiser`] on the positions, so each frame of a run is the molecule further relaxed and
/// the `energy` reading never rises from one frame to the next. It stops moving when the largest
/// force on any atom is at most the tolerance ([`Molecule::DEFAULT_TOLERANCE`], or
/// [`Molecule::with_tolerance`]) or when no step lowers the energy, and says which. The step's
/// `dt` is not used: this is a minimisation, not dynamics, so [`Domain::max_stable_dt`] is
/// infinite. The ledger stays empty — the energy that leaves is not a flow to anywhere, and
/// claiming one would be inventing a sink.
///
/// [`Molecule::minimise`] runs the same iterations to convergence in one call. On a molecule in
/// motion it moves the atoms and leaves the velocities, so the books jump by the energy it took
/// out: minimise first, then thermalise.
///
/// # Readings
///
/// `atoms`, `heavy atoms`, `bonds` and `formal charge`, and the energy at the current positions in
/// **kcal/mol** — the paper's unit and the one every reader of a force field compares in — as
/// `energy`, `bond stretch`, `angle bend`, `torsion`, `inversion`, `van der Waals` and
/// `electrostatic`. `energy` is UFF's total, the sum of the six. `max force` and `rms force`
/// (kcal mol⁻¹ Å⁻¹, per-atom force vectors), `converged` (1 or 0) and `minimiser steps`. For a
/// molecule [`ForceField::new`] refuses, the energy and force readings are `NaN` and
/// [`Molecule::force_field`] says why. The last four are its [`Domain::diagnostics`]: they
/// describe the minimisation, not the molecule.
///
/// In dynamics, `converged` and `minimiser steps` give way to `temperature` (°C, as every
/// temperature reading in this workspace is), `kinetic energy`, `thermostat work`, `conserved
/// energy` (kinetic + the potential's change since the start − thermostat work, kcal/mol: the
/// ledger's books) and `dynamics steps`; the diagnostics
/// are then `max force`, `rms force`, `conserved energy` and `dynamics steps`. The temperature is
/// read right after the bath's `O` under a Langevin bath, where BAOAB's momenta are canonical, and
/// at the whole step without one (see [`dynamics`]).
#[derive(Clone, Debug)]
pub struct Molecule {
    name: String,
    component: Component,
    types: Vec<UffType>,
    force_field: Result<ForceField, Unsupported>,
    at: Vec<[f64; 3]>,
    minimiser: Minimiser,
    dynamics: Option<MolecularDynamics>,
    time_step: f64,
    saved: Option<Saved>,
}

/// What [`Domain::checkpoint`] keeps: the positions, the minimiser and the dynamics.
type Saved = (Vec<[f64; 3]>, Minimiser, Option<MolecularDynamics>);

impl Molecule {
    /// The tolerance a new molecule minimises to: the largest force on any atom at most
    /// 1e-4 kcal mol⁻¹ Å⁻¹, in newtons. A bond of `k` ~ 700 kcal mol⁻¹ Å⁻² is then within 1.4e-7 Å
    /// of where the force would put it.
    pub const DEFAULT_TOLERANCE: f64 = 1e-4 * minimise::KCAL_PER_MOL_ANGSTROM;

    /// The longest dynamics step a molecule takes, seconds: 0.5 fs, chosen from the measurement in
    /// [`dynamics`]. Aspirin's fastest bond vibrates with a period near 11 fs; NVE aspirin at 300 K
    /// blows up at 3 fs, and at 0.5 fs its energy stays within 0.11 kcal/mol (a fifth of `k_BT`) of
    /// its start over a picosecond, where 1 fs allows 0.45. Under BAOAB the whole-step kinetic
    /// temperature of a mode is low by `(ωh)²/4`: 1.9% for a C–H stretch at 0.5 fs, 7.5% at 1.
    pub const DEFAULT_TIME_STEP: f64 = 0.5e-15;

    /// A relative tolerance for [`quantity::ENERGY`](pantometry_core::conserved::quantity::ENERGY)
    /// that a [`Simulation`](pantometry_core::Simulation) running a molecule in motion can be given
    /// ([`Simulation::conservation_tolerance_for`](pantometry_core::Simulation::conservation_tolerance_for)):
    /// how far one frame may move `kinetic + (potential − its start) − thermostat work`, relative to
    /// the largest of those three and of the totals, as the kernel's audit judges it.
    ///
    /// **Earned, at [`Molecule::DEFAULT_TIME_STEP`]** (`tests/a_molecule_in_motion.rs`,
    /// `the_energy_tolerance_measured`). Over 20 ps of aspirin at 300 K in a 1 ps⁻¹ bath, four
    /// seeds, the largest change was 2.63e-3 for a frame of one step, 3.20e-3 for 50 fs and
    /// 4.35e-3 for 1 ps. On the same runs a bath whose work was half counted would move them by at
    /// least 1.06e-2, 6.7e-2 and 0.146. 7e-3 sits between them, 1.6 times the worst correct frame
    /// and 1.5 times below the least wrong one. The band scales with the temperature as the
    /// books' scale does, so the ratio does not. It scales as the step squared, so a molecule
    /// given a longer step ([`Molecule::with_time_step`]) needs a tolerance larger by
    /// `(dt / 0.5 fs)²`, and with it loses that margin.
    ///
    /// Before the potential was counted from its start, the scale held UFF's arbitrary zero —
    /// 29.6 kcal/mol at aspirin's minimum. Half-counted work then moved a one-step frame by 5.0e-3
    /// of it, inside the 1e-2 this was set to, and passed.
    pub const ENERGY_TOLERANCE: f64 = 7e-3;

    /// A molecule named `name`, at the component's coordinates, typed by [`uff::assign`], with
    /// its [`ForceField`] built and every partial charge zero.
    ///
    /// Never fails: a molecule the force field refuses is still read, typed and drawn, and
    /// [`Molecule::force_field`] carries the refusal.
    pub fn new(name: impl Into<String>, component: Component) -> Molecule {
        let types = uff::assign(&component);
        let force_field = ForceField::new(&component, &types);
        let at = component.atoms().iter().map(|a| a.at).collect();
        Molecule {
            name: name.into(),
            component,
            types,
            force_field,
            at,
            minimiser: Minimiser::new(Molecule::DEFAULT_TOLERANCE),
            dynamics: None,
            time_step: Molecule::DEFAULT_TIME_STEP,
            saved: None,
        }
    }

    /// The same molecule minimising to `tolerance` newtons instead.
    pub fn with_tolerance(mut self, tolerance: f64) -> Molecule {
        self.minimiser = Minimiser::new(tolerance);
        self
    }

    /// Runs minimiser iterations until converged, stalled, or `max_steps` more have been taken,
    /// and says where it stopped; `None` for a molecule the force field refuses.
    pub fn minimise(&mut self, max_steps: usize) -> Option<Progress> {
        let ff = self.force_field.as_ref().ok()?;
        let mut last = None;
        for _ in 0..max_steps.max(1) {
            let p = self.minimiser.step(ff, &[], &mut self.at);
            last = Some(p);
            if p.status != Status::Running {
                break;
            }
        }
        last
    }

    /// The minimiser's status after the last step.
    pub fn status(&self) -> Status {
        self.minimiser.status()
    }

    /// The same molecule with partial charges `charges`, in elementary charges, one per atom.
    ///
    /// # Panics
    ///
    /// If `charges` is not one per atom.
    pub fn with_charges(mut self, charges: Vec<f64>) -> Molecule {
        assert_eq!(charges.len(), self.at.len(), "one charge per atom");
        if let Ok(ff) = self.force_field {
            self.force_field = Ok(ff.with_charges(charges));
        }
        // A different energy: the minimiser's history and status belong to the old one, and the
        // dynamics' forces were the old potential's — its cache knows the positions, not the
        // potential, so it is told.
        self.minimiser = Minimiser::new(self.minimiser.tolerance());
        if let Some(d) = &mut self.dynamics {
            d.forget_forces();
            if let Ok(ff) = &self.force_field {
                d.prepare(ff, &self.at);
            }
        }
        self
    }

    /// The same molecule **in molecular dynamics**: Maxwell–Boltzmann velocities at `temperature`
    /// kelvin drawn from `seed`, with the rigid motion removed (see [`MolecularDynamics::thermalised`]),
    /// the masses of its elements, and `bath`. From here [`Domain::step`] advances time. See
    /// [`Molecule`] and [`dynamics`].
    pub fn thermalised(self, temperature: f64, bath: Bath, seed: u64) -> Molecule {
        let dynamics =
            MolecularDynamics::for_elements(self.component.atoms().iter().map(|a| a.element))
                .with_bath(bath)
                .thermalised(&self.at, temperature, seed);
        self.with_dynamics(dynamics)
    }

    /// The same molecule moving under `dynamics` — built by hand, for frozen atoms or given
    /// velocities. Its forces are computed here, at the current positions.
    ///
    /// # Panics
    ///
    /// If `dynamics` is not one mass per atom.
    pub fn with_dynamics(mut self, mut dynamics: MolecularDynamics) -> Molecule {
        assert_eq!(dynamics.masses().len(), self.at.len(), "one mass per atom");
        if let Ok(ff) = &self.force_field {
            dynamics.prepare(ff, &self.at);
        }
        self.dynamics = Some(dynamics);
        self
    }

    /// The same molecule taking dynamics steps of at most `dt` seconds instead of
    /// [`Molecule::DEFAULT_TIME_STEP`].
    ///
    /// # Panics
    ///
    /// If `dt` is not positive and finite.
    pub fn with_time_step(mut self, dt: f64) -> Molecule {
        assert!(dt.is_finite() && dt > 0.0, "a time step must be positive");
        self.time_step = dt;
        self
    }

    /// The dynamics, if the molecule is in motion; `None` when it minimises.
    pub fn dynamics(&self) -> Option<&MolecularDynamics> {
        self.dynamics.as_ref()
    }

    /// The longest dynamics step it takes, seconds.
    pub fn time_step(&self) -> f64 {
        self.time_step
    }

    /// Its force field, or why there is none.
    pub fn force_field(&self) -> Result<&ForceField, &Unsupported> {
        self.force_field.as_ref()
    }

    /// The energy and forces at the current positions, or why there are none.
    pub fn evaluate(&self) -> Result<Evaluation, &Unsupported> {
        self.force_field().map(|ff| ff.evaluate(&self.at))
    }

    /// The component it was built from.
    pub fn component(&self) -> &Component {
        &self.component
    }

    /// Every atom's UFF type, in the component's atom order.
    pub fn types(&self) -> &[UffType] {
        &self.types
    }

    /// Where atom `i` is now, in metres.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn at(&self, i: usize) -> [f64; 3] {
        self.at[i]
    }
}

impl Bodies for Molecule {
    fn count(&self) -> usize {
        self.at.len()
    }

    fn position(&self, i: usize) -> LengthVec {
        LengthVec::new(
            Qty::from_si(self.at[i][0]),
            Qty::from_si(self.at[i][1]),
            Qty::from_si(self.at[i][2]),
        )
    }

    /// The atomic number, so a picture is coloured by element — the one per-atom quantity every
    /// reader of a ball-and-stick model already knows how to read.
    fn value(&self, i: usize) -> f64 {
        f64::from(self.component.atoms()[i].element.atomic_number())
    }

    /// Dimensionless: an atomic number is a count.
    fn value_unit(&self) -> &'static str {
        ""
    }

    /// The atom's name in the dictionary, `C7`, `HO1` — what the entry and every bond in it call
    /// the atom, and so what a reader can match against the file.
    fn label(&self, i: usize) -> Option<String> {
        self.component.atoms().get(i).map(|a| a.name.clone())
    }

    /// Every bond the dictionary states, as index pairs. Real by construction: the reader refuses
    /// a bond that names an atom the entry does not have.
    fn bonds(&self) -> Vec<[u32; 2]> {
        self.component
            .bonds()
            .iter()
            .map(|b| [b.atoms[0] as u32, b.atoms[1] as u32])
            .collect()
    }
}

impl Domain for Molecule {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// Minimising, infinite: a step is a minimiser iteration, and its length in time means
    /// nothing. In dynamics, [`Molecule::time_step`]. See [`Molecule`].
    fn max_stable_dt(&self, _now: Time) -> Time {
        if self.dynamics.is_some() {
            Qty::from_si(self.time_step)
        } else {
            Qty::from_si(f64::INFINITY)
        }
    }

    /// Minimising, one minimiser iteration. In dynamics, `dt` of time, in as many equal steps of
    /// at most [`Molecule::time_step`] as it takes ([`substeps_for`]): one, when the kernel has
    /// already cut the frame to [`Domain::max_stable_dt`]. Nothing, for a molecule the force field
    /// refuses. See [`Molecule`].
    fn step(&mut self, _t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let Ok(ff) = &self.force_field else {
            return Ok(());
        };
        match &mut self.dynamics {
            None => {
                self.minimiser.step(ff, &[], &mut self.at);
            }
            Some(d) => {
                if dt.to_si() > 0.0 {
                    let n = substeps_for(dt, Qty::from_si(self.time_step));
                    let h = dt.to_si() / f64::from(n);
                    for _ in 0..n {
                        d.step(ff, &mut self.at, h);
                    }
                }
            }
        }
        Ok(())
    }

    /// Minimising, empty: a minimisation lowers the energy without sending it anywhere, so there
    /// is no flow to record. In dynamics, `energy` = kinetic + potential − thermostat work
    /// ([`MolecularDynamics::ledger`]), which moves only by the integrator's error — and that is
    /// `O(h²)`, not the 1e-9 a [`Simulation`](pantometry_core::Simulation) checks by default, so a
    /// simulation running a molecule in motion has to set an energy tolerance; see [`Molecule`].
    fn ledger(&self) -> Ledger {
        match &self.dynamics {
            Some(d) => d.ledger(),
            None => Ledger::new(),
        }
    }

    /// True: it takes nothing from the bus, so its books are its own. Minimising, it holds
    /// nothing; in dynamics, its energy, which the audit holds to the tolerance it is given.
    fn books_balance(&self) -> bool {
        true
    }

    /// The positions, the minimiser's whole state and the dynamics' — velocities, step count,
    /// thermostat work — so a restored run repeats bit for bit, noise included.
    fn checkpoint(&mut self) {
        self.saved = Some((
            self.at.clone(),
            self.minimiser.clone(),
            self.dynamics.clone(),
        ));
    }

    fn restore(&mut self) {
        if let Some((at, minimiser, dynamics)) = &self.saved {
            self.at.clone_from(at);
            self.minimiser = minimiser.clone();
            self.dynamics = dynamics.clone();
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    fn readings(&self) -> Vec<Reading> {
        let atoms = self.component.atoms();
        let heavy = atoms.iter().filter(|a| a.element != Element::H).count();
        let charge: i32 = atoms.iter().map(|a| a.charge).sum();
        let ev = self.evaluate().ok();
        let e = ev.as_ref().map(|ev| ev.energy).unwrap_or(Energy::NAN);
        let kcal = |joules: f64| joules / uff::KCAL_PER_MOL;
        let (max_force, rms_force) = ev.as_ref().map_or((f64::NAN, f64::NAN), |ev| {
            let sq: Vec<f64> = ev
                .forces
                .iter()
                .map(|f| f[0] * f[0] + f[1] * f[1] + f[2] * f[2])
                .collect();
            let max = sq.iter().fold(0.0f64, |m, v| m.max(*v)).sqrt();
            let rms = (sq.iter().sum::<f64>() / sq.len().max(1) as f64).sqrt();
            (
                max / minimise::KCAL_PER_MOL_ANGSTROM,
                rms / minimise::KCAL_PER_MOL_ANGSTROM,
            )
        });
        let converged = if self.minimiser.status() == Status::Converged {
            1.0
        } else {
            0.0
        };
        let mut readings = vec![
            Reading::new(&self.name, "atoms", atoms.len() as f64, ""),
            Reading::new(&self.name, "heavy atoms", heavy as f64, ""),
            Reading::new(&self.name, "bonds", self.component.bonds().len() as f64, ""),
            Reading::new(&self.name, "formal charge", f64::from(charge), "e"),
            Reading::new(&self.name, "energy", kcal(e.total), "kcal/mol"),
            Reading::new(&self.name, "bond stretch", kcal(e.bond), "kcal/mol"),
            Reading::new(&self.name, "angle bend", kcal(e.angle), "kcal/mol"),
            Reading::new(&self.name, "torsion", kcal(e.torsion), "kcal/mol"),
            Reading::new(&self.name, "inversion", kcal(e.inversion), "kcal/mol"),
            Reading::new(
                &self.name,
                "van der Waals",
                kcal(e.van_der_waals),
                "kcal/mol",
            ),
            Reading::new(
                &self.name,
                "electrostatic",
                kcal(e.electrostatic),
                "kcal/mol",
            ),
            Reading::new(&self.name, "max force", max_force, "kcal/mol/Å"),
            Reading::new(&self.name, "rms force", rms_force, "kcal/mol/Å"),
        ];
        match &self.dynamics {
            None => {
                readings.push(Reading::new(&self.name, "converged", converged, ""));
                readings.push(Reading::new(
                    &self.name,
                    "minimiser steps",
                    self.minimiser.steps() as f64,
                    "",
                ));
            }
            Some(d) => {
                // The temperature the ensemble samples best: after the bath's `O` under BAOAB, the
                // whole step without a bath. See `dynamics`.
                let kelvin = if matches!(d.bath(), Bath::Langevin { .. }) && d.steps() > 0 {
                    d.half_step_temperature()
                } else {
                    d.temperature()
                };
                let kinetic = d.kinetic_energy();
                // The books exactly as the ledger keeps them, the potential counted from its start.
                let books = d
                    .ledger()
                    .get(pantometry_core::conserved::quantity::ENERGY)
                    .unwrap_or(f64::NAN);
                readings.extend([
                    Reading::new(&self.name, "temperature", kelvin - 273.15, "C"),
                    Reading::new(&self.name, "kinetic energy", kcal(kinetic), "kcal/mol"),
                    Reading::new(
                        &self.name,
                        "thermostat work",
                        kcal(d.thermostat_work()),
                        "kcal/mol",
                    ),
                    Reading::new(&self.name, "conserved energy", kcal(books), "kcal/mol"),
                    Reading::new(&self.name, "dynamics steps", d.steps() as f64, ""),
                ]);
            }
        }
        readings
    }

    /// Minimising, the four that describe the minimisation rather than the molecule: `max force`
    /// and `rms force` are its gradient, the residual a minimiser drives towards zero, and
    /// `converged` and `minimiser steps` are its state. A sweep that compared them across runs as
    /// though they converged to something would be measuring the stopping rule.
    ///
    /// In dynamics, `max force` and `rms force` — instantaneous now, and no more an answer than
    /// before — `conserved energy`, whose only movement is the integrator's error, and `dynamics
    /// steps`. The temperature, the kinetic energy and the thermostat work describe the molecule
    /// and its bath, and are not here.
    fn diagnostics(&self) -> &'static [&'static str] {
        if self.dynamics.is_some() {
            &[
                "max force",
                "rms force",
                "conserved energy",
                "dynamics steps",
            ]
        } else {
            &["max force", "rms force", "converged", "minimiser steps"]
        }
    }

    fn as_bodies(&self) -> Option<&dyn Bodies> {
        Some(self)
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}
