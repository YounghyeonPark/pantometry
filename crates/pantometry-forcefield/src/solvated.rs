//! A protein–ligand complex in a box of rigid TIP3P water, with counter-ions: step W4 of the
//! explicit-water track. W1 built the box and Ewald ([`crate::periodic`]), W2 the water
//! ([`crate::water`]), W3 benzene's hydration ([`crate::shake`], [`PeriodicDecoupling`]); this is
//! the complex leg's system.
//!
//! # What is built, and from what
//!
//! [`SolvatedComplex::new`] takes a [`Binding`] — the protein's residues within its cutoff and
//! the ligand, typed on the whole system, with QEq charges on the protein and on the ligand
//! separately ([`Binding::new`]) — and puts it, as it stands, in a box. **A binding whose cutoff
//! reaches every residue is the whole protein**, in the system's own atom order, and the same
//! code takes a pocket cut out of it as a small version for the tests. Nothing about the complex
//! is new here: its bonded terms are the vacuum force field's on whole molecules
//! ([`PeriodicForceField::whole`] walks the protein's bond graph, peptide bonds and disulfides
//! included, as it walks a water's), its 1-2 and 1-3 pairs are excluded from van der Waals and
//! Ewald alike, and 1-4 pairs are whole, as in vacuum. **No periodic version of a bonded term was
//! needed**: the walk makes the chain whole wherever the box cuts it, so every bend, torsion and
//! inversion reads the positions it reads in vacuum.
//!
//! # The box
//!
//! **Orthorhombic, aligned with the coordinates as given**, each edge the complex's extent along
//! that axis plus [`Solvation::margin`] on both sides, rounded up to whole lattice spacings: the
//! complex is at least `2 × margin` from each of its images. The default margin is 12 Å, the
//! cutoff (9 Å) and 3 Å more. Turning the complex to the box that holds it most tightly was
//! measured and not done: over 20 000 uniform rotations of 181L the smallest such box is 10.5% less
//! volume, where a truncated octahedron, which this crate does not have, would save about 30%.
//!
//! **The water** is [`WaterBox::lattice_box`] at [`Solvation::density`] — W3's 33.00 nm⁻³, where Yeh
//! and Hummer measure TIP3P's pressure under Ewald at −2.9 ± 2.5 bar — **less every water with an
//! atom within [`Solvation::clearance`] of a complex atom**, W3's 1.4 Å, **and every water whose
//! oxygen is within [`Solvation::oxygen_clearance`] of two heavy atoms of the complex on opposite
//! sides of it** — `(O − a)·(O − b) < 0`, the angle `a–O–b` past 90°, the oxygen inside the sphere
//! that has `ab` for its diameter — 2.6 Å being about the shortest O···N or O···O hydrogen bond;
//! every distance by minimum image ([`Solvation::removed`]).
//!
//! **The second rule is there because the first let a water into the protein.** With 1.4 Å alone,
//! 181L's box as phase 1 built it lost 780 waters, about the protein's volume — its mass, 18 419
//! g/mol, at the usual partial specific volume of a protein, 0.73 cm³/g (a typical value, not T4
//! lysozyme's measured one), is 22 330 Å³, 737 waters — and lattice water 4663 sat in a crevice
//! between Arg95 and Trp126, its oxygen within 2.6 Å of nine of their heavy atoms. Melted at 0.5 fs
//! with the protein frozen it went in further, its oxygen 1.77 Å from Arg95's Nε and a hydrogen
//! 0.83 Å from it: **a TIP3P hydrogen has no van der Waals term**, so nothing but its oxygen's keeps
//! it off an atom whose charge pulls it, and against a frozen protein the oxygen's was not enough.
//! The first 2 fs step after the melt put 415 kcal mol⁻¹ Å⁻¹ on that oxygen, the fourteenth turned
//! a water past what SETTLE can solve.
//!
//! **Why between two, and not near one.** Phase 1 removed every water whose oxygen was within 2.6 Å
//! of any heavy atom: 911 of 10 080, 174 more than the protein's volume, so the box was about 1.7%
//! under-dense, some −300 bar by TIP3P's compressibility (57.4 × 10⁻⁶ bar⁻¹, Izadi et al.). But a
//! water too close to one side of the protein is pushed straight back out by the melt, as at the
//! rest of the surface; **a water too close to two atoms on opposite sides of it has no direction
//! that relieves both**, and a frozen protein cannot open the gap. Lattice water 4663's widest pair
//! is at cos −0.947. The rule needs no tuning — 90° is where "between" begins — and is checked by
//! brute force over 27 images. On 181L with its hydrogens relaxed (the leg's box) the clearance
//! removes 788, the oxygen's rule 19 more, 807 in all, against phase 1's 909 for the same box; and
//! **water far from the protein is then within about 0.1–0.4% of the lattice's 33.00 nm⁻³,
//! depending on the shell**: 33.035 nm⁻³ beyond 10 Å, 33.122 beyond 6, over 1 ps after the leg's
//! preparation, twenty frames, while the potential energy was still falling
//! (`tests/benzene_bound_in_tip3p.rs`, `the_cost_of_the_mobile_leg_measured`). The shells differ
//! by more than each one's standard error, so that error is not the uncertainty of the density.
//! With the hydrogens relaxed, the clearance alone removes water 4663 too; with them as placed it
//! does not, and the oxygen's rule does, which a test holds. There is no barostat.
//!
//! # Counter-ions
//!
//! **The box is neutral**: one counter-ion of [`Solvation::counter_ion`]'s charge for each
//! elementary charge of the complex's total formal charge — nine chlorides for 181L, whose charge
//! is +9 ([`crate::System::from_pdb`]). A neutral box needs no correction for a net charge, and the
//! Ewald background, which is what makes a charged box's energy finite, then has nothing to
//! neutralise; QEq's charges sum to the formal charges, so the total is zero to the rounding of
//! the solve. **Each ion takes the place of a water**: of the waters left, the one whose oxygen is
//! farthest, by minimum image, from every complex atom and every ion already placed, the first on
//! a tie; the ion goes where that oxygen was, and the water is removed. Deterministic, and it puts
//! the ions in the bulk and apart. A complex whose charge an ion of this sign cannot cancel is
//! refused by name ([`SolvationError::NoCounterIon`]): this crate has a chloride and no cation,
//! because [`Element`] has no sodium.
//!
//! **The chloride is Joung and Cheatham's for TIP3P** (*J. Phys. Chem. B* **112**, 9020 (2008)),
//! **read secondarily**: σ = 4.477 656 957 Å and ε = 0.148 912 744 kJ/mol as OpenMM's
//! `amber14/tip3p.xml` (commit `7f4bb348`, which cites that paper) carries them, and as AMBER's
//! `frcmod.ionsjc_tip3p` prints them in its own form, `R_min/2` = 2.513 Å and ε = 0.035 591 0
//! kcal/mol (ambermini, commit `f7c421c4`), the two agreeing to the file's digits; the paper was
//! not opened. **A judgement, as water's own mixing is** ([`crate::water`]): they were fitted with
//! Lorentz–Berthelot combination against TIP3P, and here an ion combines with every other atom by
//! UFF's geometric rule, as a water oxygen does. Against the oxygen the contact distance is 4.216 Å
//! by the geometric rule and 4.281 Å by the arithmetic one, 1.5% closer; the depths are the same.
//!
//! # The potential and the constraints
//!
//! The force field is [`PeriodicForceField::solvated`] on the complex, the waters after it and the
//! ions after them ([`PeriodicForceField::with_ions`]), the complex's charges the binding's, and
//! the reciprocal sum through the particle mesh ([`PmeParameters::for_accuracy`] at
//! [`Solvation::accuracy`]). The dynamics ([`SolvatedComplex::dynamics`]) holds every water rigid
//! by SETTLE and **every bond to a hydrogen of the complex at UFF's natural length by SHAKE**
//! ([`Shake::to_hydrogen`]), which is what a 2 fs step needs. **Unlike benzene's six C–H, these
//! share atoms** — a methyl's three, an amine's — and are solved iteratively; and for bonds that
//! share an atom the mass-weighted Gram matrix of their gradients depends on the angle between
//! them, so the constrained distribution is not the stiff-bond one without a Fixman correction.
//! It is left out, as molecular-dynamics codes leave it out.
//!
//! # How much of the protein moves
//!
//! [`Flexibility`]: the protein frozen whole, residues beyond a radius of the ligand frozen
//! ([`Complex::zone`]'s rule: a residue moves when a heavy atom of it is within the radius of a
//! ligand atom at the crystal pose), or nothing frozen. The ligand, the water and the ions always
//! move. A frozen atom keeps its bits ([`MolecularDynamics::with_frozen`]); a held bond with one
//! end frozen holds the other at its length, and one frozen at both is skipped. Phase 1 measured
//! all three at the same cost a step; the complex leg moves the whole protein.
//!
//! # Preparing a leg
//!
//! [`SolvatedComplex::preparation`] is what a leg runs before its first window, as [`Stage`]s that
//! [`SolvatedComplex::prepare`] runs: **the melt** ([`SolvatedComplex::MELT`]), 0.5 ps at 0.5 fs
//! in a 50 ps⁻¹ bath and 1 ps at 2 fs in a 5 ps⁻¹ one with the protein frozen, which takes out the
//! lattice and the clearances' contacts; then, for a protein that moves, **the release**
//! ([`SolvatedComplex::release`]), 0.2 ps at 0.5 fs in a 50 ps⁻¹ bath and 0.5 ps at 2 fs in a
//! 5 ps⁻¹ one. **Released straight into 2 fs** from the frozen melt, residues within 8 Å of benzene
//! free drifted +149 kcal/mol/ps in NVE and within 10 Å SHAKE failed on a protein C–H: the crystal's
//! contacts — a backbone amide H 1.50 Å from an Asp Oδ1 among them, which the frozen protein carried
//! with 3 235 kcal mol⁻¹ Å⁻¹ on it — let go at the first step. Each stage that frees what the one
//! before held starts from Maxwell–Boltzmann velocities; the others carry theirs on.
//!
//! # Not here
//!
//! No barostat, no cation, no triclinic box, no rotation of the complex. The neighbour list is
//! opt-in ([`SolvatedComplex::with_neighbour_list`]), as it is for the force field.

use crate::binding::Binding;
use crate::ccd::{Element, ANGSTROM};
use crate::complex::Complex;
use crate::dynamics::{Bath, MolecularDynamics, Potential};
use crate::energy::Unsupported;
use crate::ewald::EwaldParameters;
use crate::free_energy::window_seed;
use crate::periodic::{PeriodicBox, PeriodicDecoupling, PeriodicForceField};
use crate::pme::PmeParameters;
use crate::shake::Shake;
use crate::uff::KCAL_PER_MOL;
use crate::water::{self, Settle, WaterBox};
use std::fmt;
use std::ops::Range;

/// Joung and Cheatham's TIP3P chloride σ, metres: 4.477 656 957 373 345 Å, as OpenMM's
/// `amber14/tip3p.xml` carries it. See the module documentation.
pub const CHLORIDE_SIGMA: f64 = 4.477_656_957_373_345e-10;

/// Joung and Cheatham's TIP3P chloride ε, kJ/mol: 0.148 912 744, 0.035 591 0 kcal/mol, as OpenMM's
/// `amber14/tip3p.xml` carries it.
pub const CHLORIDE_EPSILON_KJ: f64 = 0.148_912_744;

/// A monatomic ion: a point charge with a van der Waals term in UFF's form `D[(x/r)¹² − 2(x/r)⁶]`,
/// which is Lennard-Jones's with `x = 2^⅙ σ` and `D = ε`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ion {
    /// A name for reports, `Cl-`.
    pub name: &'static str,
    /// The charge, elementary charges.
    pub charge: f64,
    /// The mass, kilograms.
    pub mass: f64,
    /// UFF's `x`, the distance of the minimum, metres.
    pub vdw_distance: f64,
    /// UFF's `D`, the well depth, joules per molecule.
    pub vdw_energy: f64,
}

impl Ion {
    /// Chloride, −1, with Joung and Cheatham's TIP3P σ and ε ([`CHLORIDE_SIGMA`],
    /// [`CHLORIDE_EPSILON_KJ`]) and the crate's mass for chlorine.
    pub fn chloride() -> Ion {
        Ion {
            name: "Cl-",
            charge: -1.0,
            mass: Element::Cl.mass(),
            // 2^⅙ = √(∛2), by `+ − × ÷` and `sqrt`.
            vdw_distance: water::cube_root(2.0).sqrt() * CHLORIDE_SIGMA,
            vdw_energy: CHLORIDE_EPSILON_KJ / 4.184 * KCAL_PER_MOL,
        }
    }
}

/// How a complex is put in water: see the module documentation for each choice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Solvation {
    /// The least distance from the complex's extent to each face of the box, metres.
    pub margin: f64,
    /// The lattice water's density, kg m⁻³.
    pub density: f64,
    /// A water with an atom within this of a complex atom is removed, metres.
    pub clearance: f64,
    /// A water whose oxygen is within this of two heavy atoms of the complex on opposite sides of
    /// it — the angle at the oxygen between them more than 90° — is removed too, metres.
    pub oxygen_clearance: f64,
    /// The real-space cutoff, metres.
    pub cutoff: f64,
    /// The Ewald accuracy, in `k_e Q / r_c` ([`EwaldParameters::for_accuracy`]), and the mesh's.
    pub accuracy: f64,
    /// The seed of the waters' orientations.
    pub seed: u64,
    /// The counter-ion that neutralises the complex.
    pub counter_ion: Ion,
}

impl Solvation {
    /// W4's: a 12 Å margin, 33.00 waters per nm³ (W3's), a 1.4 Å clearance (W3's) and 2.6 Å from
    /// a water oxygen to two heavy atoms either side of it, `r_c` = 9 Å, δ = 10⁻⁶ through the
    /// mesh, and chloride.
    pub fn w4() -> Solvation {
        Solvation {
            margin: 12.0 * ANGSTROM,
            density: 33.0e27 * water::molecular_mass(),
            clearance: 1.4 * ANGSTROM,
            oxygen_clearance: 2.6 * ANGSTROM,
            cutoff: 9.0 * ANGSTROM,
            accuracy: 1e-6,
            seed: 0x4_4A7E,
            counter_ion: Ion::chloride(),
        }
    }
}

impl Solvation {
    /// The lattice a complex whose atoms are at `positions` (metres) is put in — waters along each
    /// axis, each edge the complex's extent and both margins rounded up to whole spacings — and the
    /// shift that moves the complex rigidly to the box's centre, metres. See the module
    /// documentation, "The box".
    ///
    /// # Panics
    ///
    /// If `positions` is empty.
    pub fn lattice_for(&self, positions: &[[f64; 3]]) -> ([usize; 3], [f64; 3]) {
        assert!(!positions.is_empty(), "a complex has atoms");
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in positions {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let spacing = water::cube_root(water::molecular_mass() / self.density);
        let lattice =
            [0, 1, 2].map(|k| ((hi[k] - lo[k] + 2.0 * self.margin) / spacing).ceil() as usize);
        let l = WaterBox::lattice_box(lattice, self.density, self.seed)
            .cell()
            .lengths();
        let shift = [0, 1, 2].map(|k| 0.5 * l[k] - 0.5 * (lo[k] + hi[k]));
        (lattice, shift)
    }

    /// One mark per water of `lattice`, `true` for each this solvation removes around a complex
    /// whose atoms are at `complex` (in the box, metres) and of `elements`, every distance by
    /// minimum image: **every water with an atom within [`Solvation::clearance`] of a complex
    /// atom, and every water whose oxygen is within [`Solvation::oxygen_clearance`] of two heavy
    /// atoms of the complex that lie on opposite sides of it** — `(O − a)·(O − b) < 0`, the angle
    /// `a–O–b` more than 90°, the oxygen inside the sphere on `ab` as its diameter. See the module
    /// documentation for why the second rule is that and not every heavy atom.
    ///
    /// # Panics
    ///
    /// If `elements` is not one per atom of `complex`.
    pub fn removed(
        &self,
        lattice: &WaterBox,
        complex: &[[f64; 3]],
        elements: &[Element],
    ) -> Vec<bool> {
        assert_eq!(complex.len(), elements.len(), "one element per atom");
        let cell = lattice.cell();
        let d =
            |a: [f64; 3], b: [f64; 3]| cell.minimum_image([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let c2 = self.clearance * self.clearance;
        let o2 = self.oxygen_clearance * self.oxygen_clearance;
        let heavy: Vec<[f64; 3]> = complex
            .iter()
            .zip(elements)
            .filter(|(_, &e)| e != Element::H)
            .map(|(p, _)| *p)
            .collect();
        let at = lattice.positions();
        let mut close = Vec::new();
        (0..lattice.count())
            .map(|w| {
                let atoms = &at[3 * w..3 * w + 3];
                if atoms.iter().any(|&p| {
                    complex.iter().any(|&c| {
                        let v = d(p, c);
                        dot(v, v) < c2
                    })
                }) {
                    return true;
                }
                close.clear();
                close.extend(
                    heavy
                        .iter()
                        .map(|&h| d(atoms[0], h))
                        .filter(|&v| dot(v, v) < o2),
                );
                close
                    .iter()
                    .enumerate()
                    .any(|(i, &a)| close[i + 1..].iter().any(|&b| dot(a, b) < 0.0))
            })
            .collect()
    }
}

/// Why a complex could not be solvated.
#[derive(Clone, Debug, PartialEq)]
pub enum SolvationError {
    /// The complex's formal charge is not a whole number of counter-ions of the opposite sign.
    NoCounterIon {
        /// The complex's total formal charge, e.
        charge: i32,
        /// The counter-ion's charge, e.
        ion: f64,
    },
    /// Fewer waters are left than ions are needed.
    NoRoomForIons {
        /// Ions needed.
        ions: usize,
        /// Waters left.
        waters: usize,
    },
    /// The complex has a term the force field does not support.
    Unsupported(Unsupported),
}

impl fmt::Display for SolvationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SolvationError::NoCounterIon { charge, ion } => write!(
                f,
                "a complex of charge {charge:+} cannot be neutralised by ions of charge {ion:+}"
            ),
            SolvationError::NoRoomForIons { ions, waters } => {
                write!(f, "{ions} ions are needed and {waters} waters are left")
            }
            SolvationError::Unsupported(u) => write!(f, "{u}"),
        }
    }
}

impl std::error::Error for SolvationError {}

/// How much of the protein moves. The ligand, the water and the ions always do.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Flexibility {
    /// Every protein atom frozen.
    Frozen,
    /// Every residue with a heavy atom within this many metres of a ligand atom at the crystal
    /// pose moves ([`Complex::zone`]); every other is frozen.
    Zone(f64),
    /// Nothing frozen.
    Mobile,
}

/// One stage of preparing a box for a leg ([`SolvatedComplex::prepare`]): `steps` steps of
/// `time_step` seconds in a Langevin bath of `friction` per second, with the protein held as
/// `flexibility` says.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stage {
    /// How much of the protein moves.
    pub flexibility: Flexibility,
    /// The time step, seconds.
    pub time_step: f64,
    /// How many steps.
    pub steps: usize,
    /// The bath's friction, per second.
    pub friction: f64,
}

impl Stage {
    /// `steps × time_step`, seconds.
    pub fn duration(&self) -> f64 {
        self.steps as f64 * self.time_step
    }
}

/// A complex in a box of rigid TIP3P water with its counter-ions: see the module documentation.
/// Atoms in the order the complex's ([`Binding`]'s complex order), the waters, oxygen first, then
/// the ions.
#[derive(Clone, Debug)]
pub struct SolvatedComplex {
    binding: Binding,
    field: PeriodicForceField,
    positions: Vec<[f64; 3]>,
    masses: Vec<f64>,
    complex: usize,
    waters: usize,
    ions: usize,
    lattice: [usize; 3],
    overlapping: usize,
    shake: Shake,
}

impl SolvatedComplex {
    /// `binding`, at its positions moved rigidly to the centre of the box, in water as
    /// `solvation` says. See the module documentation.
    ///
    /// # Errors
    ///
    /// [`SolvationError`], naming what failed.
    pub fn new(
        binding: &Binding,
        solvation: &Solvation,
    ) -> Result<SolvatedComplex, SolvationError> {
        let component = binding.component();
        let n_complex = binding.positions().len();
        let (pocket_charge, ligand_charge) = binding.formal_charges();
        let charge = pocket_charge + ligand_charge;
        let ion = solvation.counter_ion;
        let per_ion = -f64::from(charge) / ion.charge;
        let ions = if charge == 0 {
            0
        } else if per_ion > 0.0 && per_ion.fract() == 0.0 {
            per_ion as usize
        } else {
            return Err(SolvationError::NoCounterIon {
                charge,
                ion: ion.charge,
            });
        };

        // The box: the extent and the margin, in whole lattice spacings.
        let at = binding.positions();
        let (lattice, shift) = solvation.lattice_for(at);
        let full = WaterBox::lattice_box(lattice, solvation.density, solvation.seed);
        let cell = full.cell();
        let mut solute: Vec<[f64; 3]> = at
            .iter()
            .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
            .collect();
        // Every held bond at its length from the start, each hydrogen moved along its own bond and
        // every heavy atom where it was: a start off the constraints is a step in the energy at
        // the first step, which no step size removes.
        let parameters = EwaldParameters::for_accuracy(&cell, solvation.cutoff, solvation.accuracy);
        let alone = PeriodicForceField::new(component, binding.types(), cell, parameters)
            .map_err(SolvationError::Unsupported)?;
        let bonded = alone.bonded().expect("a complex has bonded terms");
        let shake = Shake::to_hydrogen(component, bonded);
        let complex_masses: Vec<f64> = binding.elements().iter().map(|e| e.mass()).collect();
        let fixed: Vec<bool> = binding
            .elements()
            .iter()
            .map(|&e| e != Element::H)
            .collect();
        let before = solute.clone();
        shake.constrain_positions(&complex_masses, &fixed, &before, &mut solute);
        let kept = full.without_waters(&solvation.removed(&full, &solute, binding.elements()));
        let overlapping = full.count() - kept.count();

        // The ions, each where the oxygen farthest from the complex and the ions placed was.
        if ions > kept.count() {
            return Err(SolvationError::NoRoomForIons {
                ions,
                waters: kept.count(),
            });
        }
        let oxygen = |w: usize| kept.positions()[3 * w];
        let distance2 = |a: [f64; 3], b: [f64; 3]| {
            let d = cell.minimum_image([a[0] - b[0], a[1] - b[1], a[2] - b[2]]);
            d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
        };
        let mut nearest: Vec<f64> = (0..kept.count())
            .map(|w| {
                solute
                    .iter()
                    .map(|&s| distance2(oxygen(w), s))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect();
        let mut taken = vec![false; kept.count()];
        let mut sites = Vec::with_capacity(ions);
        for _ in 0..ions {
            let mut best = None;
            for (w, &d) in nearest.iter().enumerate() {
                let better = match best {
                    None => true,
                    Some((_, b)) => d > b,
                };
                if !taken[w] && better {
                    best = Some((w, d));
                }
            }
            let (w, _) = best.expect("a water is left");
            taken[w] = true;
            let site = oxygen(w);
            sites.push(site);
            for (v, d) in nearest.iter_mut().enumerate() {
                *d = d.min(distance2(oxygen(v), site));
            }
        }
        let solvent = kept.without_waters(&taken);

        let n = n_complex + 3 * solvent.count() + ions;
        let mesh = PmeParameters::for_accuracy(&cell, &parameters, solvation.accuracy, n);
        let field = PeriodicForceField::solvated(
            component,
            binding.types(),
            solvent.count(),
            cell,
            parameters,
        )
        .expect("the complex was accepted alone")
        .with_solute_charges(binding.charges().to_vec())
        .with_ions(&vec![ion; ions])
        .with_mesh(mesh);
        let mut positions = solute;
        positions.extend_from_slice(solvent.positions());
        positions.extend_from_slice(&sites);
        let mut masses = complex_masses;
        masses.extend(solvent.masses());
        masses.resize(masses.len() + ions, ion.mass);
        Ok(SolvatedComplex {
            binding: binding.clone(),
            field,
            positions,
            masses,
            complex: n_complex,
            waters: solvent.count(),
            ions,
            lattice,
            overlapping,
            shake,
        })
    }

    /// The binding the complex came from, at the positions it was built with (not moved into the
    /// box).
    pub fn binding(&self) -> &Binding {
        &self.binding
    }

    /// The force field: the complex, the waters and the ions, the mesh in use.
    pub fn field(&self) -> &PeriodicForceField {
        &self.field
    }

    /// The box.
    pub fn cell(&self) -> PeriodicBox {
        self.field.cell()
    }

    /// Every atom's position as built, metres: the complex moved to the box's centre with each
    /// hydrogen moved along its bond to the length SHAKE holds it at ([`SolvatedComplex::shake`]),
    /// the lattice waters that were kept, the ions on the oxygens they replaced.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Every atom's mass, kilograms.
    pub fn masses(&self) -> &[f64] {
        &self.masses
    }

    /// How many atoms the complex has: the first ones.
    pub fn complex_atoms(&self) -> usize {
        self.complex
    }

    /// The ligand's atoms.
    pub fn ligand(&self) -> Range<usize> {
        self.binding.ligand_range()
    }

    /// One mark per atom, the ligand's set: the group a decoupling takes out.
    pub fn ligand_mask(&self) -> Vec<bool> {
        let ligand = self.ligand();
        (0..self.positions.len())
            .map(|i| ligand.contains(&i))
            .collect()
    }

    /// How many waters there are.
    pub fn waters(&self) -> usize {
        self.waters
    }

    /// How many counter-ions there are: the last atoms.
    pub fn ions(&self) -> usize {
        self.ions
    }

    /// The lattice's waters along each axis before any was removed.
    pub fn lattice(&self) -> [usize; 3] {
        self.lattice
    }

    /// How many lattice waters the complex overlapped and were removed; the counter-ions took the
    /// place of [`SolvatedComplex::ions`] more.
    pub fn overlapping(&self) -> usize {
        self.overlapping
    }

    /// Every bond to a hydrogen of the complex, at UFF's natural length.
    pub fn shake(&self) -> &Shake {
        &self.shake
    }

    /// One mark per atom, `true` for each that moves under `flexibility`.
    pub fn mobile(&self, flexibility: Flexibility) -> Vec<bool> {
        let n = self.positions.len();
        let ligand = self.ligand();
        let complex: Vec<bool> = match flexibility {
            Flexibility::Frozen => (0..self.complex).map(|i| ligand.contains(&i)).collect(),
            Flexibility::Zone(radius) => Complex::zone(&self.binding, radius),
            Flexibility::Mobile => vec![true; self.complex],
        };
        let mut mobile = complex;
        mobile.resize(n, true);
        mobile
    }

    /// Molecular dynamics of the box at rest, in no bath: every water rigid by SETTLE, every bond
    /// to a hydrogen of the complex held by SHAKE, and every atom `flexibility` does not move
    /// frozen.
    pub fn dynamics(&self, flexibility: Flexibility) -> MolecularDynamics {
        let frozen = self.mobile(flexibility).iter().map(|m| !m).collect();
        MolecularDynamics::new(self.masses.clone())
            .with_constraints(Settle::tip3p(self.field.rigid_waters().to_vec()))
            .with_bond_constraints(self.shake.clone())
            .with_frozen(frozen)
    }

    /// The same box with the force field's real-space pairs kept in a Verlet neighbour list of
    /// skin `skin` metres ([`PeriodicForceField::with_neighbour_list`]): the same bits at every
    /// evaluation, and [`SolvatedComplex::decoupling`] keeps a list of its own.
    ///
    /// # Panics
    ///
    /// If the skin is not positive and finite.
    pub fn with_neighbour_list(mut self, skin: f64) -> SolvatedComplex {
        self.field = self.field.with_neighbour_list(skin);
        self
    }

    /// **The melt**, phase 1's, with the protein frozen: 0.5 ps at 0.5 fs in a 50 ps⁻¹ bath —
    /// W2's, which takes out the lattice's strain and the clearance's contacts without a step a
    /// hydrogen could cross a neighbour in — then 1 ps at 2 fs in a 5 ps⁻¹ one.
    pub const MELT: [Stage; 2] = [
        Stage {
            flexibility: Flexibility::Frozen,
            time_step: 0.5e-15,
            steps: 1000,
            friction: 50e12,
        },
        Stage {
            flexibility: Flexibility::Frozen,
            time_step: 2e-15,
            steps: 500,
            friction: 5e12,
        },
    ];

    /// **The release** of a protein that moves after a frozen melt: 0.2 ps at 0.5 fs in a
    /// 50 ps⁻¹ bath, then 0.5 ps at 2 fs in a 5 ps⁻¹ one, at `flexibility`. See the module
    /// documentation for what it is for.
    pub fn release(flexibility: Flexibility) -> [Stage; 2] {
        [
            Stage {
                flexibility,
                time_step: 0.5e-15,
                steps: 400,
                friction: 50e12,
            },
            Stage {
                flexibility,
                time_step: 2e-15,
                steps: 250,
                friction: 5e12,
            },
        ]
    }

    /// **What a leg at `flexibility` runs before its first window**: [`SolvatedComplex::MELT`],
    /// and then, unless the protein is [`Flexibility::Frozen`], [`SolvatedComplex::release`] at
    /// `flexibility`.
    pub fn preparation(flexibility: Flexibility) -> Vec<Stage> {
        let mut stages = SolvatedComplex::MELT.to_vec();
        if flexibility != Flexibility::Frozen {
            stages.extend(SolvatedComplex::release(flexibility));
        }
        stages
    }

    /// The key mixed into a preparation's seed, so that its stages draw from other streams than
    /// the windows of a campaign seeded alike ([`SolvatedComplex::stage_seed`]).
    pub const PREPARATION_STREAM: u64 = 0x5D58_8B65_6C07_8965;

    /// Stage `k`'s seed in a preparation seeded `seed`:
    /// [`window_seed`]`(seed ^ PREPARATION_STREAM, k)`.
    pub fn stage_seed(seed: u64, k: usize) -> u64 {
        window_seed(seed ^ SolvatedComplex::PREPARATION_STREAM, k)
    }

    /// Runs `stages` in order on `potential` from `at`, in place, at `temperature` kelvin, each
    /// stage's dynamics [`SolvatedComplex::dynamics`] at its flexibility in its own bath, seeded
    /// [`SolvatedComplex::stage_seed`]`(seed, k)` for stage `k`. **The first stage, and any whose flexibility is
    /// not the one before it, starts from Maxwell–Boltzmann velocities** drawn from that seed —
    /// atoms just released have none — and every other carries the velocities on. `on_stage`
    /// is given each stage's index, the stage, its dynamics and the positions as it ended; the
    /// dynamics returned is the last stage's.
    ///
    /// # Panics
    ///
    /// If `stages` is empty, or `at` is not one position per atom.
    pub fn prepare(
        &self,
        potential: &(impl Potential + ?Sized),
        stages: &[Stage],
        temperature: f64,
        seed: u64,
        at: &mut [[f64; 3]],
        mut on_stage: impl FnMut(usize, &Stage, &MolecularDynamics, &[[f64; 3]]),
    ) -> MolecularDynamics {
        assert!(!stages.is_empty(), "a preparation has a stage");
        let mut last: Option<(Flexibility, MolecularDynamics)> = None;
        for (k, stage) in stages.iter().enumerate() {
            let s = SolvatedComplex::stage_seed(seed, k);
            let md = self.dynamics(stage.flexibility).with_bath(Bath::Langevin {
                temperature,
                friction: stage.friction,
                seed: s,
            });
            let mut md = match last {
                Some((f, ref before)) if f == stage.flexibility => {
                    md.with_velocities(before.velocities().to_vec())
                }
                _ => md.thermalised(at, temperature, s),
            };
            md.run(potential, at, stage.time_step, stage.steps);
            on_stage(k, stage, &md, at);
            last = Some((stage.flexibility, md));
        }
        last.expect("a stage ran").1
    }

    /// The ligand decoupled from everything else in the box, with no restraint.
    ///
    /// # Panics
    ///
    /// Never for a binding's complex: its ligand is bonded to nothing else.
    pub fn decoupling(&self) -> PeriodicDecoupling {
        PeriodicDecoupling::new(&self.field, &self.ligand_mask()).expect("the ligand decouples")
    }
}
