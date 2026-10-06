//! A ligand and its pocket in motion: a mobile zone of whole residues around the ligand, a frozen
//! buffer of protein around that, and molecular dynamics on the mobile atoms, recording whether
//! and how the ligand stays bound.
//!
//! [`Complex::new`] takes a [`Binding`] — the pocket cut at the buffer's radius, its QEq charges,
//! its force fields, and its positions, hydrogen-relaxed if the binding was — and marks a mobile
//! zone inside it. [`Complex::minimise`] relaxes the mobile zone, [`Complex::dynamics`] gives a
//! [`MolecularDynamics`] that holds the buffer to the bit, and [`Complex::run`] advances it,
//! writing a [`Frame`] of observables every so many steps into a [`Record`], which also keeps
//! every atom's mean-square fluctuation. [`Estimate::of`] turns a correlated series into a mean
//! and a standard error from its autocorrelation time.
//!
//! # The choices, and why
//!
//! **The mobile zone is whole residues**, by the pocket's own rule: a protein residue is mobile when
//! any of its heavy atoms is within the mobile cutoff of any ligand atom at the crystal pose, and
//! then every atom of it is free, hydrogens included. Whole residues so that no side chain is held
//! at one end and free at the other; the crystal pose so that the zone is chosen once and does not
//! follow the ligand. The ligand is always mobile.
//!
//! **The buffer is the rest of the binding's pocket, frozen.** A binding cut at a larger radius than
//! the mobile zone's gives a shell of protein that the mobile atoms feel — its van der Waals walls
//! and its charges — and that does not move: its atoms are frozen in [`MolecularDynamics`], so their
//! positions keep their bits for the whole run. **A mobile atom bonded across the cut is refused**
//! ([`ComplexError::MobileAtTheCut`]): a mobile residue whose peptide partner the buffer left out
//! would have a backbone end free in vacuum. Whole residues reaching out from a small mobile radius
//! can reach far, so the buffer's radius has to be chosen with that in mind, and this is how a
//! choice too small says so.
//!
//! **The frame is the frozen protein's.** Nothing is superposed: the buffer is where the crystal put
//! it, so the ligand's RMSD from its crystal pose is taken in the laboratory frame, and a ligand that
//! slid or turned in the cavity has moved by exactly that.
//!
//! **Charges: the binding's — QEq on the pocket and on the ligand separately, fixed**, at the
//! crystal positions, for every atom of the simulated system, mobile and frozen. The reason is
//! [`Binding`]'s: fixed-charge practice, in which the ligand's charges are its own, and a binding
//! energy and a trajectory that use the same charges are the same model. The buffer's charges are
//! the pocket's QEq at the buffer's cut, so they are the charges of a fragment, and they include its
//! cut surface; the mobile atoms see them through the vacuum Coulomb term or through generalized
//! Born. Fixed in time: QEq at each step would move the charges with the geometry, which the
//! crate's [`ForceField`] does not do anywhere.
//!
//! **The potential is the complex's force field less the terms among frozen atoms alone**
//! ([`Complex::potential`]): those are a constant while the buffer does not move, and the force on a
//! mobile atom is the whole complex's, term for term. **With [`Solvent::GeneralizedBorn`]**, OBC II
//! over **every** atom of the simulated system is added: generalized Born is many-body through the
//! Born radii, so a frozen atom's radius changes when a mobile one moves past it, and none of the
//! sum is a constant. What is constant is the descreening integral between two frozen atoms, and
//! the model is given the buffer as frozen ([`GeneralizedBorn::with_frozen`]) so that it computes
//! those once; the rest — the radii and the pairs — is the whole system's `O(N²)` every step, where
//! the vacuum step costs mobile × all. The result is the direct evaluation's to the bit. Vacuum is
//! the default.
//!
//! # Observables
//!
//! Each [`Frame`] records, at the positions of one step:
//!
//! - the ligand's heavy-atom RMSD from its crystal pose, without superposition ([`rmsd`]), and
//!   the same blind to which atom is on which site ([`site_rmsd`]) — a benzene turned by 60° in
//!   its plane has an RMSD of its radius and a site RMSD of zero;
//! - the ligand's turn about the normal of its crystal plane ([`in_plane_turn`]), which tells,
//!   where the site RMSD cannot, whether a ring sits on its sites or between them;
//! - the displacement of the ligand's heavy-atom centroid from its crystal position;
//! - the distance from that centroid to the **cavity centre** — the centroid of the cavity-lining
//!   atoms, the protein heavy atoms within [`LINING_DISTANCE`] of a crystal ligand heavy atom,
//!   at their positions now. It moves when the mobile side chains do, so it measures the ligand
//!   against the cavity rather than against the laboratory;
//! - the ligand–protein heavy-atom contacts below [`CONTACT_DISTANCE`] ([`contacts`]), counted as
//!   pairs;
//! - the vacuum interaction energy, the protein–ligand van der Waals and Coulomb cross terms
//!   ([`Binding::cross_terms_at`]) — under generalized Born too, since the solvated interaction
//!   does not decompose into pairs;
//! - the temperature of the mobile atoms, from the kinetic energy after the bath's `O` on `3 N_mobile`
//!   degrees of freedom ([`MolecularDynamics::half_step_temperature`]);
//! - the books, `kinetic + (potential − its start) − thermostat work`, whose change is the
//!   integrator's error.
//!
//! **Leaving** is a displacement of the centroid above a threshold the caller states: there is no
//! one right number, and a ligand of benzene's size rattling in a cavity of its own size moves
//! about an ångström. The tests state theirs.
//!
//! # Against experiment: B-factors
//!
//! A crystal's isotropic B-factor is `B = 8π² ⟨u_x²⟩`, with `u_x` an atom's displacement along one
//! axis, so the three-dimensional mean-square displacement is
//! `⟨|u|²⟩ = 3B/(8π²)` ([`mean_square_displacement`]). The simulation's is the mean-square
//! fluctuation about the atom's mean position, [`Record::mean_square_fluctuation`]. **They are not the
//! same quantity**, and the comparison means only what survives the differences: a B-factor holds
//! static disorder across the crystal's unit cells and the lattice's own motion as well as
//! thermal motion, the refinement's restraints and its errors, and was measured in a crystal at
//! its temperature, while the simulation has a frozen buffer that cannot move, no lattice and no
//! water. A simulated fluctuation well above the B-factor's is a ligand that moves more than the
//! crystal allows; one well below it is not a contradiction.

//!
//! # What it measures, on benzene in T4 lysozyme L99A (PDB 181L)
//!
//! `tests/the_complex_in_motion.rs`, release, one core of the machine it was written on. Each run
//! starts from the hydrogen-relaxed binding, minimises the mobile zone to
//! [`Binding::HYDROGEN_TOLERANCE`], and draws 300 K velocities. It runs at 0.5 fs in a 1 ps⁻¹ bath,
//! discards 10 ps (2 ps under OBC II), and records a frame every 10 fs.
//!
//! **The cost of a step**, and what it buys:
//!
//! | zone in binding | atoms (mobile) | ms/step | ns/day |
//! | --- | --- | --- | --- |
//! | 6 Å in 10 Å, vacuum | 987 (322) | 1.63 | 26.5 |
//! | 6 Å in 14 Å, vacuum | 1486 (322) | 2.56 | 16.9 |
//! | 8 Å in 14 Å, vacuum | 1486 (719) | 5.00 | 8.6 |
//! | 6 Å in 10 Å, OBC II | 987 (322) | 87.2 | 0.50 |
//! | 6 Å in 10 Å, OBC II, since 3c-2 | 987 (322) | 45.5 | 0.95 |
//! | 6 Å in 10 Å, OBC II, with the kernel's `exp` and `ln` | 987 (322) | 18.2 | 2.37 |
//!
//! Generalized Born over every atom cost 53 times the vacuum step, so 100 ps of it would have taken
//! nearly five hours. It was run for 8 ps, which is one autocorrelation time of the vacuum run's
//! slowest observable. Step 3c-2 halved it without changing a bit of the result — the radii and
//! pairs are the same sums, and the trajectory below is the same one — and it was still 28 times
//! the vacuum step. The kernel's own `exp` and `ln` brought it to 11.4 times (see
//! [`crate::solvation`], "The cost"), and **that did change bits**: a few forces by a few ulps,
//! which a run of dynamics amplifies, so the OBC II column below is 3b's trajectory with the
//! platform's functions and a run now is a different trajectory of the same model.
//!
//! **Benzene stays bound in every run.** Its centroid never moved more than 1.16 Å from the crystal
//! position, against a threshold of 3 Å. Means, with standard errors from each series'
//! autocorrelation time τ:
//!
//! | | 6 Å in 10 Å, 100 ps | 6 Å in 14 Å, 50 ps | 8 Å in 14 Å, 50 ps | OBC II, 8 ps |
//! | --- | --- | --- | --- | --- |
//! | RMSD (Å) | 2.03 ± 0.33 (τ 8 ps) | 2.14 ± 0.32 | 2.065 ± 0.007 | 0.79 ± 0.03 |
//! | site RMSD (Å) | 0.579 ± 0.014 | 0.593 ± 0.015 | 0.656 ± 0.003 | 0.604 ± 0.014 |
//! | centroid displacement (Å) | 0.500 ± 0.016 | 0.528 ± 0.013 | 0.557 ± 0.007 | 0.518 ± 0.026 |
//! | distance to cavity centre (Å) | 1.197 ± 0.019 | 1.234 ± 0.024 | 1.200 ± 0.009 | 1.289 ± 0.027 |
//! | contacts < 4 Å | 22.76 ± 0.10 | 22.87 ± 0.11 | 22.22 ± 0.12 | 21.78 ± 0.20 |
//! | interaction (kcal/mol) | −23.42 ± 0.03 | −23.58 ± 0.04 | −24.21 ± 0.05 | −22.85 ± 0.10 |
//! | mobile temperature (K) | 299.96 ± 1.37 | 301.7 ± 1.7 | 300.9 ± 0.9 | 299.6 ± 4.2 |
//!
//! **The RMSD is the ring turning in its plane, not the ligand leaving, and the ring does not sit
//! where the crystal puts it.** In both 6 Å runs the RMSD wanders between 0.14 and 2.91 Å, with
//! τ = 8 ps, while the centroid stays near 0.5 Å. A ring of radius `a` = 1.39 Å turned by 180° has
//! an RMSD of `2a` = 2.78 Å. [`Frame::turn`] measures the orientation directly. Folded into the
//! 60° sector it sits in, the turn has an RMS of 22.6°, 23.1° and 25.4° in the three vacuum runs,
//! against 17.3° for a ring that fills its sectors evenly. It **peaks 20–30° from the crystal
//! orientation**, half a sector away, between the crystal's sites. Only 13%, 10% and 1.2% of
//! frames lie within ±10° of the crystal orientation, where a uniform turn would put a third. The
//! ring also turns through sectors: −300° net over 100 ps in the 6 Å zone in 10 Å, −101° in 14 Å,
//! and −30° in the 8 Å zone, which stayed in two neighbouring sectors. So each atom's fluctuation
//! about its mean is 1.94 Å² in the 10 Å binding, about `a²` = 1.93, the value for a ring that
//! turns all the way round. **A 100 ps run holds about six of the slow RMSD's correlation
//! times**, so its ±0.33 is itself uncertain by roughly a factor of two. The site RMSD, the
//! centroid and the interaction decorrelate within a picosecond and are measured well.
//!
//! **The zone and the buffer move the numbers little.** A wider buffer around the same 6 Å zone
//! changes every mean by about its standard error, or less. The 8 Å zone binds benzene 0.8 kcal/mol
//! more strongly and holds it 0.08 Å further from its sites. Under OBC II the vacuum interaction
//! energy is 0.6 kcal/mol weaker and the ring did not turn in 8 ps. **What the solvated run can
//! say is limited**: one τ of the vacuum RMSD, and a model that over-solvates small molecules
//! with QEq charges ([`crate::solvation`]) and whose polar term A-2 found not to rank the congener
//! series. It is reported beside vacuum as a
//! check of how far the solvent moves the pose, not as the better answer.
//!
//! **The books** (`kinetic + ΔU − work`) moved between −0.85 and +1.06 kcal/mol over the four
//! runs, a random walk under the bath as `tests/the_dynamics.rs` explains, against a kinetic
//! energy of several hundred kcal/mol.
//!
//! **Against the B-factors: the simulated ring disagrees with the crystal.** Benzene's six carbons
//! in 181L have B = 20.05–26.80 Å², a mean `3B/(8π²)` of 0.894 Å², an RMS displacement of
//! 0.95 Å. That is a ring held at one orientation with its atoms smeared by about an ångström. The
//! simulation's ring sits half a sector from that orientation and turns, so its per-atom
//! fluctuation is 1.94 Å² in the 6 Å zone in 10 Å and 1.09 for the same zone in 14 Å. That is not
//! the density the crystal measured. **This is a finding about the model**: vacuum, the frozen
//! buffer and UFF's description of the cavity are each a candidate, and nothing here separates
//! them. No fluctuation with 60° jumps folded out is reported, because the ring does not jump
//! between the crystal's sites. It sits between them.
//!
//! **The site RMSD is not a fluctuation to set beside the B-factor.** It is a nearest-site
//! distance ([`site_rmsd`]) with a ceiling. A ring of the crystal's radius (1.367 Å) turning
//! uniformly reads `2a²(1 − 3/π)` = 0.168 Å², however freely it turns. A Gaussian with the
//! crystal's own `⟨|u|²⟩` reads back 0.828 of it, by Monte Carlo over the crystal ring. The runs'
//! 0.350, 0.365 and 0.433 Å² exceed the free-turning ceiling because the ring is offset: at 30°
//! from its sites each atom is `2a sin 15°` = 0.71 Å from the nearest, 0.50 Å². The earlier
//! reading of those numbers as "0.39–0.48 of the crystal's, so no looser than the crystal" was not
//! earned, and is withdrawn.
//!
//! **The mobile protein's heavy atoms fluctuate by a tenth of their B-factors** (0.08–0.12 Å²
//! against 0.97). Their fluctuations correlate with the B-factors at r = +0.67 atom by atom, 95%
//! interval [+0.57, +0.75], and +0.79 by residue, [+0.52, +0.92] (6 Å in 10 Å). The surface
//! residues Lys83, Lys85 and Asp89 move most in both. **What this can mean**: the simulation
//! ranks the same atoms as mobile as the crystal does. **What it cannot**: a B-factor holds the
//! crystal's static disorder and its lattice's motion, and here the buffer is frozen, there is no
//! lattice and no water, and the run is 100 ps. A protein fluctuation below the B-factor is
//! therefore expected, and its size is not a test of the force field. The ligand is the opposite
//! case: it moves more than the crystal allows.
//!
//! **n-butylbenzene in 186L**, the series' best binder, 50 ps in vacuum. Its 6 Å zone needs a 12 Å
//! binding (1388 atoms, 439 mobile): at 10 Å Phe153 is bonded across the cut and the zone is
//! refused. A step costs 3.09 ms (14.0 ns/day). It stays bound: the centroid moves at most 1.32 Å,
//! a mean of 0.764 ± 0.007, and the interaction is −40.74 ± 0.06 kcal/mol, of which van der Waals
//! is −38.79. Minimising the zone moved the ligand 0.99 Å from the crystal, and it stays there:
//! RMSD 1.190 ± 0.005 Å, site RMSD 0.962 ± 0.002. Each heavy atom fluctuates by 0.095 Å², against
//! a mean `3B/(8π²)` of 1.73. Two of its butyl carbons carry B = 100.00 Å² in the entry, which
//! reads as a refinement's ceiling rather than a measurement. The protein's fluctuations correlate
//! with their B-factors less well here: r = +0.35 by atom, [+0.22, +0.47], and +0.29 by residue,
//! [−0.11, +0.61].

// Three-vectors as `[f64; 3]`, indexed by component as the formulas are written.
#![allow(clippy::needless_range_loop)]

use crate::binding::Binding;
use crate::ccd::{Element, ANGSTROM};
use crate::dynamics::{Bath, MolecularDynamics};
use crate::energy::ForceField;
use crate::minimise::{Minimiser, Progress, Status};
use crate::solvation::GeneralizedBorn;
use pantometry_core::conserved::quantity;
use std::fmt;

/// The distance below which a ligand heavy atom and a protein heavy atom are in contact, metres:
/// 4 Å, a common heavy-atom contact cutoff, about the sum of two carbons' van der Waals radii
/// plus half an ångström.
pub const CONTACT_DISTANCE: f64 = 4.0 * ANGSTROM;

/// How near a crystal ligand heavy atom a protein heavy atom must be to line the cavity, metres:
/// 5 Å, the first shell of contacts and a little more, so that every wall of a small cavity has
/// atoms in the set.
pub const LINING_DISTANCE: f64 = 5.0 * ANGSTROM;

/// The solvent the complex moves in.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Solvent {
    /// None: the force field alone. The default.
    #[default]
    Vacuum,
    /// Implicit water: OBC II generalized Born over every atom of the simulated system, mobile and
    /// frozen, added to the energy and the force ([`GeneralizedBorn::new`]).
    GeneralizedBorn,
}

/// Why a mobile zone could not be marked.
#[derive(Clone, Debug, PartialEq)]
pub enum ComplexError {
    /// The mobile cutoff is not positive, or not below the binding's own: there would be no frozen
    /// buffer.
    NoBuffer {
        /// The mobile cutoff asked for, metres.
        mobile: f64,
        /// The binding's cutoff, metres.
        binding: f64,
    },
    /// A mobile residue is bonded to an atom the binding left out: its peptide partner is beyond
    /// the buffer, so a free backbone end would move in vacuum. A larger binding cutoff is the fix.
    MobileAtTheCut {
        /// The residue, `A:LYS135`.
        residue: String,
    },
}

impl fmt::Display for ComplexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ComplexError::NoBuffer { mobile, binding } => write!(
                f,
                "a mobile zone of {} Å leaves no frozen buffer in a binding cut at {} Å",
                mobile / ANGSTROM,
                binding / ANGSTROM
            ),
            ComplexError::MobileAtTheCut { residue } => write!(
                f,
                "mobile residue {residue} is bonded to an atom beyond the buffer; cut the binding \
                 wider"
            ),
        }
    }
}

impl std::error::Error for ComplexError {}

/// What one step of a trajectory looked like: see the module documentation, "Observables".
/// Lengths in metres, energies in joules per molecule, the temperature in kelvin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Frame {
    /// The dynamics' step count at this frame.
    pub step: u64,
    /// The ligand's heavy-atom RMSD from its crystal pose, without superposition.
    pub rmsd: f64,
    /// The ligand's heavy-atom RMS distance from the nearest crystal site of the same element
    /// ([`site_rmsd`]): blind to a permutation of equivalent atoms, so a benzene turned by 60° in
    /// its plane is at zero, as its electron density is. Nearest-site, so it has a ceiling: see
    /// [`site_rmsd`].
    pub site_rmsd: f64,
    /// How far the ligand has turned about the normal of its crystal heavy atoms' best-fit plane,
    /// radians in (−π, π] ([`in_plane_turn`]): for benzene, the ring's in-plane orientation, a
    /// multiple of 60° where it sits on its own sites.
    pub turn: f64,
    /// How far the ligand's heavy-atom centroid is from its crystal position.
    pub centroid_displacement: f64,
    /// How far the ligand's heavy-atom centroid is from the cavity centre now.
    pub cavity_distance: f64,
    /// Ligand–protein heavy-atom pairs closer than [`CONTACT_DISTANCE`].
    pub contacts: usize,
    /// The protein–ligand van der Waals cross terms.
    pub van_der_waals: f64,
    /// The protein–ligand Coulomb cross terms, ε = 1.
    pub electrostatic: f64,
    /// The mobile atoms' temperature after the bath's `O`; `NaN` without a bath or before a step.
    pub temperature: f64,
    /// `kinetic + (potential − its start) − thermostat work`.
    pub books: f64,
}

impl Frame {
    /// `van_der_waals + electrostatic`: the vacuum interaction energy.
    pub fn interaction(&self) -> f64 {
        self.van_der_waals + self.electrostatic
    }
}

/// A trajectory's frames, and every atom's fluctuation about its mean position over them.
///
/// A frame is taken after each step whose count is a multiple of [`Record::every`] — the
/// dynamics' own count, so a run cut into calls records the same frames as one call.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    every: u64,
    frames: Vec<Frame>,
    sum: Vec<[f64; 3]>,
    square: Vec<f64>,
}

impl Record {
    /// An empty record for `complex`, taking a frame every `every` steps.
    ///
    /// # Panics
    ///
    /// If `every` is zero.
    pub fn new(complex: &Complex, every: u64) -> Record {
        assert!(every > 0, "a frame every zero steps");
        let n = complex.at.len();
        Record {
            every,
            frames: Vec::new(),
            sum: vec![[0.0; 3]; n],
            square: vec![0.0; n],
        }
    }

    /// The interval between frames, in steps.
    pub fn every(&self) -> u64 {
        self.every
    }

    /// The frames, in order.
    pub fn frames(&self) -> &[Frame] {
        &self.frames
    }

    /// One field of every frame, in order: `record.series(|f| f.rmsd)`.
    pub fn series(&self, field: impl Fn(&Frame) -> f64) -> Vec<f64> {
        self.frames.iter().map(field).collect()
    }

    /// Atom `i`'s mean-square fluctuation about its mean position over the frames, `⟨|r − ⟨r⟩|²⟩`,
    /// m²: zero for a frozen atom, `NaN` with no frame. Accumulated as displacements from the
    /// crystal position, so the subtraction does not lose the fluctuation to the size of the
    /// coordinates.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the atom count.
    pub fn mean_square_fluctuation(&self, i: usize) -> f64 {
        let n = self.frames.len() as f64;
        let s = self.sum[i];
        let mean_sq = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]) / (n * n);
        (self.square[i] / n - mean_sq).max(0.0)
    }

    fn take(&mut self, complex: &Complex, frame: Frame) {
        for (i, (p, c)) in complex.at.iter().zip(&complex.crystal).enumerate() {
            let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            for k in 0..3 {
                self.sum[i][k] += d[k];
            }
            self.square[i] += d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        }
        self.frames.push(frame);
    }
}

/// A mean and its standard error from a correlated series: see [`Estimate::of`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Estimate {
    /// The mean.
    pub mean: f64,
    /// The standard error of the mean, `√(2 τ_int Var / N)`.
    pub error: f64,
    /// The integrated autocorrelation time, in samples: ½ for independent samples.
    pub tau: f64,
    /// The number of samples.
    pub samples: usize,
}

impl Estimate {
    /// The mean of `xs`, and its standard error from the integrated autocorrelation time
    /// `τ_int = ½ + Σ ρ(t)`, summed over a self-consistent window — Sokal's rule, stopping at the
    /// first `t ≥ 6 τ_int(t)` (A. D. Sokal, "Monte Carlo methods in statistical mechanics",
    /// Cargèse lectures 1996, §3) — the same estimator `tests/the_dynamics.rs` earns its bounds
    /// with. For a series with no variance the error is zero and `τ_int` ½.
    ///
    /// # Panics
    ///
    /// If `xs` is empty.
    pub fn of(xs: &[f64]) -> Estimate {
        assert!(!xs.is_empty(), "an estimate of nothing");
        let n = xs.len();
        let mean = xs.iter().sum::<f64>() / n as f64;
        let var = xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / n as f64;
        let mut tau = 0.5;
        if var > 0.0 {
            for t in 1..n / 2 {
                let c: f64 = (0..n - t)
                    .map(|i| (xs[i] - mean) * (xs[i + t] - mean))
                    .sum::<f64>()
                    / ((n - t) as f64 * var);
                tau += c;
                if t as f64 >= 6.0 * tau {
                    break;
                }
            }
        }
        Estimate {
            mean,
            error: (2.0 * tau.max(0.5) * var / n as f64).sqrt(),
            tau,
            samples: n,
        }
    }
}

/// `⟨|u|²⟩ = 3B/(8π²)`: the three-dimensional mean-square displacement an isotropic B-factor `b`
/// stands for, in the units of `b` (m² in, m² out). `B = 8π² ⟨u_x²⟩` along one axis, and an
/// isotropic displacement has three equal axes. See the module documentation for what it can
/// and cannot be compared with.
pub fn mean_square_displacement(b: f64) -> f64 {
    3.0 * b / (8.0 * std::f64::consts::PI * std::f64::consts::PI)
}

/// The unweighted centroid of `at`.
///
/// # Panics
///
/// If `at` is empty.
pub fn centroid(at: &[[f64; 3]]) -> [f64; 3] {
    assert!(!at.is_empty(), "the centroid of nothing");
    let mut c = [0.0; 3];
    for p in at {
        for k in 0..3 {
            c[k] += p[k];
        }
    }
    let n = at.len() as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

/// The root-mean-square distance between `a[i]` and `b[i]`, **without superposition**: a
/// translation by `t` gives `|t|` exactly, up to rounding.
///
/// # Panics
///
/// If the two are not the same length, or are empty.
pub fn rmsd(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    assert_eq!(a.len(), b.len(), "one position against one");
    assert!(!a.is_empty(), "the RMSD of nothing");
    let sum: f64 = a.iter().zip(b).map(|(p, q)| distance_squared(*p, *q)).sum();
    (sum / a.len() as f64).sqrt()
}

/// The root-mean-square, over the atoms of `at`, of each one's distance to the nearest of `sites`
/// of the same kind — `kinds[i]` is atom `i`'s and site `i`'s, the sites being the same atoms
/// somewhere else. **Blind to which atom sits on which site**: a ring turned onto itself is at
/// zero, where [`rmsd`] gives the ring's radius.
///
/// **Nearest, not one-to-one, and so it has a ceiling.** Two atoms on one site both read zero.
/// A ring of radius `a` turning uniformly in its plane reads `⟨site RMSD²⟩ = 2a²(1 − 3/π)`,
/// 0.174 Å² at 1.39 Å, however freely it turns. A large displacement is read as the distance to a
/// neighbour's site: an isotropic Gaussian with benzene's crystal `⟨|u|²⟩` of 0.894 Å² reads back
/// 0.83 of itself. **It is not a mean-square fluctuation to set beside a B-factor** without those
/// two numbers beside it.
///
/// # Panics
///
/// If the three are not the same length, or are empty.
pub fn site_rmsd(at: &[[f64; 3]], sites: &[[f64; 3]], kinds: &[Element]) -> f64 {
    assert_eq!(at.len(), sites.len(), "one site per atom");
    assert_eq!(at.len(), kinds.len(), "one kind per atom");
    assert!(!at.is_empty(), "the RMSD of nothing");
    let sum: f64 = at
        .iter()
        .zip(kinds)
        .map(|(p, k)| {
            sites
                .iter()
                .zip(kinds)
                .filter(|(_, j)| *j == k)
                .map(|(q, _)| distance_squared(*p, *q))
                .fold(f64::INFINITY, f64::min)
        })
        .sum();
    (sum / at.len() as f64).sqrt()
}

/// How far `at` has turned from `crystal` (the same atoms, in the same order) about the normal of
/// `crystal`'s best-fit plane, radians in (−π, π], right-handed about that normal as the plane's
/// smallest principal axis orients it. Each atom's offset from its own centroid is projected onto
/// the plane in both, and the angles between the pairs are averaged on the circle
/// (`atan2(Σ sin, Σ cos)`), so a rigid turn by `θ` about the normal gives `θ`, whatever the
/// translation. For a ring this is its in-plane orientation; for a ligand that is not flat, the
/// turn of its heavy atoms about the axis along which they spread least.
///
/// # Panics
///
/// If the two are not the same length, or have fewer than three atoms.
pub fn in_plane_turn(at: &[[f64; 3]], crystal: &[[f64; 3]]) -> f64 {
    assert_eq!(at.len(), crystal.len(), "one position against one");
    assert!(at.len() >= 3, "a plane needs three atoms");
    let (c, c0) = (centroid(at), centroid(crystal));
    let mut cov = [[0.0; 3]; 3];
    for p in crystal {
        let d = [p[0] - c0[0], p[1] - c0[1], p[2] - c0[2]];
        for a in 0..3 {
            for b in 0..3 {
                cov[a][b] += d[a] * d[b];
            }
        }
    }
    let (values, vectors) = crate::dynamics::symmetric_eigen(cov);
    let k = (0..3)
        .min_by(|&i, &j| values[i].total_cmp(&values[j]))
        .expect("three");
    let n = [vectors[0][k], vectors[1][k], vectors[2][k]];
    let project = |v: [f64; 3]| {
        let along = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
        [
            v[0] - along * n[0],
            v[1] - along * n[1],
            v[2] - along * n[2],
        ]
    };
    let (mut sin, mut cos) = (0.0, 0.0);
    for (p, q) in at.iter().zip(crystal) {
        let u = project([p[0] - c[0], p[1] - c[1], p[2] - c[2]]);
        let u0 = project([q[0] - c0[0], q[1] - c0[1], q[2] - c0[2]]);
        let x = [
            u0[1] * u[2] - u0[2] * u[1],
            u0[2] * u[0] - u0[0] * u[2],
            u0[0] * u[1] - u0[1] * u[0],
        ];
        let (s, co) = (
            x[0] * n[0] + x[1] * n[1] + x[2] * n[2],
            u0[0] * u[0] + u0[1] * u[1] + u0[2] * u[2],
        );
        let r = (s * s + co * co).sqrt();
        if r > 0.0 {
            sin += s / r;
            cos += co / r;
        }
    }
    sin.atan2(cos)
}

/// How many pairs, one from `a` and one from `b`, are closer than `cutoff`.
pub fn contacts(a: &[[f64; 3]], b: &[[f64; 3]], cutoff: f64) -> usize {
    let c2 = cutoff * cutoff;
    a.iter()
        .map(|p| b.iter().filter(|q| distance_squared(*p, **q) < c2).count())
        .sum()
}

fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    distance_squared(a, b).sqrt()
}

/// A binding with a mobile zone: see the module documentation.
///
/// Atoms are in the binding's complex order, the pocket's then the ligand's.
#[derive(Clone, Debug)]
pub struct Complex {
    binding: Binding,
    mobile_cutoff: Option<f64>,
    mobile: Vec<bool>,
    mobile_residues: Vec<String>,
    solvent: Solvent,
    potential: ForceField,
    ligand_heavy: Vec<usize>,
    protein_heavy: Vec<usize>,
    lining: Vec<usize>,
    crystal: Vec<[f64; 3]>,
    at: Vec<[f64; 3]>,
}

impl Complex {
    /// `binding` with every residue that has a heavy atom within `mobile_cutoff` metres of a ligand
    /// atom at the crystal pose mobile, the ligand mobile, and the rest of the binding's pocket a
    /// frozen buffer; at the binding's positions — its relaxed complex's, for a binding whose
    /// hydrogens were relaxed — in `solvent`.
    ///
    /// # Errors
    ///
    /// [`ComplexError::NoBuffer`] for a mobile cutoff that is not positive or not below the
    /// binding's, and [`ComplexError::MobileAtTheCut`] for a mobile residue bonded across the
    /// binding's cut.
    pub fn new(
        binding: &Binding,
        mobile_cutoff: f64,
        solvent: Solvent,
    ) -> Result<Complex, ComplexError> {
        if !(mobile_cutoff > 0.0 && mobile_cutoff < binding.cutoff()) {
            return Err(ComplexError::NoBuffer {
                mobile: mobile_cutoff,
                binding: binding.cutoff(),
            });
        }
        let mobile = Complex::zone(binding, mobile_cutoff);
        let mut complex = Complex::with_mobile(binding, mobile, solvent)?;
        complex.mobile_cutoff = Some(mobile_cutoff);
        Ok(complex)
    }

    /// The mask [`Complex::new`] frees, complex order: the ligand, and every atom of each pocket
    /// residue with a heavy atom within `cutoff` metres of a ligand atom — hydrogens included on
    /// the ligand's side — at the crystal pose. The pocket's own rule ([`Binding::new`]), so a
    /// zone at a radius is the pocket cut at that radius. Nothing is checked: a zone across the cut
    /// is refused by [`Complex::with_mobile`], not here.
    pub fn zone(binding: &Binding, cutoff: f64) -> Vec<bool> {
        let elements = binding.elements();
        let crystal = binding.crystal_positions();
        let ligand = binding.ligand_range();
        let mut mobile = vec![false; crystal.len()];
        for k in ligand.clone() {
            mobile[k] = true;
        }
        for span in binding.residue_spans() {
            if span.clone().any(|k| {
                elements[k] != Element::H
                    && crystal[ligand.clone()]
                        .iter()
                        .any(|&l| distance(crystal[k], l) < cutoff)
            }) {
                for k in span.clone() {
                    mobile[k] = true;
                }
            }
        }
        mobile
    }

    /// `binding` with the atoms `mobile` marks (complex order) mobile and every other atom frozen,
    /// in `solvent`: any zone, for a caller whose zone is not a cutoff's — whole residues or not,
    /// with the ligand or without. [`Complex::new`] is this with the cutoff's mask.
    ///
    /// # Errors
    ///
    /// [`ComplexError::MobileAtTheCut`] when a marked atom is bonded across the binding's cut.
    ///
    /// # Panics
    ///
    /// If `mobile` is not one entry per atom.
    pub fn with_mobile(
        binding: &Binding,
        mobile: Vec<bool>,
        solvent: Solvent,
    ) -> Result<Complex, ComplexError> {
        let elements = binding.elements();
        let crystal = binding.crystal_positions().to_vec();
        assert_eq!(mobile.len(), crystal.len(), "one mark per atom");
        let ligand = binding.ligand_range();
        let mut mobile_residues = Vec::new();
        for (span, label) in binding.residue_spans().iter().zip(binding.residues()) {
            if span.clone().any(|k| mobile[k]) {
                if span
                    .clone()
                    .any(|k| mobile[k] && binding.bonded_outside()[k])
                {
                    return Err(ComplexError::MobileAtTheCut {
                        residue: label.clone(),
                    });
                }
                mobile_residues.push(label.clone());
            }
        }
        let mut potential = binding.force_field().touching(&mobile);
        if solvent == Solvent::GeneralizedBorn {
            let frozen: Vec<bool> = mobile.iter().map(|m| !m).collect();
            potential = potential.with_solvation(
                GeneralizedBorn::new(elements).with_frozen(&frozen, binding.positions()),
            );
        }
        let heavy = |k: &usize| elements[*k] != Element::H;
        let ligand_heavy: Vec<usize> = ligand.clone().filter(heavy).collect();
        let protein_heavy: Vec<usize> = (0..binding.pocket_len()).filter(heavy).collect();
        let crystal_ligand: Vec<[f64; 3]> = ligand_heavy.iter().map(|&k| crystal[k]).collect();
        let lining = protein_heavy
            .iter()
            .copied()
            .filter(|&k| {
                crystal_ligand
                    .iter()
                    .any(|&l| distance(crystal[k], l) < LINING_DISTANCE)
            })
            .collect();
        Ok(Complex {
            at: binding.positions().to_vec(),
            binding: binding.clone(),
            mobile_cutoff: None,
            mobile,
            mobile_residues,
            solvent,
            potential,
            ligand_heavy,
            protein_heavy,
            lining,
            crystal,
        })
    }

    /// The binding it was built from, at the positions it was built with.
    pub fn binding(&self) -> &Binding {
        &self.binding
    }

    /// The mobile cutoff, metres; `None` for a zone given as a mask ([`Complex::with_mobile`]).
    pub fn mobile_cutoff(&self) -> Option<f64> {
        self.mobile_cutoff
    }

    /// Which atoms move, complex order: the ligand and every atom of a mobile residue.
    pub fn mobile(&self) -> &[bool] {
        &self.mobile
    }

    /// How many atoms move.
    pub fn mobile_count(&self) -> usize {
        self.mobile.iter().filter(|&&m| m).count()
    }

    /// The residues with a mobile atom, `A:LEU84`, in the binding's order.
    pub fn mobile_residues(&self) -> &[String] {
        &self.mobile_residues
    }

    /// The solvent.
    pub fn solvent(&self) -> Solvent {
        self.solvent
    }

    /// What the mobile atoms move on: the complex's force field less the terms among frozen atoms
    /// alone, and, with [`Solvent::GeneralizedBorn`], OBC II over every atom. See the module
    /// documentation.
    pub fn potential(&self) -> &ForceField {
        &self.potential
    }

    /// Every atom's position now, metres, complex order.
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.at
    }

    /// The ligand's heavy atoms, complex indices.
    pub fn ligand_heavy_atoms(&self) -> &[usize] {
        &self.ligand_heavy
    }

    /// The protein heavy atoms that line the cavity: within [`LINING_DISTANCE`] of a crystal
    /// ligand heavy atom. Complex indices.
    pub fn lining(&self) -> &[usize] {
        &self.lining
    }

    /// The frozen mask a [`MolecularDynamics`] of this complex must have: every atom that is not
    /// mobile.
    pub fn frozen(&self) -> Vec<bool> {
        self.mobile.iter().map(|m| !m).collect()
    }

    /// Minimises [`Complex::potential`] over the mobile atoms, the buffer frozen, from where they
    /// are, until the largest force on a mobile atom is at most `tolerance` newtons or `max_steps`
    /// steps have been taken. Says where it stopped.
    pub fn minimise(&mut self, max_steps: usize, tolerance: f64) -> Progress {
        let mut m = Minimiser::new(tolerance).with_frozen(self.frozen());
        let mut p = m.step(&self.potential, &[], &mut self.at);
        while p.status == Status::Running && p.steps < max_steps {
            p = m.step(&self.potential, &[], &mut self.at);
        }
        p
    }

    /// Molecular dynamics of this complex in `bath`: the elements' masses, the buffer frozen, and
    /// Maxwell–Boltzmann velocities at `temperature` kelvin for the mobile atoms from `seed`. With
    /// atoms frozen no momentum is removed, so the degrees of freedom are `3 N_mobile`.
    pub fn dynamics(&self, bath: Bath, temperature: f64, seed: u64) -> MolecularDynamics {
        MolecularDynamics::for_elements(self.binding.elements().iter().copied())
            .with_bath(bath)
            .with_frozen(self.frozen())
            .thermalised(&self.at, temperature, seed)
    }

    /// `steps` steps of `dt` seconds of `md` on [`Complex::potential`], a frame into `record` after
    /// each step whose count is a multiple of its interval. The same, to the bit, however the
    /// steps are cut into calls.
    ///
    /// # Panics
    ///
    /// If `md` is not this complex's — one mass per atom and the buffer frozen — or `record` is
    /// not one entry per atom.
    pub fn run(&mut self, md: &mut MolecularDynamics, dt: f64, steps: usize, record: &mut Record) {
        assert_eq!(md.masses().len(), self.at.len(), "one mass per atom");
        assert!(
            md.frozen().iter().zip(&self.mobile).all(|(f, m)| f != m),
            "the dynamics must freeze exactly the buffer"
        );
        assert_eq!(record.sum.len(), self.at.len(), "one record entry per atom");
        // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
        // builds on 1.78.
        #[allow(clippy::manual_is_multiple_of)]
        for _ in 0..steps {
            md.step(&self.potential, &mut self.at, dt);
            if md.steps() % record.every == 0 {
                let frame = self.observe(md);
                record.take(self, frame);
            }
        }
    }

    /// The observables at the current positions, with `md`'s step count, temperature and books.
    /// See [`Frame`].
    pub fn observe(&self, md: &MolecularDynamics) -> Frame {
        let ligand: Vec<[f64; 3]> = self.ligand_heavy.iter().map(|&k| self.at[k]).collect();
        let reference: Vec<[f64; 3]> = self.ligand_heavy.iter().map(|&k| self.crystal[k]).collect();
        let protein: Vec<[f64; 3]> = self.protein_heavy.iter().map(|&k| self.at[k]).collect();
        let lining: Vec<[f64; 3]> = self.lining.iter().map(|&k| self.at[k]).collect();
        let kinds: Vec<Element> = self
            .ligand_heavy
            .iter()
            .map(|&k| self.binding.elements()[k])
            .collect();
        let c = centroid(&ligand);
        let (van_der_waals, electrostatic) = self.binding.cross_terms_at(&self.at);
        let temperature = match md.bath() {
            Bath::Langevin { .. } if md.steps() > 0 => md.half_step_temperature(),
            _ => f64::NAN,
        };
        Frame {
            step: md.steps(),
            rmsd: rmsd(&ligand, &reference),
            site_rmsd: site_rmsd(&ligand, &reference, &kinds),
            turn: in_plane_turn(&ligand, &reference),
            centroid_displacement: distance(c, centroid(&reference)),
            cavity_distance: distance(c, centroid(&lining)),
            contacts: contacts(&ligand, &protein, CONTACT_DISTANCE),
            van_der_waals,
            electrostatic,
            temperature,
            books: md.ledger().get(quantity::ENERGY).unwrap_or(f64::NAN),
        }
    }
}
