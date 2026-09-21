//! pantometry-protein: how a protein's shape moves, and what a bound ligand does to it.
//!
//! A coarse-grained elastic network model. One node per residue at its alpha carbon, one spring
//! between every pair within a cutoff, and the eigenvectors of the resulting Hessian are the
//! collective motions the fold has. [`Protein`] then lets the structure move along them at a
//! temperature, so a viewer can watch a hinge open.
//!
//! Nothing in it is fitted. There are two parameters — a cutoff and a spring constant — and the
//! second cancels out of every shape the crate reports, so what is left is one number taken from
//! the literature.
//!
//! # What it answers
//!
//! - **Which parts of this structure move, and together with what.** [`Modes::shape`] is a
//!   direction in `3n` space; [`Modes::fluctuations`] is how far each residue goes.
//! - **What binding does to that.** [`Network::with_ligand`] adds the ligand as more nodes, and
//!   the answer is a theorem before it is a computation: **no residue can become more mobile on
//!   binding**. Binding stiffens, and the question worth asking is where, by how much, and how
//!   far from the site the effect reaches. The theorem is about the protein's coordinates after
//!   the ligand is integrated out, and not about the spectrum — which gains modes rather than
//!   rising, and whose fifth entry this crate measured going *down*. [`Network::with_ligand`]
//!   has the argument.
//! - **Whether any of it is true.** Two independent measurements, and they say different things
//!   about how much this model has earned.
//!
//! # What it was checked against, and what each check is worth
//!
//! | | measured |
//! | --- | --- |
//! | The lowest mode of open adenylate kinase against the motion it performs on closing | **0.799**, where a random direction gives `0.039` |
//! | The same, over ten modes | **0.966** |
//! | Predicted fluctuations against deposited B-factors | `0.659` lysozyme, `0.571` ubiquitin, `0.395` crambin |
//! | The distance from the centroid against the same B-factors | `0.699`, **`0.804`**, `0.402` |
//!
//! The second pair of rows is the important one and it is not flattering. **A B-factor is one
//! scalar per residue, and one line of geometry predicts them as well as this model does** — on
//! every structure in `structures/`. So the fluctuation comparison is necessary and is not
//! evidence that the collective-motion physics is right.
//!
//! What is evidence is the first pair. Adenylate kinase closes two lids over its substrates, and
//! both ends of that `7.1 Å` motion are deposited. Given only the open structure, the **single
//! softest mode** of this model overlaps what the enzyme actually does by `0.799` in a
//! 642-dimensional space where chance is `0.039`. That is a claim about a *direction* — which
//! residues move together and which way — and no per-residue scalar can produce it by accident.
//!
//! # What it does not answer, and would be wrong to be asked
//!
//! - **Binding affinity.** There is no solvent, no entropy, no electrostatics and no chemistry
//!   here. Nothing in this crate says how *tightly* anything binds, and a number that looked like
//!   it did would be an artefact of the spring constant.
//! - **Where the ligand goes.** The pose is an input. This is not docking; there is no search and
//!   no scoring function.
//! - **Large conformational change.** This is a harmonic model about one structure, so it
//!   describes small motions around the coordinates it was given. The low modes often *point*
//!   along a large transition, which is itself checkable against a second structure of the same
//!   protein in the other state — but the model cannot travel along one.
//! - **Anything at atomic detail.** No side chains, no hydrogens, no bond breaking, no folding.
//!
//! All-atom molecular dynamics answers some of those and is not what this is. It needs a force
//! field's parameter tables, particle-mesh electrostatics, constraint algorithms and a GPU; this
//! workspace resolves twelve external crates, compiles every one of them to `wasm32` and Rust
//! 1.78, and keeps no GPU in the library. What is here instead is the thing that can be checked.
//!
//! # The one thing to remember when reading a mode
//!
//! **A mode's sign means nothing.** `v` and `−v` are the same motion. [`spectrum`] fixes a sign
//! so that results are bit-identical everywhere, and that rule is a coin toss whenever the
//! structure is symmetric — measured, five of twenty modes of a twenty-bead chain. Anything read
//! off a mode has to be invariant under `v → −v`: [`Modes::overlap`] is an absolute value for
//! this reason and not as a convenience.

#![deny(missing_docs)]

pub mod modes;
pub mod network;
pub mod spectrum;
pub mod structure;

use pantometry_core::conserved::quantity;
use pantometry_core::{Bodies, Domain, Exchange, Kind, Ledger, Reading, Rng, Violation};
use pantometry_units::{Length, LengthVec, Mass, Qty, Temperature, Time, BOLTZMANN};

pub use modes::{correlation, Modes};
pub use network::Network;
pub use structure::{Atom, ParseError, Residue, Structure};

/// A structure moving along its own normal modes at a temperature.
///
/// # How it moves
///
/// Every non-rigid mode is an independent harmonic oscillator, so the motion is written down
/// rather than integrated:
///
/// ```text
/// r_i(t) = r⁰_i + Σ_k A_k u_k,i cos(ω_k t + φ_k),   A_k = √(2k_BT/λ_k),   ω_k = √(λ_k/m)
/// ```
///
/// There is no integrator and therefore no integrator error, no stability limit, and no drift:
/// [`Domain::max_stable_dt`] is infinite and it means it. The amplitude is set by equipartition,
/// so each mode carries `k_BT` — half kinetic, half potential — and the total energy is exactly
/// `(3n − 6) k_BT` at every instant. That is what [`Domain::ledger`] reports and what makes
/// [`Domain::books_balance`] true rather than approximately true.
///
/// The phases come from [`Rng::for_index`] on a seed, so a run is reproducible bit for bit and
/// two runs with different seeds are two members of the same thermal ensemble.
#[derive(Clone, Debug)]
pub struct Protein {
    name: String,
    rest: Vec<[f64; 3]>,
    at: Vec<[f64; 3]>,
    modes: Modes,
    amplitude: Vec<f64>,
    omega: Vec<f64>,
    phase: Vec<f64>,
    fluctuation: Vec<f64>,
    temperature: f64,
    energy: f64,
}

impl Protein {
    /// A protein shaking at `temperature`, with `mass` on every node and phases drawn from
    /// `seed`.
    ///
    /// `mass` is the effective mass of a residue — around 110 daltons, `1.8e-25` kg — and it sets
    /// the *rate* of the motion and nothing else: every amplitude, every fluctuation and every
    /// energy here is independent of it. A heavier node makes the same picture run slower.
    ///
    /// # Panics
    ///
    /// If the temperature or the mass is not positive, or if the network is in more than one
    /// piece — a disconnected network's lowest mode is one piece drifting away from another, and
    /// every fluctuation it predicts is that drift rather than the protein.
    pub fn new(
        name: impl Into<String>,
        network: &Network,
        temperature: Temperature,
        mass: Mass,
        seed: u64,
    ) -> Protein {
        let pieces = network.components();
        assert_eq!(
            pieces,
            1,
            "the network falls into {pieces} pieces at a cutoff of {:e} m; each of them has six \
             rigid-body modes of its own and the softest motion is one drifting off",
            network.cutoff().to_si()
        );
        let (t, m) = (temperature.to_si(), mass.to_si());
        assert!(t > 0.0, "a temperature of {t} K");
        assert!(m > 0.0, "a mass of {m} kg");

        let modes = Modes::of(network);
        let kt = BOLTZMANN.to_si() * t;
        let count = modes.len();
        let mut amplitude = Vec::with_capacity(count);
        let mut omega = Vec::with_capacity(count);
        let mut phase = Vec::with_capacity(count);
        for k in 0..count {
            let lambda = modes.eigenvalue(k).to_si();
            amplitude.push((2.0 * kt / lambda).sqrt());
            omega.push((lambda / m).sqrt());
            phase.push(Rng::for_index(seed, k as u64).unit() * std::f64::consts::TAU);
        }
        let rest: Vec<[f64; 3]> = (0..network.nodes()).map(|i| network.at(i)).collect();
        let mut protein = Protein {
            name: name.into(),
            at: rest.clone(),
            rest,
            fluctuation: modes.fluctuations(temperature),
            modes,
            amplitude,
            omega,
            phase,
            temperature: t,
            // Equipartition: `k_BT` in each mode, and it does not move.
            energy: count as f64 * kt,
        };
        protein.place(0.0);
        protein
    }

    /// The modes underneath.
    pub fn modes(&self) -> &Modes {
        &self.modes
    }

    /// `√⟨Δr²⟩` for each node, in metres: the equilibrium spread, not the current displacement.
    pub fn spread(&self) -> Vec<f64> {
        self.fluctuation.iter().map(|x| x.sqrt()).collect()
    }

    /// How far node `i` is from where it started, right now, in metres.
    ///
    /// # Panics
    ///
    /// If `i` is at or past the node count.
    pub fn displacement(&self, i: usize) -> Length {
        Qty::from_si(
            (0..3)
                .map(|k| (self.at[i][k] - self.rest[i][k]).powi(2))
                .sum::<f64>()
                .sqrt(),
        )
    }

    /// Put every node where time `t` says it is.
    fn place(&mut self, t: f64) {
        self.at.copy_from_slice(&self.rest);
        for k in 0..self.modes.len() {
            let swing = self.amplitude[k] * (self.omega[k] * t + self.phase[k]).cos();
            let shape = self.modes.shape(k);
            for (i, node) in self.at.iter_mut().enumerate() {
                node[0] += swing * shape[3 * i];
                node[1] += swing * shape[3 * i + 1];
                node[2] += swing * shape[3 * i + 2];
            }
        }
    }
}

impl Bodies for Protein {
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

    /// The equilibrium spread, in metres, so a run is coloured by which residues are mobile
    /// rather than by where each one happens to be at this frame.
    fn value(&self, i: usize) -> f64 {
        self.fluctuation[i].sqrt()
    }

    fn value_unit(&self) -> &'static str {
        "m"
    }
}

impl Domain for Protein {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// Infinite, and not as a placeholder: the motion is evaluated in closed form at whatever
    /// time it is asked for, so no step size is unstable and none is more accurate than another.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Qty::from_si(f64::INFINITY)
    }

    fn step(&mut self, t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        self.place(t.to_si() + dt.to_si());
        Ok(())
    }

    /// `(3n − 6) k_BT`, which does not change. Half of it is the potential in the springs and
    /// half is kinetic, and the split moves between them while the total does not.
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, self.energy)
    }

    fn books_balance(&self) -> bool {
        true
    }

    fn readings(&self) -> Vec<Reading> {
        let softest = if self.modes.is_empty() {
            0.0
        } else {
            self.omega[0] / std::f64::consts::TAU
        };
        let nodes = self.at.len() as f64;
        let spread = (self.fluctuation.iter().sum::<f64>() / nodes).sqrt();
        let now = (0..self.at.len())
            .map(|i| self.displacement(i).to_si().powi(2))
            .sum::<f64>()
            / nodes;
        vec![
            Reading::new(&self.name, "softest mode", softest, "Hz"),
            Reading::new(&self.name, "spread", spread, "m"),
            Reading::new(&self.name, "displacement", now.sqrt(), "m"),
            Reading::new(&self.name, "energy", self.energy, "J"),
            Reading::new(&self.name, "temperature", self.temperature - 273.15, "C"),
            Reading::new(
                &self.name,
                "rigid modes",
                self.modes.rigid_body_modes() as f64,
                "",
            ),
            Reading::new(&self.name, "mode separation", self.modes.separation(), ""),
        ]
    }

    /// Neither is an answer about the protein. `rigid modes` is six when the null space came out
    /// the size the geometry says, and `mode separation` is the evidence that calling those six
    /// zero was a measurement rather than a threshold — see [`Modes::separation`].
    fn diagnostics(&self) -> &'static [&'static str] {
        &["rigid modes", "mode separation"]
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
