//! The Lennard-Jones parameters of substances the potential actually describes.
//!
//! [`LennardJones::reduced`](crate::LennardJones::reduced) sets `σ = ε = m = 1`, which is the
//! right thing for reproducing a paper's state point and the wrong thing for a file: every length
//! in the run is then a number of `σ` written into a format whose positions are **metres**. The
//! shipped fluid scene declared a box `5.0388 x 5.0388 x 5.0388 m` for 108 atoms of argon, which
//! is 1.7 nanometres — the viewer drew it under a scale bar reading `2 M`, and the bar was right
//! about the numbers it was given.
//!
//! A substance is three numbers, and they settle every unit in the run:
//!
//! ```text
//! length σ     energy ε     mass m
//! temperature  ε/k_B        time  σ√(m/ε)        density m/σ³        speed √(ε/m)
//! ```
//!
//! # Why noble gases, and why not all of them
//!
//! The Lennard-Jones potential is a sphere with a well. That is what a closed-shell atom is, and
//! it is **not** what a diatomic, a polar molecule or a metal is: nitrogen has a shape, water has
//! a dipole, and copper's bonding is not pairwise at all. Parameters exist in the literature for
//! all of them and every one is an effective fit that this potential cannot earn.
//!
//! **Two gases were in this catalogue and are not**, for two different reasons, and both are in
//! [`UNSHIPPED`] with the measurement that put them there.
//!
//! | | `Λ*` | triple-point T | triple-point ρ | why it is out |
//! | --- | --- | --- | --- | --- |
//! | argon | 0.186 | 0.8% | 0.4% | — |
//! | krypton | 0.104 | 1.7% | 2.2% | — |
//! | neon | **0.594** | 0.6% | **8.7%** | the model |
//! | xenon | 0.063 | **4.5%** | **10.8%** | the parameters |
//!
//! **Neon's parameters are fine and the potential is wrong for it.** The de Boer parameter
//! `Λ* = h / (σ √(m ε))` is how much of an atom's position is quantum spread compared with the
//! width of its own well, and neon's is three times argon's: at its triple point the thermal de
//! Broglie wavelength is a quarter of `σ`. No choice of `σ` and `ε` fixes that, and helium and
//! hydrogen are further still.
//!
//! **Xenon is the opposite**, and that is how it was caught. It is the *most* classical of the
//! four — `Λ*` of 0.063 — so a 10.8% miss cannot be the potential; it was the pair written down
//! here, which could not be sourced. Solving for the `σ` that reproduces the measured density
//! gives 3.94 Å against the 4.10 that was in the table. **Substituting one that passes is the one
//! thing this check exists to stop**, so xenon waits for a sourced pair rather than a fitted one.
//!
//! # What checks them, and what the check is worth
//!
//! `σ` and `ε` are fitted to a gas's **second virial coefficient and viscosity**. The Lennard-Jones
//! model's own triple point is a simulation result for the potential, known to three figures:
//! `T* = 0.694`, `ρ* = 0.84`. Those two facts together predict each gas's triple point — and a gas's
//! triple point is measured by putting it in a cell and watching it, which has nothing to do with
//! either fit.
//!
//! | | predicted | measured | out by |
//! | --- | --- | --- | --- |
//! | argon, triple-point temperature | 83.1 K | 83.81 K | 0.8% |
//! | argon, triple-point liquid density | 1412 kg/m³ | 1417 kg/m³ | 0.4% |
//!
//! `the_parameters_predict_each_gases_triple_point` holds both shipped gases against both
//! numbers, and `a_gas_this_model_cannot_describe` holds the two that are out against the same
//! arithmetic, so the boundary of the catalogue is a measurement anybody can re-run.
//!
//! **And the same arithmetic says where the model stops.** The critical temperature comes out
//! `T*_c = 1.316` times `ε/k_B`, which for argon is 157.6 K against a measured 150.7 K — **4.6%
//! out**. That is the potential's error and not the parameters', a sphere-with-a-well standing in
//! for an atom near a critical point where fluctuations are everything. The test prints it and
//! asserts only that it stays in the band it has always been in, because tightening that bound
//! would be asserting the model is better than it is.

/// A substance a Lennard-Jones sphere describes, with the three numbers that settle its units.
///
/// The names are the ones a scene writes: `"substance": "argon"`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Substance {
    /// What a scene calls it.
    pub name: &'static str,
    /// The distance at which the potential crosses zero, in metres.
    pub sigma: f64,
    /// The well depth divided by Boltzmann's constant, in kelvin — the form every table quotes.
    pub epsilon_over_k: f64,
    /// One atom's mass, in kilograms.
    pub mass: f64,
    /// The measured triple-point temperature, in kelvin. Not used by the model: this is what the
    /// model is checked *against*.
    pub triple_point_k: f64,
    /// The measured liquid density at the triple point, in kg/m³. Also a check and not an input.
    pub triple_point_density: f64,
    /// The measured critical temperature, in kelvin, which is where this potential is worst.
    pub critical_k: f64,
}

/// Planck's constant, in joule seconds.
const PLANCK: f64 = 6.626_070_15e-34;

/// The Lennard-Jones model's own triple point, in reduced units.
///
/// A simulation result for the potential — not a closed form and not a measurement of anything in
/// the world. It is the bridge: `σ`, `ε` and this give a prediction about a real gas, and the gas
/// has been measured independently.
pub const LJ_TRIPLE_T: f64 = 0.694;
/// The liquid branch of that triple point, as `ρσ³`.
pub const LJ_TRIPLE_DENSITY: f64 = 0.84;
/// The model's critical temperature, as `k_BT/ε`. This is where it is furthest from a real gas.
pub const LJ_CRITICAL_T: f64 = 1.316;

/// Boltzmann's constant, in joules per kelvin.
const BOLTZMANN: f64 = 1.380_649e-23;
/// The atomic mass unit, in kilograms.
const DALTON: f64 = 1.660_539_068_92e-27;

impl Substance {
    /// The well depth in joules, which is what [`LennardJones`](crate::LennardJones) wants.
    pub fn epsilon(&self) -> f64 {
        self.epsilon_over_k * BOLTZMANN
    }

    /// Look one up by the name a scene writes, or `None`.
    ///
    /// `None` rather than a default, for the reason the thermal catalogue gives: a name that
    /// arrived as text is a name that can be wrong, and a substance chosen for somebody because
    /// their spelling missed is worse than a refusal.
    pub fn named(name: &str) -> Option<Substance> {
        CATALOGUE.iter().copied().find(|s| s.name == name)
    }

    /// The temperature this substance's `ε` corresponds to a reduced `t` at.
    pub fn temperature(&self, reduced: f64) -> f64 {
        reduced * self.epsilon_over_k
    }

    /// The mass density, in kg/m³, that a reduced number density `ρσ³` corresponds to.
    pub fn density(&self, reduced: f64) -> f64 {
        reduced * self.mass / (self.sigma * self.sigma * self.sigma)
    }

    /// The de Boer parameter `h / (σ √(m ε))`: how quantum this substance is.
    ///
    /// **The number that decides whether a classical sphere describes it.** It is the thermal de
    /// Broglie wavelength at the well depth, divided by the well's own width, so it says directly
    /// how much of the atom's position is spread out compared with how far apart the atoms are.
    /// Xenon is 0.063 and argon 0.186; neon is 0.594 and misses its triple-point density by 8.66%
    /// for exactly that reason, which is why it is not in the catalogue.
    pub fn de_boer(&self) -> f64 {
        PLANCK / (self.sigma * (self.mass * self.epsilon()).sqrt())
    }
}

/// Every substance this potential is honest about.
///
/// `σ` and `ε/k_B` are the standard Lennard-Jones fits; the triple and critical points beside them
/// are measurements of the gas, carried so the parameters can be checked rather than believed.
pub const CATALOGUE: [Substance; 2] = [
    Substance {
        name: "argon",
        sigma: 3.405e-10,
        epsilon_over_k: 119.8,
        mass: 39.948 * DALTON,
        triple_point_k: 83.81,
        triple_point_density: 1417.0,
        critical_k: 150.69,
    },
    Substance {
        name: "krypton",
        sigma: 3.60e-10,
        epsilon_over_k: 164.0,
        mass: 83.798 * DALTON,
        triple_point_k: 115.78,
        triple_point_density: 2451.0,
        critical_k: 209.48,
    },
];

/// The names, for an error message that lists what a caller could have meant.
pub const NAMES: [&str; 2] = ["argon", "krypton"];

/// The two that did not pass, kept so the catalogue's boundary is a measurement and not a
/// sentence.
///
/// Neon because the potential is wrong for it, xenon because the parameters written down here
/// could not be sourced — see this module's own docs for both, and
/// `a_gas_this_model_cannot_describe` for the arithmetic.
pub const UNSHIPPED: [Substance; 2] = [
    Substance {
        name: "neon",
        sigma: 2.749e-10,
        epsilon_over_k: 35.6,
        mass: 20.1797 * DALTON,
        triple_point_k: 24.56,
        triple_point_density: 1247.0,
        critical_k: 44.49,
    },
    Substance {
        name: "xenon",
        sigma: 4.10e-10,
        epsilon_over_k: 222.0,
        mass: 131.293 * DALTON,
        triple_point_k: 161.40,
        triple_point_density: 2978.0,
        critical_k: 289.73,
    },
];
