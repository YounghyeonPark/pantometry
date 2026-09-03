//! pantometry-pharmacokinetic: where a drug goes in a body, as a domain on the
//! `pantometry-core` kernel.
//!
//! A [`CompartmentModel`] is *n* well-stirred compartments joined by intercompartmental
//! clearances, with elimination out of the ones that eliminate. A compartment has a volume and
//! holds a mass of drug; its concentration is the one over the other. Between two compartments an
//! intercompartmental clearance `Q` carries `Q·(C_i − C_j)` kilograms a second, and on a
//! compartment an elimination clearance `CL` removes `CL·C` from the body altogether.
//!
//! ```text
//! V_i dC_i/dt = Σ_j Q_ij (C_j − C_i) − CL_i C_i + R_i(t)
//! ```
//!
//! That is the same algebra a thermal network runs — `Q·ΔC` where a joint carries `UA·ΔT` — and
//! this crate copies its structure deliberately. Two things are new.
//!
//! **Elimination.** A conductance moves heat between two nodes; a clearance moves drug out of the
//! system. The mass is gone from the body and must not be gone from the books, which is what
//! [`CompartmentModel::ledger`] is about below.
//!
//! **Dosing.** Heat arrives on a bus from another domain. A dose arrives from a clinician, so it
//! is declared rather than exchanged: [`CompartmentModel::bolus`] puts an amount somewhere now,
//! and [`CompartmentModel::infuse`] schedules one over a window.
//!
//! # What it is checked against
//!
//! Compartmental pharmacokinetics is one of the few pieces of physiology with exact solutions,
//! which is the reason this is the half of "medical simulation" that arrived first. Every claim
//! here is checked against one of them and none against a second implementation:
//!
//! - **One compartment, IV bolus.** `A(t) = A₀ e^{−kt}` with `k = CL/V`. Exact.
//! - **One compartment, constant infusion.** `A(t) = (R/k)(1 − e^{−kt})`, and a steady state of
//!   `R·V/CL` that explicit Euler reproduces *exactly* rather than approximately — the fixed
//!   point of the discrete map is the fixed point of the differential equation, at any step.
//! - **Two compartments, IV bolus.** The bi-exponential, with `α` and `β` the roots of
//!   `λ² − (k10 + k12 + k21)λ + k10 k21 = 0`. This is the one that catches a transposed index in
//!   a link, which the conservation audit structurally cannot — see below.
//! - **Mass balance.** Dose = still in the compartments + cleared + still in the syringe, to
//!   machine precision.
//! - **First order.** Halving the step halves the error. Measured as a rate, so the tolerances
//!   above are consequences rather than choices.
//!
//! # What the audit cannot see here, and what covers it instead
//!
//! A link contributes `+q` to one compartment and `−q` to another **in the same sum**. They
//! cancel identically, so the ledger is blind to links by construction: a sign error, a
//! transposed index, or a link dropped altogether passes the conservation audit at machine
//! precision. `pantometry-thermal`'s network module found this first and the lesson transfers
//! whole — every link check in this crate is per compartment or against a closed form, never on
//! the total.
//!
//! It decided the API the same way, too. A compartment is addressed by a [`Compartment`] handle
//! rather than by name, because a link naming a compartment that does not exist would be exactly
//! that invisible case. A handle can only come from [`CompartmentModel::compartment`] or
//! [`CompartmentModel::eliminating`], so a dangling reference is not representable, and
//! [`CompartmentModel::compartment_named`] is the bridge for a caller holding names off a file.
//!
//! # The ledger holds what was cleared, and that is a decision
//!
//! A ledger reports what a domain is **holding**. Drug removed by a clearance has left the body,
//! so the compartments no longer hold it — and if the books said only what the compartments hold,
//! a perfectly correct model would report a leak every step.
//!
//! There are two ways out and this crate takes the second. It could publish the cleared mass on
//! the [`Exchange`], but there is nothing on the other side: excretion leaves the *model*, not
//! just the body, and a channel nobody consumes is refused by the bus's own audit. So the cleared
//! mass stays on the books as a running total reported beside the compartments, and the sum is
//! the dose. That sum is exactly what a mass-balance test checks, which is the second reason.
//!
//! The two ways to get this wrong are both on record elsewhere in this workspace, and neither can
//! happen here for the same reason: **nothing is published at all.** A thermal domain once
//! subtracted what it absorbed *and* reported what it stored, so the entry cancelled itself; a
//! mechanics domain kept counting joules it had already handed to the bus, and energy grew 63%.
//! This domain's drug never leaves its books — only the body's — so there is no traffic to
//! double-count against.
//!
//! Scheduled infusions are on the books too, as [`CompartmentModel::pending_dose`]. A syringe
//! with 400 mg still in it is 400 mg the model is holding, and counting it is what makes mass
//! balance checkable *during* an infusion rather than only after one.
//!
//! # What is deliberately not in it
//!
//! **No absorption compartment and no oral dosing.** First-order absorption is `dA/dt = −kₐA`, a
//! rate constant acting on an *amount*, and everything here is a clearance acting on a
//! *concentration*. They are different mechanisms, not one mechanism with a different number, and
//! a depot faked as a compartment with a large volume would be a fiction the closed forms would
//! not check.
//!
//! **No Michaelis–Menten elimination.** `V_max C/(K_m + C)` is the interesting case for ethanol
//! and phenytoin, and it is exactly where the exact solutions stop: there is no closed form for a
//! saturable two-compartment model, and this workspace's README says what it thinks of a domain
//! with nothing to check against.
//!
//! **No pharmacodynamics.** An `E_max` curve, an effect compartment, a receptor occupancy — those
//! answer "what does the concentration *do*", which is a response model rather than a mass
//! balance, and it has no conserved quantity for the audit to hold.
//!
//! **No population variability.** Inter-individual spread is the reason a real PK study exists,
//! and it is a sampling problem: draw a thousand subjects from a distribution of `CL` and `V` and
//! look at the spread. That is a different execution axis from the one this kernel has — it
//! advances `dt` and audits conservation, where a Monte Carlo run samples — and putting it here
//! would be a loop over a domain rather than a domain.
//!
//! **No protein binding, no metabolites, no AUC as a typed quantity.** The first two are more
//! compartments a caller can build out of what is here. The third is `Σ eliminated / CL` and
//! would need a `Qty` this workspace has no name for; it is one division away in a caller.
//!
//! # The twelfth domain, and what it had to add to the kernel
//!
//! Nothing, which is the claim the crate split exists to test and the twelfth time it has held.
//! It did need one thing from `pantometry-units`, which is a different crate and a different
//! claim: a clearance is a **volume per time** and there was no name for that dimension, so
//! [`pantometry_units::VolumetricFlow`] is new and [`Clearance`] is this domain's name for it.

// Every public item carries a doc comment. Denied rather than warned: a public physics API
// whose `Length::mm` shows a blank summary in rustdoc is documented in the sense that a
// paragraph exists somewhere, and not in the sense a reader needs.
#![deny(missing_docs)]

use std::any::Any;

use pantometry_core::conserved::quantity;
use pantometry_core::{Domain, Exchange, Kind, Ledger, Reading, Violation};
use pantometry_units::{Concentration, Frequency, Mass, MassFlow, Time, Volume};

/// A clearance: m³·s⁻¹, a volume of plasma emptied of drug per unit time.
///
/// The same dimension as [`pantometry_units::VolumetricFlow`], which is what this is an alias
/// for, and deliberately a distinct name — the way [`Concentration`] is a distinct name for
/// [`Density`](pantometry_units::Density). A clearance is *not* a rate constant, and the
/// difference is the most common confusion in the subject: `k = CL/V` has units of one over time
/// and depends on the volume it is quoted against, while `CL` does not. The type system holds the
/// two apart here, and `CL·C` coming out as a [`MassFlow`] is that statement made mechanically.
///
/// Constructed in the units a clinician uses — `Clearance::l_per_h(10.0)`,
/// `Clearance::ml_per_min(120.0)` — because m³/s is nobody's unit for this and the gap between
/// L/h and m³/s is a factor of 3.6 million.
pub type Clearance = pantometry_units::VolumetricFlow;

/// A compartment in one particular model.
///
/// Carries the model's identity as well as the index, so a handle from one model used on another
/// is refused rather than silently addressing whatever sits at that position. The identity is a
/// hash of the model's name, which is deterministic: no counter, no clock, no global state.
///
/// **The limit of that, stated because it is easy to assume otherwise.** Two models given the same
/// *name* hash to the same identity, and their handles are then interchangeable — the check sees
/// two models called `"subject"` as one. Nothing above prevents it:
/// [`Simulation`](pantometry_core::Simulation) accepts two domains with one name and
/// [`Simulation::domain`](pantometry_core::Simulation::domain) returns the first, measured rather
/// than assumed. So a duplicate name is already a mistake at the layer above, and this is one more
/// thing it costs rather than a case this type can catch. Name your models distinctly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Compartment {
    index: u32,
    model: u64,
}

/// One well-stirred volume holding some drug, and possibly clearing it.
struct CompartmentState {
    label: String,
    /// m³. A compartment is an apparent volume, not a place: the "central compartment" is plasma
    /// plus everything that equilibrates with it faster than the sampling can see.
    volume: f64,
    /// kg of drug in it.
    amount: f64,
    /// m³/s out of the *system*. Zero for a compartment that only distributes.
    clearance: f64,
}

/// An intercompartmental clearance, in m³/s.
struct Link {
    a: usize,
    b: usize,
    q: f64,
}

/// A constant-rate infusion over a window.
struct Infusion {
    into: usize,
    /// kg/s.
    rate: f64,
    /// Seconds on the simulation clock.
    start: f64,
    end: f64,
    /// What the whole window is worth, in kg. Held rather than recomputed so the last slice can
    /// deliver the exact remainder — see [`CompartmentModel::step`].
    total: f64,
    delivered: f64,
}

/// *n* compartments, intercompartmental clearances between them, elimination out of some of them.
///
/// # Why one domain and not one domain per compartment
///
/// A link carries `Q·(C_i − C_j)`: it needs **both** concentrations. Domains in this workspace
/// never read each other — they meet on an [`Exchange`], which carries *amounts* and not state —
/// so neither side could publish a concentration and neither could compute the flux alone. Any
/// `distributing_to(peer)` would break the property the rest of the design rests on.
///
/// So the model is one domain holding many compartments, which is also what a compartmental model
/// physically is: a single coupled system of ODEs, not independent volumes posting parcels to
/// each other. One `ledger`, one stability limit, and the conservation audit unchanged. The same
/// argument, in the same words, as `pantometry-thermal`'s network — because it is the same shape.
///
/// ```
/// use pantometry_pharmacokinetic::{Clearance, CompartmentModel};
/// use pantometry_core::{Domain, Exchange};
/// use pantometry_units::{Mass, Time, Volume};
///
/// // One compartment: 10 L of distribution volume, 5 L/h of clearance. k = 0.5 per hour.
/// let mut body = CompartmentModel::new("subject");
/// let plasma = body.eliminating(
///     "plasma",
///     Volume::litres(10.0),
///     Clearance::l_per_h(5.0),
/// );
/// body.bolus(plasma, Mass::g(0.1)).unwrap();       // 100 mg IV
///
/// let dt = Time::s(30.0);
/// let mut bus = Exchange::new();
/// let mut t = Time::ZERO;
/// while t.to_si() < 4.0 * 3600.0 {
///     body.step(t, dt, &mut bus).unwrap();
///     t = t + dt;
/// }
///
/// // Four hours is two elimination time constants: e^-2 of the dose is left.
/// let left = body.amount(plasma).to_si();
/// assert!((left / 1e-4 - (-2.0f64).exp()).abs() < 1e-2, "{left}");
/// // And what is not in the body is on the books as cleared.
/// let books = left + body.eliminated_mass().to_si();
/// assert!((books / 1e-4 - 1.0).abs() < 1e-12, "{books}");
/// ```
pub struct CompartmentModel {
    name: String,
    id: u64,
    compartments: Vec<CompartmentState>,
    links: Vec<Link>,
    infusions: Vec<Infusion>,
    /// kg removed from the body by clearance over the run.
    eliminated: f64,
    /// kg put into the body over the run by **bolus**. The infused half is not accumulated
    /// here: it is `Σ inf.delivered`, which the last slice of each window snaps to the scheduled
    /// total exactly. A second running sum of the same quantity would drift from the first by
    /// `N·ε` over `N` slices — measured at 4.2e-13 of the dose over 2 928 of them — and then two
    /// numbers that must agree would not.
    bolused: f64,
    saved: Option<Saved>,
}

/// Everything `ledger` reads, which is everything `checkpoint` has to save.
///
/// All of it, because saving the state and not the running totals is how a rewound sweep comes to
/// report drug it never cleared. `LumpedMass` did exactly that with its `lost` until an iterative
/// coupling finally reached the branch, and an infusion's `delivered` is the same hazard wearing
/// different clothes: restore without it and the syringe refills.
type Saved = (Vec<f64>, Vec<f64>, f64, f64);

/// A deterministic identity for a model, from its name.
fn identity(name: &str) -> u64 {
    // FNV-1a. Not for security; for telling two models apart without a counter or a clock, which
    // the determinism rule forbids.
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for b in name.as_bytes() {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

impl CompartmentModel {
    /// An empty model. Add compartments, then link them, then dose them.
    pub fn new(name: impl Into<String>) -> CompartmentModel {
        let name = name.into();
        let id = identity(&name);
        CompartmentModel {
            name,
            id,
            compartments: Vec::new(),
            links: Vec::new(),
            infusions: Vec::new(),
            eliminated: 0.0,
            bolused: 0.0,
            saved: None,
        }
    }

    /// A compartment that only distributes: drug reaches it and leaves it through links, and
    /// nothing is eliminated from it.
    ///
    /// The peripheral compartments of every classical model are this. Giving one a clearance of
    /// zero would be the same thing spelled worse, and absence of an elimination path is the
    /// absence of a thing — the way `pantometry-thermal` spells an interior node.
    pub fn compartment(&mut self, label: impl Into<String>, volume: Volume) -> Compartment {
        self.push(label, volume, 0.0)
    }

    /// A compartment that also eliminates: `CL·C` kilograms a second leave the body from it.
    ///
    /// Almost always the central one. Renal and hepatic clearance both act on plasma, and a
    /// peripheral compartment with its own clearance is a modelling choice worth making
    /// explicitly rather than by default.
    pub fn eliminating(
        &mut self,
        label: impl Into<String>,
        volume: Volume,
        clearance: Clearance,
    ) -> Compartment {
        self.push(label, volume, clearance.to_si().max(0.0))
    }

    fn push(&mut self, label: impl Into<String>, volume: Volume, clearance: f64) -> Compartment {
        let index = self.compartments.len() as u32;
        self.compartments.push(CompartmentState {
            label: label.into(),
            volume: volume.to_si(),
            amount: 0.0,
            clearance,
        });
        Compartment {
            index,
            model: self.id,
        }
    }

    /// Join two compartments by an intercompartmental clearance `Q`, in volume per time.
    ///
    /// Refuses a self-link, a negative or non-finite `Q`, and a handle from a different model.
    /// Two links between the same pair **accumulate**, because parallel clearances add — the same
    /// convention [`Exchange::publish`](pantometry_core::Exchange::publish) follows for repeated
    /// offers on one channel.
    ///
    /// The flux is symmetric in the two concentrations and antisymmetric in sign, so the order of
    /// the arguments does not change the physics. It does change nothing else either, which is
    /// worth saying because the micro rate constants it implies are *not* symmetric:
    /// `k₁₂ = Q/V₁` and `k₂₁ = Q/V₂` differ whenever the volumes do.
    pub fn link(&mut self, a: Compartment, b: Compartment, q: Clearance) -> Result<(), Violation> {
        let (i, j) = (self.resolve(a)?, self.resolve(b)?);
        if i == j {
            return Err(Violation::at(
                format!("{}/{}", self.name, self.compartments[i].label),
                "a compartment cannot exchange with itself",
                0.0,
            ));
        }
        let w = q.to_si();
        // `!is_finite` first, so NaN is rejected by the branch that reads as rejecting it rather
        // than by a negated comparison that happens to be false.
        if !w.is_finite() || w < 0.0 {
            return Err(Violation::at(
                format!(
                    "{}/{}-{}",
                    self.name, self.compartments[i].label, self.compartments[j].label
                ),
                "an intercompartmental clearance must be finite and not negative",
                w,
            ));
        }
        if let Some(existing) = self
            .links
            .iter_mut()
            .find(|l| (l.a == i && l.b == j) || (l.a == j && l.b == i))
        {
            existing.q += w;
        } else {
            self.links.push(Link { a: i, b: j, q: w });
        }
        Ok(())
    }

    /// An intravenous bolus: an amount placed in a compartment **now**.
    ///
    /// A dose is not something another domain hands over, so it does not come off the bus. It
    /// comes from a caller, between steps, and that is the one place it is safe: the conservation
    /// audit compares the books either side of an [`advance`](pantometry_core::Simulation::advance),
    /// so an addition made outside one is seen as the new starting state rather than as mass from
    /// nowhere. Inside a step it would be exactly that. Reach a running model through
    /// [`Domain::as_any_mut`], which is what that method is for.
    ///
    /// Refuses a negative or non-finite dose. A zero one is allowed and does nothing, because a
    /// scheduler stepping through a regimen should not have to special-case the empty dose.
    pub fn bolus(&mut self, into: Compartment, dose: Mass) -> Result<(), Violation> {
        let i = self.resolve(into)?;
        let m = dose.to_si();
        if !m.is_finite() || m < 0.0 {
            return Err(Violation::at(
                format!("{}/{}", self.name, self.compartments[i].label),
                "a dose must be finite and not negative",
                m,
            ));
        }
        self.compartments[i].amount += m;
        self.bolused += m;
        Ok(())
    }

    /// An infusion of `dose` delivered at a constant rate over `[start, end)`, on the simulation
    /// clock.
    ///
    /// The clinical ordering — "500 mg over thirty minutes" — and the rate is `dose/(end − start)`.
    /// [`CompartmentModel::infuse_at`] is the same thing said as a rate.
    ///
    /// The whole `dose` is on the books from the moment it is scheduled, as
    /// [`pending_dose`](CompartmentModel::pending_dose), and moves to the compartment as it is
    /// delivered. That is what makes mass balance checkable *during* an infusion: without it the
    /// books would grow every step and a correct model would look like it was inventing drug.
    ///
    /// Refuses a window that does not run forwards, a non-finite bound, and a negative dose. An
    /// infinite window is refused too, and not as an oversight: an infusion with no end has an
    /// unbounded reserve, and a ledger entry of `inf` is a ledger that can no longer be audited.
    pub fn infuse(
        &mut self,
        into: Compartment,
        dose: Mass,
        start: Time,
        end: Time,
    ) -> Result<(), Violation> {
        let (t0, t1) = (start.to_si(), end.to_si());
        let span = t1 - t0;
        if !t0.is_finite() || !t1.is_finite() || span <= 0.0 {
            return Err(Violation::at(
                self.name.clone(),
                "an infusion window must be finite and run forwards",
                span,
            ));
        }
        let m = dose.to_si();
        if !m.is_finite() || m < 0.0 {
            return Err(Violation::at(
                self.name.clone(),
                "an infused dose must be finite and not negative",
                m,
            ));
        }
        let i = self.resolve(into)?;
        self.infusions.push(Infusion {
            into: i,
            rate: m / span,
            start: t0,
            end: t1,
            total: m,
            delivered: 0.0,
        });
        Ok(())
    }

    /// The same, given a rate rather than a total. `dose = rate·(end − start)`.
    pub fn infuse_at(
        &mut self,
        into: Compartment,
        rate: MassFlow,
        start: Time,
        end: Time,
    ) -> Result<(), Violation> {
        let span = end.to_si() - start.to_si();
        if !span.is_finite() || span <= 0.0 {
            return Err(Violation::at(
                self.name.clone(),
                "an infusion window must be finite and run forwards",
                span,
            ));
        }
        self.infuse(into, Mass::from_si(rate.to_si() * span), start, end)
    }

    fn resolve(&self, c: Compartment) -> Result<usize, Violation> {
        if c.model != self.id {
            return Err(Violation::at(
                self.name.clone(),
                "a compartment handle from a different model",
                c.index as f64,
            ));
        }
        let i = c.index as usize;
        if i >= self.compartments.len() {
            return Err(Violation::at(
                self.name.clone(),
                "a compartment handle from a later state of this model",
                c.index as f64,
            ));
        }
        Ok(i)
    }

    /// How much drug is in a compartment. Total, because a handle cannot name one that is absent.
    pub fn amount(&self, c: Compartment) -> Mass {
        Mass::from_si(self.compartments[c.index as usize].amount)
    }

    /// A compartment's concentration, `A/V` — the quantity a blood sample measures and the one
    /// every closed form here is written in.
    pub fn concentration(&self, c: Compartment) -> Concentration {
        let s = &self.compartments[c.index as usize];
        Concentration::from_si(if s.volume > 0.0 {
            s.amount / s.volume
        } else {
            0.0
        })
    }

    /// A compartment's apparent volume.
    pub fn volume(&self, c: Compartment) -> Volume {
        Volume::from_si(self.compartments[c.index as usize].volume)
    }

    /// A compartment's own elimination clearance. Zero for one that only distributes.
    pub fn clearance(&self, c: Compartment) -> Clearance {
        Clearance::from_si(self.compartments[c.index as usize].clearance)
    }

    /// The label a compartment was given, for a violation or a legend.
    pub fn label(&self, c: Compartment) -> &str {
        &self.compartments[c.index as usize].label
    }

    /// A compartment by label, for a caller that built the model from a file and has names rather
    /// than handles. The seam between a name-shaped format and a handle-shaped API.
    pub fn compartment_named(&self, label: &str) -> Option<Compartment> {
        self.compartments
            .iter()
            .position(|c| c.label == label)
            .map(|i| Compartment {
                index: i as u32,
                model: self.id,
            })
    }

    /// Every compartment, in the order they were added, with its label.
    ///
    /// A [`Compartment`] can only come from a constructor on this model, which is what makes a
    /// dangling link unrepresentable — and would also leave a caller holding a model it did not
    /// build unable to walk it at all.
    pub fn handles(&self) -> impl Iterator<Item = (Compartment, &str)> + '_ {
        let model = self.id;
        self.compartments.iter().enumerate().map(move |(i, c)| {
            (
                Compartment {
                    index: i as u32,
                    model,
                },
                c.label.as_str(),
            )
        })
    }

    /// How many compartments.
    pub fn compartments(&self) -> usize {
        self.compartments.len()
    }

    /// Drug in the body right now, summed over the compartments.
    pub fn body_burden(&self) -> Mass {
        Mass::from_si(self.compartments.iter().map(|c| c.amount).sum())
    }

    /// Drug removed from the body by clearance over the run.
    ///
    /// Still on this domain's books — see the module docs. `body_burden + eliminated_mass +
    /// pending_dose` is the dose, exactly.
    pub fn eliminated_mass(&self) -> Mass {
        Mass::from_si(self.eliminated)
    }

    /// Drug given over the run, by bolus and by delivered infusion.
    ///
    /// Summed from the infusions rather than tallied as it goes, so it agrees exactly with
    /// `scheduled − pending` instead of to within an accumulation of roundings.
    pub fn administered_mass(&self) -> Mass {
        Mass::from_si(self.bolused + self.infusions.iter().map(|i| i.delivered).sum::<f64>())
    }

    /// Scheduled infusion that has not been delivered yet — what is still in the syringe.
    pub fn pending_dose(&self) -> Mass {
        Mass::from_si(
            self.infusions
                .iter()
                .map(|i| (i.total - i.delivered).max(0.0))
                .sum(),
        )
    }

    /// Steady-state volume of distribution, `Σ V_i`.
    ///
    /// The volume that relates the body burden to the plasma concentration once everything has
    /// equilibrated, which is the only time it relates them at all. Early after a bolus the
    /// central volume is the relevant one and this is an overestimate by however much has not
    /// distributed yet.
    pub fn volume_of_distribution(&self) -> Volume {
        Volume::from_si(self.compartments.iter().map(|c| c.volume).sum())
    }

    /// Total elimination clearance, `Σ CL_i`.
    ///
    /// At steady state under a constant infusion `Σ CL_i C_i = R`, whatever the distribution
    /// between compartments — so with elimination from one compartment only, its steady-state
    /// concentration is `R/CL` and does not depend on the intercompartmental clearances at all.
    pub fn total_clearance(&self) -> Clearance {
        Clearance::from_si(self.compartments.iter().map(|c| c.clearance).sum())
    }

    /// The elimination rate constant of a compartment, `k = CL/V`.
    ///
    /// For a one-compartment model this is *the* rate constant and `ln2/k` is the half-life. For
    /// anything larger it is a micro constant: the terminal half-life is `ln2/β` where `β` is the
    /// smaller root of the characteristic polynomial, and the two are not close. This crate
    /// deliberately offers no `half_life`, because which half-life a reader means depends on how
    /// many compartments there are and a method could only pick one.
    pub fn elimination_rate(&self, c: Compartment) -> Frequency {
        let s = &self.compartments[c.index as usize];
        Frequency::from_si(if s.volume > 0.0 {
            s.clearance / s.volume
        } else {
            0.0
        })
    }

    /// The micro rate constant `k_ij = Q/V_i` for the link out of `from` towards `to`.
    ///
    /// Zero if they are not linked, which is a real answer rather than a missing one. Not
    /// symmetric: the same `Q` gives a different constant in each direction whenever the volumes
    /// differ, and that asymmetry is what makes the bi-exponential bi-exponential.
    pub fn transfer_rate(&self, from: Compartment, to: Compartment) -> Frequency {
        let (Ok(i), Ok(j)) = (self.resolve(from), self.resolve(to)) else {
            return Frequency::from_si(0.0);
        };
        let q = self
            .links
            .iter()
            .find(|l| (l.a == i && l.b == j) || (l.a == j && l.b == i))
            .map(|l| l.q)
            .unwrap_or(0.0);
        let v = self.compartments[i].volume;
        Frequency::from_si(if v > 0.0 { q / v } else { 0.0 })
    }

    /// Mass crossing the link between two compartments right now, positive from `a` to `b`.
    ///
    /// Zero if they are not linked. The per-link number the conservation audit cannot see, which
    /// is why it is public: a test asserting it is the only kind that would notice a transposed
    /// index.
    pub fn transfer(&self, a: Compartment, b: Compartment) -> MassFlow {
        let (Ok(i), Ok(j)) = (self.resolve(a), self.resolve(b)) else {
            return MassFlow::from_si(0.0);
        };
        let q = self
            .links
            .iter()
            .find(|l| (l.a == i && l.b == j) || (l.a == j && l.b == i))
            .map(|l| l.q)
            .unwrap_or(0.0);
        MassFlow::from_si(q * (self.conc(i) - self.conc(j)))
    }

    fn conc(&self, i: usize) -> f64 {
        let s = &self.compartments[i];
        if s.volume > 0.0 {
            s.amount / s.volume
        } else {
            0.0
        }
    }

    /// Everything carrying drug *out* of compartment `i`, as a rate constant: `(ΣQ + CL)/V`.
    ///
    /// The magnitude of the diagonal of the compartment matrix written in concentrations, and the
    /// only thing [`Domain::max_stable_dt`] needs. Infinite for a compartment with no volume,
    /// which `step` refuses before this is ever divided by.
    fn turnover_rate(&self, i: usize) -> f64 {
        let s = &self.compartments[i];
        // Written as a positive test rather than a negated one, so a NaN volume takes the
        // "no volume" branch by falling through the comparison rather than by a `!(v > 0.0)`
        // that reads as if it meant `v <= 0.0` and does not.
        if s.volume > 0.0 {
            let out: f64 = self
                .links
                .iter()
                .filter(|l| l.a == i || l.b == i)
                .map(|l| l.q)
                .sum();
            (out + s.clearance) / s.volume
        } else {
            f64::INFINITY
        }
    }

    /// The fraction of the fastest compartment's contents its own outflow would remove in one
    /// step of `dt`. Must not exceed 1; see [`Domain::max_stable_dt`].
    ///
    /// The counterpart of a Fourier or a Courant number for this domain, and reported rather than
    /// only enforced so a caller can see how close a chosen step is before it is refused.
    pub fn turnover_number(&self, dt: Time) -> f64 {
        let h = dt.to_si();
        (0..self.compartments.len()).fold(0.0f64, |worst, i| worst.max(h * self.turnover_rate(i)))
    }
}

impl Domain for CompartmentModel {
    fn name(&self) -> &str {
        &self.name
    }

    fn kind(&self) -> Kind {
        Kind::Evolving
    }

    /// **Nothing crosses this domain's boundary**, so its books close against zero traffic and
    /// are checked on their own rather than only in the whole-simulation sum.
    ///
    /// The claim is exact and it is easy to say why: drug enters by a declared dose and leaves
    /// only into `eliminated`, which the ledger reports. There is no bus channel, so there is
    /// nothing for the difference to be attributed to and the audit is comparing the books
    /// against themselves one step apart. That is the strictest form the check has and this
    /// domain is the cheapest possible case of it.
    fn books_balance(&self) -> bool {
        true
    }

    /// `min_i V_i/(Σ_j Q_ij + CL_i)` — the reciprocal of the fastest compartment's turnover.
    ///
    /// # Where that comes from
    ///
    /// Write the system in concentrations, `dC/dt = N C`, with
    /// `N_ii = −(Σ_j Q_ij + CL_i)/V_i = −r_i` and `N_ij = Q_ij/V_i` for `j ≠ i`. Explicit Euler
    /// on that is stable exactly when `|1 + λ h| ≤ 1` for every eigenvalue, and since a
    /// compartment matrix is similar to a symmetric one its eigenvalues are real and
    /// non-positive, so the condition is `h ≤ 2/|λ|max`.
    ///
    /// Gershgorin bounds `|λ|max` without solving for it. Row `i` puts every eigenvalue within
    /// `Σ_{j≠i} Q_ij/V_i` of `−r_i`, and that radius is at most `r_i` — it is `r_i` less the
    /// clearance term. So `|λ|max ≤ 2 max_i r_i`, and `h ≤ 1/max_i r_i` implies `h ≤ 2/|λ|max`.
    /// That is the number reported.
    ///
    /// # And it is the positivity limit as well, which is the reason not to loosen it
    ///
    /// The same expression falls out of a different requirement. A compartment's update is
    /// `A_i ← A_i(1 − h r_i) + h Σ_j Q_ij C_j + dose`, and every term after the first is
    /// non-negative, so the amount stays non-negative exactly while `h r_i ≤ 1`. A negative mass
    /// of drug is not a large error, it is a meaningless one — and unlike a divergence it does
    /// not announce itself, because the system is still stable at `h` up to twice this.
    ///
    /// # What it costs
    ///
    /// For a single compartment the bound is conservative by a factor of two: `|λ| = k` exactly,
    /// so stability allows `2/k` and this reports `1/k`. The factor is what makes one expression
    /// safe for any number of compartments, and it is the price of not running an eigensolver
    /// every step. For two compartments of equal volume joined by `Q` and nothing else the bound
    /// is *tight* — `|λ|max = 2Q/V` and `2/|λ|max = V/Q = 1/r` — which is the case that shows the
    /// factor is not slack but the actual answer where the Gershgorin disc touches.
    ///
    /// Infinite for a model with no clearances and no links: nothing is moving, so nothing limits
    /// the step. Honest rather than defensive, and the same answer a still medium gives in
    /// `pantometry-acoustic`.
    fn max_stable_dt(&self, _now: Time) -> Time {
        let mut limit = f64::INFINITY;
        for i in 0..self.compartments.len() {
            let r = self.turnover_rate(i);
            if r > 0.0 {
                limit = limit.min(1.0 / r);
            }
        }
        Time::from_si(limit)
    }

    /// One explicit Euler step, with every flux taken from the same snapshot.
    ///
    /// Jacobi rather than Gauss-Seidel: the answer must not depend on the order the links were
    /// declared in, and updating in place would be a different scheme rather than a rounding
    /// difference.
    ///
    /// The bus is not touched. This domain neither produces nor consumes anything another domain
    /// could want — see the module docs — and publishing on a channel nobody consumes is refused
    /// by the bus's own audit, correctly.
    fn step(&mut self, t: Time, dt: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        let h = dt.to_si();
        if h <= 0.0 {
            return Ok(());
        }
        if self.compartments.is_empty() {
            return Err(Violation::at(
                self.name.clone(),
                "a model with no compartments",
                0.0,
            ));
        }

        // Refuse rather than diverge, and name the compartment. In a model with several the fast
        // one is not obvious — a small central volume on a large intercompartmental clearance is
        // far faster than the peripheral compartment it feeds — so the site says which it is.
        for i in 0..self.compartments.len() {
            let v = self.compartments[i].volume;
            if !v.is_finite() || v <= 0.0 {
                return Err(Violation::at(
                    format!("{}/{}", self.name, self.compartments[i].label),
                    "a compartment must have a positive volume",
                    v,
                ));
            }
            let turnover = h * self.turnover_rate(i);
            if turnover > 1.0 + 1e-12 {
                return Err(Violation {
                    quantity: "compartment turnover number".to_string(),
                    site: format!(
                        "{}/{} (explicit Euler)",
                        self.name, self.compartments[i].label
                    ),
                    before: 1.0,
                    after: turnover,
                    scale: 1.0,
                    tolerance: 1e-12,
                });
            }
        }

        let n = self.compartments.len();
        let conc: Vec<f64> = (0..n).map(|i| self.conc(i)).collect();
        let mut delta = vec![0.0; n];

        // Infusion into the same right-hand side as everything else. **The snapshot above is
        // what enforces that**, not the order of the loops: every flux below reads `conc`, taken
        // before anything moved, so the arriving dose cannot raise the concentration its own
        // compartment is then drained at.
        //
        // Getting it the other way round is the obvious reading and it is wrong in a way
        // conservation cannot see — every total stays exactly right and the steady state sits
        // low by exactly the turnover number `h·CL/V`. `pantometry-thermal` measured that same
        // mistake at 0.31% on a three-node ladder with the audit silent throughout; here it is
        // 0.5% at a 36 s step, and two tests below are shaped to see it.
        //
        // The window is intersected with the step exactly, so an infusion that starts or stops
        // part-way through a step delivers exactly its share and the answer does not depend on
        // the step happening to land on the boundary. The final slice takes the *remainder* of
        // the scheduled dose rather than `rate·span`, which is the same reason
        // `Exchange::take_share` hands the last substep whatever is left: `n` slices of a window
        // do not sum to the window in floating point, and the residue would sit in the ledger
        // forever.
        let now = t.to_si();
        let until = now + h;
        for inf in self.infusions.iter_mut() {
            let lo = now.max(inf.start);
            let hi = until.min(inf.end);
            // A positive guard, for the reason `turnover_rate` gives: a NaN bound leaves the
            // window closed instead of taking a branch whose spelling denies it.
            if hi > lo {
                let amount = if until >= inf.end {
                    (inf.total - inf.delivered).max(0.0)
                } else {
                    (inf.rate * (hi - lo))
                        .min(inf.total - inf.delivered)
                        .max(0.0)
                };
                inf.delivered += amount;
                delta[inf.into] += amount;
            }
        }

        // Each link's flux computed **once** and applied twice with opposite signs. Computing it
        // from each side separately gives two values differing in the last bit, and the model then
        // leaks about 1e-16 per link per step — invisible per step and not invisible over a run
        // that has to reach a terminal phase five half-lives out.
        for l in &self.links {
            let q = l.q * (conc[l.a] - conc[l.b]) * h;
            delta[l.a] -= q;
            delta[l.b] += q;
        }

        // Elimination. Accumulated with exactly the value subtracted, so what left the body and
        // what the books gained are the same number rather than two roundings of one.
        let mut cleared = 0.0;
        for (i, c) in self.compartments.iter().enumerate() {
            let out = c.clearance * conc[i] * h;
            delta[i] -= out;
            cleared += out;
        }
        self.eliminated += cleared;

        for (c, d) in self.compartments.iter_mut().zip(&delta) {
            c.amount += d;
        }
        Ok(())
    }

    /// Every compartment's holding, plus what has been cleared and what is still in the syringe.
    ///
    /// One `add` per compartment rather than one for the sum, so [`Ledger`]'s scale is the largest
    /// single holding instead of a near-zero net — which is what the scale exists for, and the
    /// mistake `NBody::ledger` made until a test proved its momentum audit was inert.
    ///
    /// Absolute masses rather than differences from a reference, unlike the thermal domains: a
    /// compartment starts empty and the natural datum is zero, so there is no rounding floor to
    /// avoid and the total *is* the dose, which is the number a reader wants the audit measured
    /// against.
    fn ledger(&self) -> Ledger {
        let mut ledger = Ledger::new();
        for c in &self.compartments {
            ledger.add(quantity::MASS, c.amount);
        }
        ledger.add(quantity::MASS, self.eliminated);
        ledger.add(quantity::MASS, self.pending_dose().to_si());
        ledger
    }

    /// Every amount, every infusion's progress, and both running totals.
    ///
    /// All of it, because `ledger` reads all of it. An infusion's `delivered` is the one easiest
    /// to forget and the one whose loss is silent: restore without it and the syringe refills,
    /// the books gain a dose, and the sweep that was meant to be re-run from the same state is
    /// re-run from a richer one.
    fn checkpoint(&mut self) {
        self.saved = Some((
            self.compartments.iter().map(|c| c.amount).collect(),
            self.infusions.iter().map(|i| i.delivered).collect(),
            self.eliminated,
            self.bolused,
        ));
    }

    fn restore(&mut self) {
        if let Some((amounts, delivered, eliminated, bolused)) = self.saved.clone() {
            for (c, a) in self.compartments.iter_mut().zip(amounts) {
                c.amount = a;
            }
            for (inf, d) in self.infusions.iter_mut().zip(delivered) {
                inf.delivered = d;
            }
            self.eliminated = eliminated;
            self.bolused = bolused;
        }
    }

    fn supports_restore(&self) -> bool {
        true
    }

    /// Every compartment's concentration, then the totals.
    ///
    /// **The whole output of this domain.** There is no field and there are no bodies, so a table
    /// is the only place a compartmental model appears at all — a domain that skipped this would
    /// run, conserve, and be silently absent from every report.
    ///
    /// Concentrations in **mg/L**, which is the unit a concentration-time curve is drawn in
    /// everywhere the subject is used, and equals µg/mL. SI would be kg/m³, a thousand times
    /// smaller; the label says which and a wrong label would be worse than none.
    fn readings(&self) -> Vec<Reading> {
        let mut out: Vec<Reading> = self
            .handles()
            .map(|(c, label)| {
                Reading::new(
                    &self.name,
                    label,
                    self.concentration(c).to_si() * 1e3,
                    "mg/L",
                )
            })
            .collect();
        out.push(Reading::new(
            &self.name,
            "in the body",
            self.body_burden().to_si() * 1e6,
            "mg",
        ));
        out.push(Reading::new(
            &self.name,
            "cleared",
            self.eliminated * 1e6,
            "mg",
        ));
        if !self.infusions.is_empty() {
            out.push(Reading::new(
                &self.name,
                "still to infuse",
                self.pending_dose().to_si() * 1e6,
                "mg",
            ));
        }
        out
    }

    fn as_any(&self) -> Option<&dyn Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn Any> {
        Some(self)
    }

    /// **`None`, and not as an oversight.**
    ///
    /// A compartment is not a place. The "central compartment" is plasma plus every tissue that
    /// equilibrates with it faster than a blood sample can resolve, and its volume is *apparent*
    /// — it is whatever number makes `A/V` match the measured concentration, and for a drug that
    /// partitions into fat it routinely exceeds the volume of the person. There is no space to
    /// sample and interpolating between two compartments would invent one, with less
    /// justification than a thermal network has for declining the same thing.
    fn as_field(&self) -> Option<&dyn pantometry_core::ScalarField> {
        None
    }

    /// `None` for the same reason. Compartments are countable, but they are not *at* anywhere,
    /// and [`Bodies`](pantometry_core::Bodies) is about things at places.
    fn as_bodies(&self) -> Option<&dyn pantometry_core::Bodies> {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use pantometry_core::{Schedule, Simulation};

    /// Seconds in an hour. Every number in pharmacokinetics is quoted per hour and every number
    /// in this workspace is stored in seconds, so the conversion is written once.
    const HOUR: f64 = 3600.0;

    fn hours(h: f64) -> Time {
        Time::s(h * HOUR)
    }

    /// A one-compartment model: 10 L, 5 L/h, so `k = 0.5/h`.
    fn one_compartment() -> (CompartmentModel, Compartment) {
        let mut m = CompartmentModel::new("subject");
        let c = m.eliminating("plasma", Volume::litres(10.0), Clearance::l_per_h(5.0));
        (m, c)
    }

    /// A two-compartment model with the micro constants worked out in the comment, so the closed
    /// form below is built from numbers a reader can check by hand.
    ///
    /// `V₁ = 5 L`, `V₂ = 20 L`, `Q = 30 L/h`, `CL = 10 L/h`:
    /// `k10 = 2/h`, `k12 = 6/h`, `k21 = 1.5/h`.
    fn two_compartment() -> (CompartmentModel, Compartment, Compartment) {
        let mut m = CompartmentModel::new("subject");
        let central = m.eliminating("central", Volume::litres(5.0), Clearance::l_per_h(10.0));
        let peripheral = m.compartment("peripheral", Volume::litres(20.0));
        m.link(central, peripheral, Clearance::l_per_h(30.0))
            .unwrap();
        (m, central, peripheral)
    }

    fn run(model: &mut CompartmentModel, from: Time, dt: Time, steps: usize) {
        let mut bus = Exchange::new();
        let mut t = from;
        for _ in 0..steps {
            model.step(t, dt, &mut bus).unwrap();
            t += dt;
        }
    }

    /// **Closed form 1: `A(t) = A₀ e^{−kt}` with `k = CL/V`.**
    ///
    /// The tolerance is not chosen. Explicit Euler on `A' = −kA` gives `A₀(1 − kh)ⁿ`, and
    /// `ln[(1 − x)ⁿ] = −nx − nx²/2 − …`, so the stepped answer is low by a factor
    /// `e^{−n x²/2} ≈ 1 − t k² h/2` — first order in `h` and growing linearly in `t`. Every
    /// tolerance below is `1.1×` that prediction, which means the test fails both if the code is
    /// wrong and if the error is *smaller* than first-order Euler can be, i.e. if the step being
    /// taken is not the step being asked for.
    #[test]
    fn a_bolus_decays_at_the_rate_the_clearance_sets() {
        let (mut m, plasma) = one_compartment();
        let dose = 1e-4; // 100 mg
        m.bolus(plasma, Mass::from_si(dose)).unwrap();

        let k = m.elimination_rate(plasma).to_si();
        assert!(
            (k * HOUR - 0.5).abs() < 1e-15,
            "k should be 0.5/h, got {}",
            k * HOUR
        );

        let h = 36.0; // seconds; k h = 5e-3, well inside the limit of 1/k = 2 h
        let mut t = 0.0;
        for target in [1.0, 2.0, 4.0, 8.0] {
            let want_at = target * HOUR;
            let steps = ((want_at - t) / h).round() as usize;
            run(&mut m, Time::s(t), Time::s(h), steps);
            t += steps as f64 * h;

            let exact = dose * (-k * t).exp();
            let got = m.amount(plasma).to_si();
            let predicted = t * k * k * h / 2.0;
            assert!(
                (1.0 - got / exact) > 0.0,
                "explicit Euler undershoots this; at {target} h it gave {got:e} against {exact:e}"
            );
            assert!(
                (got / exact - 1.0).abs() < 1.1 * predicted,
                "at {target} h: {got:e} against a closed form of {exact:e}, a relative error of \
                 {:e} where first-order Euler predicts {predicted:e}",
                (got / exact - 1.0).abs()
            );
        }

        // Eight hours is four elimination time constants; 1.83% of the dose is left. Relative
        // to the closed form and against the same `t k² h/2`, which is 1.0e-2 here — an absolute
        // 1e-4 on a value of 1.8e-2 would have been a 0.5% tolerance on a 1% error, and was.
        let left = m.amount(plasma).to_si() / dose;
        let exact = (-4.0f64).exp();
        assert!(
            (left / exact - 1.0).abs() < 1.1 * (8.0 * HOUR) * k * k * h / 2.0,
            "at eight hours {left:e} against {exact:e}"
        );
    }

    /// **Closed form 2: `A(t) = (R/k)(1 − e^{−kt})`, steady state `R·V/CL`.**
    ///
    /// Two claims with two very different tolerances, and the difference is the point.
    ///
    /// The steady state is checked at `1e-12`, which looks unearned and is not: the fixed point
    /// of the Euler map `A ← A + h(R − kA)` is `A* = R/k` **exactly**, at any step size, because
    /// the discretisation error is in how the solution approaches the fixed point and not in
    /// where the fixed point is. What is left is float accumulation over the steps taken.
    ///
    /// The transient is checked against `t k² h/2` as above, because that part *is* first order.
    #[test]
    fn a_constant_infusion_climbs_to_the_steady_state_the_clearance_sets() {
        let (mut m, plasma) = one_compartment();
        let k = m.elimination_rate(plasma).to_si();
        let dose = 1e-3; // 1 g over 100 hours
        let window = 100.0 * HOUR;
        let rate = dose / window;
        m.infuse(plasma, Mass::from_si(dose), Time::ZERO, Time::s(window))
            .unwrap();
        // The whole dose is on the books from the moment it is scheduled.
        assert_eq!(m.pending_dose().to_si(), dose);

        let h = 36.0;
        // One time constant in: the exponential is 63% of the way there.
        let t1 = 2.0 * HOUR; // k t = 1
        run(&mut m, Time::ZERO, Time::s(h), (t1 / h) as usize);
        let exact = (rate / k) * (1.0 - (-k * t1).exp());
        let got = m.amount(plasma).to_si();
        let predicted = t1 * k * k * h / 2.0;
        assert!(
            (got / exact - 1.0).abs() < 1.1 * predicted,
            "at one time constant: {got:e} against {exact:e}"
        );

        // Then out to 80 hours, which is 40 time constants: e^-40 is 4e-18 and the transient is
        // gone below the floating-point floor, so what is left is the fixed point.
        let steps = ((80.0 * HOUR - t1) / h).round() as usize;
        run(&mut m, Time::s(t1), Time::s(h), steps);
        let ss = rate * 0.01 / (5e-3 / 3600.0); // R V / CL, spelled out
        let want = Volume::litres(10.0).to_si() * rate / Clearance::l_per_h(5.0).to_si();
        assert!((ss / want - 1.0).abs() < 1e-12, "the two spellings agree");
        let got = m.amount(plasma).to_si();
        assert!(
            (got / want - 1.0).abs() < 1e-12,
            "the Euler fixed point is the exact steady state: {got:e} against {want:e}"
        );
        // And the concentration form of the same statement, R/CL.
        let c_ss = m.concentration(plasma).to_si();
        assert!((c_ss / (rate / Clearance::l_per_h(5.0).to_si()) - 1.0).abs() < 1e-12);
    }

    /// **Closed form 3: the bi-exponential, and the transposed index it catches.**
    ///
    /// With `α + β = k10 + k12 + k21` and `αβ = k10 k21`,
    ///
    /// ```text
    /// A₁(t) = D [ (α − k21)/(α − β) e^{−αt} + (k21 − β)/(α − β) e^{−βt} ]
    /// A₂(t) = D k12/(α − β) [ e^{−βt} − e^{−αt} ]
    /// ```
    ///
    /// The conservation audit cannot see a link at all — `+q` and `−q` are in the same sum — so
    /// this is the check that a link moves the drug in the right direction, in the right amount,
    /// and between the right two compartments. Swapping `k12` and `k21` leaves the *sum* of the
    /// two compartments exactly right and both curves wrong.
    ///
    /// The tolerance is earned by measuring rather than by asserting: the same comparison is run
    /// at `h` and at `h/2`, and the error is required to **halve**. A wrong scheme that happened
    /// to be within a fixed tolerance would not also halve.
    #[test]
    fn two_compartments_follow_the_bi_exponential() {
        // The micro constants, per second.
        let k10 = 2.0 / HOUR;
        let k12 = 6.0 / HOUR;
        let k21 = 1.5 / HOUR;
        let sum = k10 + k12 + k21;
        let disc = (sum * sum - 4.0 * k10 * k21).sqrt();
        let alpha = 0.5 * (sum + disc);
        let beta = 0.5 * (sum - disc);
        assert!((alpha * HOUR - 9.17298).abs() < 1e-4, "{}", alpha * HOUR);
        assert!((beta * HOUR - 0.32702).abs() < 1e-4, "{}", beta * HOUR);

        let dose = 1e-4;
        let central_exact = |t: f64| {
            dose * ((alpha - k21) / (alpha - beta) * (-alpha * t).exp()
                + (k21 - beta) / (alpha - beta) * (-beta * t).exp())
        };
        let peripheral_exact =
            |t: f64| dose * k12 / (alpha - beta) * ((-beta * t).exp() - (-alpha * t).exp());

        // Sanity on the closed form itself, independent of any code under test.
        assert!((central_exact(0.0) / dose - 1.0).abs() < 1e-14);
        assert!(peripheral_exact(0.0).abs() < 1e-20);

        let sample_at = [0.1, 0.5, 1.0, 3.0];
        let mut errors = Vec::new();
        for h in [18.0, 9.0] {
            let (mut m, central, peripheral) = two_compartment();
            m.bolus(central, Mass::from_si(dose)).unwrap();
            let mut worst: f64 = 0.0;
            let mut t = 0.0;
            for target in sample_at {
                let want_at = target * HOUR;
                let steps = ((want_at - t) / h).round() as usize;
                run(&mut m, Time::s(t), Time::s(h), steps);
                t += steps as f64 * h;

                for (got, exact, which) in [
                    (m.amount(central).to_si(), central_exact(t), "central"),
                    (
                        m.amount(peripheral).to_si(),
                        peripheral_exact(t),
                        "peripheral",
                    ),
                ] {
                    // Relative to the dose, not to the value: the central compartment falls
                    // through four orders of magnitude here and a relative-to-itself measure
                    // would make the tail dominate a comparison it should not.
                    let err = (got - exact).abs() / dose;
                    assert!(
                        err < 0.02,
                        "{which} at {target} h with h = {h} s: {got:e} against {exact:e}"
                    );
                    worst = worst.max(err);
                }
            }
            errors.push(worst);
        }
        // First order: halving the step halves the error. 2.0 within 10%, which is the residual
        // of the O(h²) term at these step sizes.
        let ratio = errors[0] / errors[1];
        assert!(
            (ratio - 2.0).abs() < 0.2,
            "the error should halve with the step; it went {:e} -> {:e}, a ratio of {ratio:.3}",
            errors[0],
            errors[1]
        );
    }

    /// **The per-link check the audit is blind to.** One step from a bolus, and the amount that
    /// crossed is `Q·C₁·h` to the last bit — out of one compartment and into the other.
    ///
    /// A transposed index here passes every conservation check ever written, because the two
    /// halves are in the same sum. This is the smallest statement that would notice.
    #[test]
    fn a_link_moves_exactly_what_the_concentration_difference_says() {
        let (mut m, central, peripheral) = two_compartment();
        let dose = 1e-4;
        m.bolus(central, Mass::from_si(dose)).unwrap();

        let q = Clearance::l_per_h(30.0).to_si();
        let cl = Clearance::l_per_h(10.0).to_si();
        let c1 = m.concentration(central).to_si();
        assert!((c1 - dose / Volume::litres(5.0).to_si()).abs() < 1e-18);
        // The rate, before anything has moved, reported by the domain.
        assert!((m.transfer(central, peripheral).to_si() - q * c1).abs() < 1e-24);
        // And the other way round is the negation, not a second value.
        assert_eq!(
            m.transfer(peripheral, central).to_si(),
            -m.transfer(central, peripheral).to_si()
        );

        let h = 10.0;
        run(&mut m, Time::ZERO, Time::s(h), 1);
        let gained = m.amount(peripheral).to_si();
        assert!(
            (gained - q * c1 * h).abs() < 1e-22,
            "the peripheral compartment should have gained exactly Q C1 h = {:e}, got {gained:e}",
            q * c1 * h
        );
        // `dose − amount` is a difference of two numbers near 1e-4 whose difference is near
        // 2e-6, so it carries `ε·dose ≈ 2e-20` of cancellation — two orders of magnitude above
        // the 1e-22 the gain above is checked to, and the reason the two tolerances differ.
        let lost = dose - m.amount(central).to_si();
        assert!(
            (lost - (q + cl) * c1 * h).abs() < 1e-19,
            "the central one should have lost the link and the clearance together: {lost:e}              against {:e}",
            (q + cl) * c1 * h
        );
        // Nothing crossed a link that does not exist.
        let mut lonely = CompartmentModel::new("lonely");
        let a = lonely.compartment("a", Volume::litres(1.0));
        let b = lonely.compartment("b", Volume::litres(1.0));
        assert_eq!(lonely.transfer(a, b).to_si(), 0.0);
        assert_eq!(lonely.transfer_rate(a, b).to_si(), 0.0);
    }

    /// **Closed form 4: dose = in the compartments + cleared + still in the syringe.**
    ///
    /// To machine precision, and with a *scale*: the residual is measured against the dose rather
    /// than against itself, because the books cancelling to nothing is what a correct model looks
    /// like and a relative tolerance on the net would be meaningless.
    ///
    /// The bound is `1e-12` and it traces to the step count. Every step does `O(n)` additions on
    /// the amounts, so the total drifts by about `ε` relative per step; 12 000 steps at
    /// `ε = 2.2e-16` is `2.6e-12` in the worst case and `~1e-14` in practice, since the errors are
    /// not all the same sign.
    ///
    /// Everything is exercised at once deliberately — two compartments, a link, clearance from
    /// both of them, a bolus and an infusion whose window neither starts nor ends on a step
    /// boundary — because mass balance is the one check that gets *easier* the simpler the model.
    #[test]
    fn the_dose_is_conserved_between_the_body_and_the_books() {
        let mut m = CompartmentModel::new("subject");
        let central = m.eliminating("central", Volume::litres(5.0), Clearance::l_per_h(10.0));
        let peripheral = m.eliminating(
            "peripheral",
            Volume::litres(20.0),
            Clearance::l_per_h(1.0), // a second elimination path, so the totals are not trivial
        );
        m.link(central, peripheral, Clearance::l_per_h(30.0))
            .unwrap();

        let bolus = 1e-4;
        let infused = 2.5e-4;
        m.bolus(central, Mass::from_si(bolus)).unwrap();
        // A window on neither a step boundary nor a whole hour: 0.37 h to 2.81 h.
        m.infuse(peripheral, Mass::from_si(infused), hours(0.37), hours(2.81))
            .unwrap();
        let dose = bolus + infused;

        let h = 3.0; // seconds
        let steps = (12.0 * HOUR / h) as usize;
        let mut bus = Exchange::new();
        let mut t = Time::ZERO;
        let mut worst: f64 = 0.0;
        for step in 0..steps {
            m.step(t, Time::s(h), &mut bus).unwrap();
            t += Time::s(h);
            let books =
                m.body_burden().to_si() + m.eliminated_mass().to_si() + m.pending_dose().to_si();
            worst = worst.max((books - dose).abs() / dose);
            // Half-way through, the infusion is over and the syringe is empty to the last bit.
            if step + 1 == steps {
                assert_eq!(m.pending_dose().to_si(), 0.0);
                assert!(
                    (m.administered_mass().to_si() / dose - 1.0).abs() < 1e-15,
                    "the whole dose was given: {:e}",
                    m.administered_mass().to_si()
                );
            }
        }
        assert!(
            worst < 1e-12,
            "mass balance drifted by {worst:e} of the dose over {steps} steps"
        );
        // And it is not vacuous: drug really did move and really was cleared.
        assert!(m.eliminated_mass().to_si() > 0.5 * dose);
        assert!(m.amount(peripheral).to_si() > 0.0);
    }

    /// **An infusion delivers its whole dose whatever the step lands on.**
    ///
    /// The window is intersected with each step exactly, so a boundary part-way through a step is
    /// apportioned rather than rounded to the nearest one. Three step sizes, none of which divide
    /// the window, and all three deliver the same total to the last bit — which a scheme that
    /// snapped the window to step boundaries would fail by up to a step's worth.
    #[test]
    fn an_infusion_delivers_its_dose_whatever_the_step_divides() {
        let infused = 3e-4;
        for h in [7.0, 13.0, 101.0] {
            let (mut m, plasma) = one_compartment();
            m.infuse(
                plasma,
                Mass::from_si(infused),
                Time::s(1234.5),
                Time::s(1234.5 + 1800.0),
            )
            .unwrap();
            run(&mut m, Time::ZERO, Time::s(h), (6.0 * HOUR / h) as usize);
            assert_eq!(
                m.pending_dose().to_si(),
                0.0,
                "the syringe should be empty at h = {h}"
            );
            assert!(
                (m.administered_mass().to_si() / infused - 1.0).abs() < 1e-15,
                "h = {h} delivered {:e} of {infused:e}",
                m.administered_mass().to_si()
            );
            // Nothing arrived before the infusion started.
            let mut early = CompartmentModel::new("early");
            let c = early.eliminating("c", Volume::litres(10.0), Clearance::l_per_h(5.0));
            early
                .infuse(c, Mass::from_si(infused), hours(2.0), hours(3.0))
                .unwrap();
            run(&mut early, Time::ZERO, Time::s(h), (HOUR / h) as usize);
            assert_eq!(early.amount(c).to_si(), 0.0);
            assert_eq!(early.pending_dose().to_si(), infused);
        }
    }

    /// **Closed form 5: explicit Euler is first order, measured as a rate.**
    ///
    /// Four step sizes each half the last, against the exact exponential. The ratio of successive
    /// errors is asserted to be 2 within 3%, which is what distinguishes "first order" from "some
    /// error that happens to be small": a second-order scheme would give 4, and a scheme with a
    /// constant bias would give 1. `CONTRIBUTING.md` records that this is the check that found
    /// the acoustic boundary defect, where the value was plausible and the *rate* was not.
    #[test]
    fn the_error_is_first_order_in_the_step() {
        let dose = 1e-4;
        let target = 2.0 * HOUR;
        let mut errors = Vec::new();
        for h in [240.0, 120.0, 60.0, 30.0] {
            let (mut m, plasma) = one_compartment();
            m.bolus(plasma, Mass::from_si(dose)).unwrap();
            let k = m.elimination_rate(plasma).to_si();
            let steps = (target / h) as usize;
            run(&mut m, Time::ZERO, Time::s(h), steps);
            let exact = dose * (-k * target).exp();
            errors.push((m.amount(plasma).to_si() - exact).abs() / exact);
        }
        for w in errors.windows(2) {
            let ratio = w[0] / w[1];
            assert!(
                (ratio - 2.0).abs() < 0.06,
                "halving the step should halve the error; {:e} -> {:e} is a ratio of {ratio:.4}",
                w[0],
                w[1]
            );
        }
        // And the errors are real rather than at the floating-point floor, which is the way a
        // convergence test most often passes without measuring anything.
        assert!(
            *errors.last().unwrap() > 1e-6,
            "the finest step still has a measurable error: {:e}",
            errors.last().unwrap()
        );
    }

    /// **The step limit is a bound on the real eigenvalue, not a guess beside it.**
    ///
    /// Checked three ways. Against the expression it claims to be; against the *actual* largest
    /// eigenvalue of the two-compartment matrix, which is `α` and available in closed form, so
    /// the Gershgorin bound is verified rather than trusted; and against the case where the bound
    /// is tight, two equal compartments joined by `Q`, where `2/|λ|max` is exactly what is
    /// reported.
    #[test]
    fn the_step_limit_is_derived_and_bounds_the_true_one() {
        // One compartment: r = k, so the limit is 1/k = 2 h. True stability is 2/k = 4 h, and the
        // factor of two is the price of a bound that holds for any number of compartments.
        let (one, plasma) = one_compartment();
        let k = one.elimination_rate(plasma).to_si();
        let limit = one.max_stable_dt(Time::ZERO).to_si();
        assert!((limit * k - 1.0).abs() < 1e-12, "{} h", limit / HOUR);
        assert!(
            limit <= 2.0 / k + 1e-9,
            "and it is inside the true 2/|lambda|"
        );

        // Two compartments: r1 = (Q + CL)/V1 = 8/h, r2 = Q/V2 = 1.5/h. The central one binds.
        let (two, _, _) = two_compartment();
        let limit = two.max_stable_dt(Time::ZERO).to_si();
        assert!(
            (limit / HOUR - 0.125).abs() < 1e-12,
            "1/8 of an hour, got {} h",
            limit / HOUR
        );
        // Against the true eigenvalue, from the characteristic polynomial.
        let (k10, k12, k21) = (2.0 / HOUR, 6.0 / HOUR, 1.5 / HOUR);
        let sum = k10 + k12 + k21;
        let alpha = 0.5 * (sum + (sum * sum - 4.0 * k10 * k21).sqrt());
        assert!(
            limit <= 2.0 / alpha,
            "the reported limit {limit} must not exceed the true 2/alpha = {}",
            2.0 / alpha
        );
        // And it is a bound rather than a shrug: within a factor of two of the true one.
        assert!(limit > 1.0 / alpha, "not needlessly conservative");

        // Where the disc touches: equal volumes, one link, no clearance. |lambda|max = 2Q/V and
        // the reported limit is exactly 2/|lambda|max.
        let mut tight = CompartmentModel::new("tight");
        let a = tight.compartment("a", Volume::litres(4.0));
        let b = tight.compartment("b", Volume::litres(4.0));
        tight.link(a, b, Clearance::l_per_h(2.0)).unwrap();
        let true_limit =
            2.0 / (2.0 * Clearance::l_per_h(2.0).to_si() / Volume::litres(4.0).to_si());
        assert!(
            (tight.max_stable_dt(Time::ZERO).to_si() / true_limit - 1.0).abs() < 1e-12,
            "the bound is exact here: {} against {true_limit}",
            tight.max_stable_dt(Time::ZERO).to_si()
        );

        // Nothing moving, no limit.
        let mut still = CompartmentModel::new("still");
        still.compartment("a", Volume::litres(1.0));
        assert!(!still.max_stable_dt(Time::ZERO).to_si().is_finite());
    }

    /// **Past the limit the step is refused, naming the compartment and by how much.**
    ///
    /// Not degraded, not clamped. A caller that asked for too much gets a [`Violation`] saying
    /// which compartment could not take it, which matters because in a model with several the
    /// fast one is not the obvious one — here it is the small central volume, not the peripheral
    /// compartment that holds four times as much.
    #[test]
    fn a_step_past_the_limit_is_refused() {
        let (mut m, central, _) = two_compartment();
        m.bolus(central, Mass::from_si(1e-4)).unwrap();
        let limit = m.max_stable_dt(Time::ZERO);
        assert!((m.turnover_number(limit) - 1.0).abs() < 1e-12);

        let mut bus = Exchange::new();
        assert!(m.step(Time::ZERO, limit, &mut bus).is_ok());
        let err = m
            .step(Time::ZERO, limit * 1.05, &mut bus)
            .expect_err("past the turnover limit must be refused");
        assert_eq!(err.quantity, "compartment turnover number");
        assert!(err.site.contains("central"), "{}", err.site);
        assert!((err.after - 1.05).abs() < 1e-9, "{}", err.after);

        // A zero or negative step changes nothing rather than being an error.
        let before = m.amount(central);
        m.step(Time::ZERO, Time::ZERO, &mut bus).unwrap();
        assert_eq!(m.amount(central), before);
    }

    /// **No compartment ever holds a negative amount of drug, right up to the limit.**
    ///
    /// The other half of what `max_stable_dt` is derived from, and the half that does not
    /// announce itself: the scheme is *stable* to twice this step, so a model run between the two
    /// oscillates about the right answer with negative concentrations in it and never diverges.
    /// Run at exactly the limit, where the diagonal factor `1 − h r` is exactly zero.
    #[test]
    fn amounts_stay_non_negative_at_the_limit() {
        let (mut m, central, peripheral) = two_compartment();
        m.bolus(central, Mass::from_si(1e-4)).unwrap();
        let dt = m.max_stable_dt(Time::ZERO);
        let mut bus = Exchange::new();
        let mut t = Time::ZERO;
        for step in 0..400 {
            m.step(t, dt, &mut bus).unwrap();
            t += dt;
            for (c, label) in m.handles() {
                assert!(
                    m.amount(c).to_si() >= 0.0,
                    "{label} went negative at step {step}: {:e}",
                    m.amount(c).to_si()
                );
            }
        }
        // Fifty hours out, essentially all of it has been cleared.
        assert!(m.eliminated_mass().to_si() / 1e-4 > 0.999);
        assert!(m.amount(peripheral).to_si() < 1e-7 * 1e-4);
    }

    /// **The domain runs under the scheduler and its books are checked on their own.**
    ///
    /// `books_balance` is the strict per-domain audit rather than the whole-simulation sum, and
    /// this domain is the cheapest case of it: no bus traffic at all, so the ledger must be
    /// constant across a step to the last bit the arithmetic allows.
    #[test]
    fn the_domain_runs_under_the_scheduler() {
        let (mut m, central, peripheral) = two_compartment();
        m.bolus(central, Mass::from_si(1e-4)).unwrap();
        // An infusion as well as a bolus, deliberately: the pending dose is a ledger entry, and
        // a ledger that dropped it would grow by the whole infusion over the run while every
        // other test here still passed. This is the only thing that reads that entry.
        m.infuse(peripheral, Mass::from_si(2e-4), Time::s(600.0), hours(4.0))
            .unwrap();
        assert!(m.books_balance());

        let mut sim = Simulation::new(Schedule::Multirate)
            // Nothing crosses the boundary and nothing accumulates but rounding, so this can be
            // as tight as the arithmetic allows over the substeps of one window.
            .conservation_tolerance(1e-12)
            .with(m);

        // A ten-minute window at a 450 s limit is two substeps.
        let report = sim.advance(Time::s(600.0)).unwrap();
        assert_eq!(report.substeps[0].0, "subject");
        assert_eq!(report.substeps[0].1, 2);
        for _ in 0..100 {
            sim.advance(Time::s(600.0)).expect("the dose should hold");
        }

        // And the concrete type comes back out, which is what `as_any` is for. Four mechanics
        // domains skipped it once and a renderer drew an empty frame with no error anywhere.
        let body: &CompartmentModel = sim.domain_as("subject").expect("as_any is implemented");
        let plasma = body.compartment_named("central").unwrap();
        assert!(body.amount(plasma).to_si() < 1e-4);
        assert_eq!(body.compartments(), 2);
        // 101 windows of 600 s is 16.8 h, so the infusion finished and its books closed.
        assert_eq!(body.pending_dose().to_si(), 0.0);
        assert!(body.eliminated_mass().to_si() / 3e-4 > 0.9);
    }

    /// **Restore puts back everything the ledger reads, the syringe included.**
    ///
    /// An infusion's `delivered` is the field easiest to leave out of a checkpoint and the one
    /// whose loss is silent: restore without it and the pending dose goes back up while the
    /// compartment keeps what it was given, so the books gain a dose out of nowhere.
    #[test]
    fn a_checkpoint_carries_the_infusion_as_well_as_the_amounts() {
        let (mut m, plasma) = one_compartment();
        m.infuse(plasma, Mass::from_si(1e-4), Time::ZERO, hours(4.0))
            .unwrap();
        assert!(m.supports_restore());
        run(&mut m, Time::ZERO, Time::s(60.0), 30);

        m.checkpoint();
        let books = m.ledger();
        let (amount, pending, cleared) = (
            m.amount(plasma).to_si(),
            m.pending_dose().to_si(),
            m.eliminated_mass().to_si(),
        );
        assert!(
            pending > 0.0 && cleared > 0.0,
            "the test needs both to move"
        );

        run(&mut m, Time::s(1800.0), Time::s(60.0), 30);
        assert!(m.pending_dose().to_si() < pending);
        m.restore();

        assert_eq!(m.amount(plasma).to_si(), amount);
        assert_eq!(m.pending_dose().to_si(), pending);
        assert_eq!(m.eliminated_mass().to_si(), cleared);
        assert_eq!(m.ledger(), books);
    }

    /// **The readings are the whole output.** No field and no bodies, so a table is where this
    /// domain appears at all — and the units on it have to be the ones the values are in.
    #[test]
    fn the_readings_name_every_compartment_and_the_totals() {
        let (mut m, central, _) = two_compartment();
        m.bolus(central, Mass::from_si(1e-4)).unwrap();
        assert!(m.as_field().is_none(), "a compartment is not a place");
        assert!(m.as_bodies().is_none(), "and it is not at one either");

        let r = m.readings();
        assert_eq!(r.len(), 4, "two compartments, the burden and the cleared");
        assert_eq!(r[0].label, "central");
        assert_eq!(r[1].label, "peripheral");
        assert_eq!(r[0].unit, "mg/L");
        // 100 mg into 5 L is 20 mg/L, and the label says mg/L.
        assert!((r[0].value - 20.0).abs() < 1e-9, "{}", r[0].value);
        assert_eq!(r[1].value, 0.0);
        assert_eq!(r[2].label, "in the body");
        assert_eq!(r[2].unit, "mg");
        assert!((r[2].value - 100.0).abs() < 1e-9, "{}", r[2].value);
        assert_eq!(r[3].label, "cleared");
        assert!(r.iter().all(|x| x.domain == "subject"));

        // A boxed domain still reports them: `Simulation` holds `Box<dyn Domain>` and a box that
        // swallowed the readings would make them vanish for exactly the caller who needs them.
        let direct = m.readings();
        let boxed: Box<dyn Domain> = Box::new(m);
        assert_eq!(direct, boxed.readings());

        // With an infusion there is a fourth total, and without one there is not — a column of
        // zeroes for a model with no pump is noise in a table.
        let (mut infusing, plasma) = one_compartment();
        infusing
            .infuse(plasma, Mass::from_si(1e-4), Time::ZERO, hours(1.0))
            .unwrap();
        let reported = infusing.readings();
        let labels: Vec<&str> = reported.iter().map(|r| &*r.label).collect();
        assert!(labels.contains(&"still to infuse"), "{labels:?}");
    }

    /// **A bad model is refused where it is built, not where it runs.**
    ///
    /// A handle from another model, a self-link, a negative clearance, a backwards infusion
    /// window, an unbounded one, and a compartment with no volume. Each is named.
    #[test]
    fn nonsense_is_refused_and_named() {
        let (mut m, central, peripheral) = two_compartment();
        // A *differently named* model. Identity is a hash of the name, deterministically and
        // with no counter, so two models called the same thing are the same model as far as a
        // handle can tell. This test tripped over that before the type said so — both fixtures
        // were called `"subject"`, and the foreign handle resolved to index 0 of this model and
        // was rejected as a self-link instead. Nothing above catches a duplicate name either:
        // `Simulation` accepts two domains with one, which was checked rather than assumed.
        let mut other = CompartmentModel::new("another subject");
        let elsewhere = other.eliminating("plasma", Volume::litres(10.0), Clearance::l_per_h(5.0));

        let err = m
            .link(central, elsewhere, Clearance::l_per_h(1.0))
            .expect_err("a handle from another model");
        assert!(err.quantity.contains("different model"), "{err}");
        assert!(m.bolus(elsewhere, Mass::from_si(1.0)).is_err());
        assert!(other.bolus(central, Mass::from_si(1.0)).is_err());

        let err = m
            .link(central, central, Clearance::l_per_h(1.0))
            .expect_err("a self-link");
        assert!(err.quantity.contains("itself"), "{err}");
        assert!(m
            .link(central, peripheral, Clearance::from_si(-1.0))
            .is_err());
        assert!(m
            .link(central, peripheral, Clearance::from_si(f64::NAN))
            .is_err());
        assert!(m.bolus(central, Mass::from_si(-1e-6)).is_err());
        // A zero dose is not an error; a regimen with a skipped dose should not need a branch.
        assert!(m.bolus(central, Mass::ZERO).is_ok());

        assert!(m
            .infuse(central, Mass::from_si(1e-4), hours(3.0), hours(1.0))
            .is_err());
        assert!(m
            .infuse(central, Mass::from_si(1e-4), hours(1.0), hours(1.0))
            .is_err());
        let err = m
            .infuse(
                central,
                Mass::from_si(1e-4),
                Time::ZERO,
                Time::from_si(f64::INFINITY),
            )
            .expect_err("an infusion with no end has an unbounded reserve");
        assert!(err.quantity.contains("finite"), "{err}");

        // A compartment with no volume is caught by `step`, which is the first moment it matters.
        let mut hollow = CompartmentModel::new("hollow");
        hollow.compartment("nothing", Volume::ZERO);
        let err = hollow
            .step(Time::ZERO, Time::s(1.0), &mut Exchange::new())
            .expect_err("no volume, no concentration");
        assert!(err.quantity.contains("positive volume"), "{err}");
        assert_eq!(hollow.max_stable_dt(Time::ZERO).to_si(), 0.0);

        // And a model with nothing in it says so rather than stepping vacuously.
        let mut empty = CompartmentModel::new("empty");
        assert!(empty
            .step(Time::ZERO, Time::s(1.0), &mut Exchange::new())
            .is_err());
    }

    /// **The derived quantities agree with the definitions they are derived from.**
    ///
    /// `Vss = ΣV`, `CL_total = ΣCL`, and the steady-state identity `Σ CL_i C_i = R` — which for
    /// elimination from one compartment says its concentration settles at `R/CL` *whatever* the
    /// intercompartmental clearance is. Run with two very different `Q` to show it does not
    /// depend on the distribution at all.
    #[test]
    fn the_steady_state_does_not_depend_on_the_distribution() {
        let dose = 5e-4;
        let window = 200.0 * HOUR;
        let rate = dose / window;
        let mut settled = Vec::new();
        for q in [3.0, 300.0] {
            let mut m = CompartmentModel::new("subject");
            let central = m.eliminating("central", Volume::litres(5.0), Clearance::l_per_h(10.0));
            let peripheral = m.compartment("peripheral", Volume::litres(20.0));
            m.link(central, peripheral, Clearance::l_per_h(q)).unwrap();
            m.infuse(central, Mass::from_si(dose), Time::ZERO, Time::s(window))
                .unwrap();

            assert!((m.volume_of_distribution().in_litres() - 25.0).abs() < 1e-12);
            assert!((m.total_clearance().in_l_per_h() - 10.0).abs() < 1e-12);

            let dt = m.max_stable_dt(Time::ZERO);
            run(&mut m, Time::ZERO, dt, (150.0 * HOUR / dt.to_si()) as usize);
            settled.push(m.concentration(central).to_si());
        }
        let want = rate / Clearance::l_per_h(10.0).to_si();
        for (c, q) in settled.iter().zip([3.0, 300.0]) {
            assert!(
                (c / want - 1.0).abs() < 1e-6,
                "at Q = {q} L/h the plateau was {c:e} against R/CL = {want:e}"
            );
        }
        // A hundredfold change in Q moves the plateau by nothing.
        assert!((settled[0] / settled[1] - 1.0).abs() < 1e-6);
    }
}
