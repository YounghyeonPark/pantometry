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
//! A tissue perfused from a blood pool is the one exception to that symmetry — it exchanges
//! `Q·(C_b − C_t/K_p)` and clears on `C_t/K_p` — and has a section of its own below.
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
//! - **One compartment, oral dose.** The Bateman function, its peak, its `kₐ = k` limit, the
//!   flip-flop terminal slope `−kₐ` when absorption is the slower process, and `AUC = F·D/CL`.
//! - **Two compartments, oral dose.** The tri-exponential, which catches a depot feeding the
//!   wrong compartment the way the bi-exponential catches a transposed link.
//! - **Perfusion-limited tissues.** Under a constant infusion every tissue settles at `K_p`
//!   times the blood and the blood at `R/(CL_r + Q_h E)`, and Euler's fixed point *is* that
//!   point; `V_ss = V_b + Σ K_p V_t` from the state and from the moments; one tissue after a
//!   bolus is the bi-exponential with `k21 = Q/(K_p V_t)`; as `Q → ∞` the model is one
//!   compartment and the deviation falls as `1/Q`; and the well-stirred liver passes `1 − E`
//!   of a dose that enters through it.
//! - **Mass balance.** Dose = still in the compartments + cleared + still in a depot + lost
//!   before absorption + still in the syringe or the bottle, to machine precision.
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
//! the total. A perfusion is the same: a `K_p` on the blood's side instead of the tissue's, or a
//! perfusion with no return flow, closes the books to the last bit and is wrong everywhere else.
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
//! # Oral doses, and why a depot is not a compartment
//!
//! First-order absorption is `dA_d/dt = −kₐA_d`, a rate constant acting on an *amount*, and
//! everything above is a clearance acting on a *concentration*. They are different mechanisms,
//! not one mechanism with a different number: a depot faked as a compartment with a large volume
//! and a link would be a fiction with a volume nobody can name, and the closed forms would not
//! check it. So an absorption route is its own thing, a [`Depot`], declared with
//! [`CompartmentModel::depot`] — a rate constant `kₐ`, a bioavailability `F` and a lag `t_lag` —
//! and dosed with [`CompartmentModel::oral`]. Any first-order extravascular route is the same
//! shape: an intramuscular or subcutaneous injection is a depot with its own `kₐ` and `F`.
//!
//! **What happens to a dose.** At its scheduled time the bioavailable part `F·D` enters the
//! depot and the rest, `(1 − F)·D`, is booked at once as
//! [lost pre-systemically](CompartmentModel::presystemic_loss). From `t + t_lag` the depot empties
//! into its compartment at `kₐ·A_d`. Booking the loss at the moment of dosing rather than as the
//! dose is absorbed is the convention every closed form here is written in — they see only `F·D`
//! — and it makes the loss a fixed number per dose rather than a time course nobody measures.
//!
//! **Integrated exactly, not by Euler.** The depot is a linear pool fed by nothing between doses,
//! so over any interval `s` it keeps exactly `e^{−kₐs}` of what it held. Each step multiplies by
//! that factor, evaluated with the kernel's [`exp`](pantometry_core::math::exp) so the bits are
//! the same on every platform, and hands the difference to the compartment. The depot therefore
//! has no truncation error at all — only `exp`'s half-ulp and one product per step — and adds
//! nothing to [`Domain::max_stable_dt`]: a `kₐ` of a hundred per hour does not shorten the step.
//! What *is* first order is the compartment the absorbed drug arrives in, which is the Euler step
//! it always was; the absorbed amount enters it the way an infusion's slice does, beside the
//! snapshot every flux is taken from.
//!
//! **On the books.** A dose not yet taken is [pending](CompartmentModel::pending_dose), exactly
//! as an infusion's syringe is. Drug in the depot is **held**, for the syringe's reason: it is in
//! the model and not yet in the body, and without it a correct model would report the dose
//! vanishing at the moment it is swallowed and reappearing as it is absorbed. It is not part of
//! the [body burden](CompartmentModel::body_burden), because the gut lumen is outside the
//! systemic circulation and a blood sample cannot see it. The pre-systemic loss is a running
//! total beside `eliminated`, for `eliminated`'s reason — it has left the model, there is nobody
//! on the bus to hand it to, and so nothing is published. The sum is the dose:
//!
//! ```text
//! bolused + infused + oral = in compartments + eliminated + in depots + lost pre-systemically + pending
//! ```
//!
//! The ledger is blind to *where* a depot delivers, for the reason it is blind to links: the
//! amount leaves one entry and joins another in the same sum. A depot feeding the wrong
//! compartment, or `F` applied twice, closes the books to the last bit. What sees them is the
//! closed forms below, and that is what they are for.
//!
//! **The closed forms for oral input.** One compartment gives the Bateman function,
//!
//! ```text
//! C(t) = F·D·kₐ / (V (kₐ − k)) · (e^{−kt} − e^{−kₐt}),   k = CL/V
//! ```
//!
//! with its peak at `t_max = ln(kₐ/k)/(kₐ − k)`, its singular limit `F·D·k·t·e^{−kt}/V` at
//! `kₐ = k`, and an area `AUC₀^∞ = F·D/CL` that does not depend on `kₐ` at all. [`bateman`] is
//! that function, written so the singular limit is a value rather than a division by zero. Two
//! compartments with oral input into the central one give a tri-exponential in `α`, `β` and `kₐ`,
//! derived in the test that checks it.
//!
//! # Perfusion-limited distribution, and the one blood pool it assumes
//!
//! A physiologically based model names its compartments after organs and fixes their volumes
//! and blood flows from anatomy, where a classical model fits them to a curve. The organ is
//! **perfusion-limited** when its membranes are fast compared with its blood flow, so the blood
//! leaving it has equilibrated with it: the tissue holds `K_p` times the concentration of that
//! venous blood, and [`CompartmentModel::perfuse`] carries
//!
//! ```text
//! V_t dC_t/dt = Q (C_b − C_t/K_p) − CL_t C_t/K_p
//! ```
//!
//! from a blood pool into a tissue. That is not the symmetric link — `C_t/K_p` and not `C_t` —
//! and it is why `perfuse` is its own method rather than a link with a coefficient: the order of
//! its arguments is the physics. At `K_p = 1` it is the symmetric link, to the bit.
//!
//! **One well-stirred blood pool, tissues in parallel off it.** Arterial blood is the pool's
//! concentration and every tissue's venous return mixes back into it. The alternative —
//! separate arterial and venous pools, a lung between them, the gut and spleen draining through
//! the portal vein into the liver — is a flow network, and its circulation is one-way: drug
//! moves *with* the blood around a loop. That is a different kind of link from this one, which
//! is an exchange with a single pool, and it would need its own closed forms. What the choice
//! implies, so nobody reads more into a result than it holds:
//!
//! - **No pulmonary first pass and no arterio-venous difference in the blood itself.** A blood
//!   sample anywhere reads `C_b`. The difference across an *organ* is kept — it is the
//!   `C_b − C_t/K_p` the flux is made of, and the liver's extraction is exactly it.
//! - **No portal circulation.** Drug absorbed from the gut reaches the liver before the blood
//!   only if it is put there: a [`Depot`] into the liver tissue is how hepatic first pass is
//!   written, and it passes `1 − E` of the dose, checked. A depot into the blood skips the liver.
//! - **No tissue is in series with another**, so `perfuse` refuses a tissue as a blood pool and a
//!   blood pool as a tissue. Each tissue has one blood flow, and is perfused once.
//!
//! **Elimination from a tissue is on its venous concentration**: an [`eliminating`] compartment
//! that is then perfused clears `CL·C_t/K_p`. That is the well-stirred liver, with an extraction
//! `E = CL/(Q + CL)` and an organ clearance `Q·E`. **Protein binding is not modelled**: `K_p` is
//! a total partition and `CL` acts on total venous concentration, so a caller with an unbound
//! fraction passes `f_u·CL_int` and has the textbook `f_u·CL_int/(Q + f_u·CL_int)`. Modelling
//! `f_u` would mean an unbound concentration that the membranes, the enzymes and `K_p` all
//! refer to, and a check of each; a multiplier on `CL` that nothing else reads would not be that.
//!
//! **No reference human.** A physiological model is only as good as its organ volumes and blood
//! flows, and those come from a published table — ICRP Publication 89, or Brown et al. 1997 —
//! read by somebody who can cite the page. This crate has the mechanism and no such table; the
//! tests use round numbers that say so.
//!
//! **Partition coefficients from composition** are the one place it does carry data, because the
//! data is printed beside the equation that uses it. [`partition`] is Rodgers and Rowland's `K_p`
//! for a neutral molecule, with the rat tissue compositions their papers tabulate, each constant
//! cited to its table and page.
//!
//! [`eliminating`]: CompartmentModel::eliminating
//!
//! # What is deliberately not in it
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
//! compartments a caller can build out of what is here — binding as an equilibrium inside a
//! perfused tissue is folded into its `K_p`, and in the liver into `f_u·CL_int`, as above, and
//! is not a separate mechanism. The third is `Σ eliminated / CL` and
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

pub mod partition;

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
    /// m³/s out of the *system*. Zero for a compartment that only distributes. It acts on the
    /// concentration *leaving* the compartment, `amount/(partition·volume)`, which is the
    /// compartment's own concentration for everything that is not a perfused tissue.
    clearance: f64,
    /// `K_p`, the tissue:blood partition coefficient, for a tissue declared by
    /// [`CompartmentModel::perfuse`]; exactly 1 for every other compartment, so that
    /// `concentration / partition` is the concentration itself, to the bit.
    partition: f64,
}

/// An intercompartmental clearance, in m³/s.
struct Link {
    a: usize,
    b: usize,
    q: f64,
}

/// A blood flow through a perfusion-limited tissue, in m³/s. The tissue's `K_p` is on the
/// tissue's [`CompartmentState`], because a tissue is perfused once and its partition is a
/// property of the tissue that elimination reads too.
struct Perfusion {
    blood: usize,
    tissue: usize,
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

/// An absorption depot in one particular model: the route an oral, or any first-order
/// extravascular, dose takes into a compartment.
///
/// A handle for the reason [`Compartment`] is one — it can only come from
/// [`CompartmentModel::depot`], so a dose into a depot that does not exist is not representable,
/// and it carries the model's identity, so one from another model is refused rather than
/// silently dosing whatever sits at that position. The same limit applies: two models with one
/// name share an identity.
///
/// **The readers are stricter than [`Compartment`]'s.** [`CompartmentModel::amount`] and its
/// siblings index by position and do not check the model, so a handle from another model reads
/// whatever compartment sits there. Every depot reader — [`CompartmentModel::depot_amount`],
/// [`absorption_rate`](CompartmentModel::absorption_rate) and the rest — instead **panics** on a
/// handle from another model or a later state, and [`CompartmentModel::oral`] refuses it with a
/// [`Violation`]. Stricter because it was new and could be; changing the compartment readers
/// to match is a change to an existing API and is not made here.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Depot {
    index: u32,
    model: u64,
}

/// A first-order absorption route: where it delivers, how fast, how much of a dose, how late.
struct DepotState {
    label: String,
    /// The compartment it empties into.
    into: usize,
    /// kₐ, per second.
    rate: f64,
    /// F, in (0, 1].
    bioavailability: f64,
    /// Seconds from a dose to the start of its absorption.
    lag: f64,
}

/// One dose given by a depot's route at a time on the simulation clock.
///
/// Each dose is its own record rather than all of a depot's doses pooled, because each has its
/// own clock: two doses an hour apart start absorbing an hour apart, and a pool with one lag
/// could not say so. The model is linear, so the records superpose.
struct OralDose {
    depot: usize,
    /// Seconds on the simulation clock.
    at: f64,
    /// kg, the whole dose.
    dose: f64,
    /// kg, `F·dose`. Held rather than recomputed so the bioavailable and the lost parts are
    /// computed once each and always the same two numbers.
    bioavailable: f64,
    /// Whether the clock has reached `at` and the dose has been taken.
    given: bool,
    /// kg still in the depot from this dose. Zero until it is given.
    remaining: f64,
}

/// `(1 − e^{−x})/x` for `x ≥ 0`, with its limit of 1 at `x = 0` and no cancellation near it.
///
/// For `x ≥ ½` the subtraction loses nothing that matters: `e^{−x} ≤ 0.61`, so `1 − e^{−x}` is
/// at least 0.39 and its relative error is under three ulps. Below that it is the Taylor series
/// `Σ (−x)ⁿ/(n+1)!`, nested, through `n = 18`: at `x = ½` the first term left out is
/// `0.5¹⁹/20! ≈ 8e-25`, far under an ulp of a value near 1. Nothing is divided by `x` on that
/// branch, which is the point.
fn one_minus_exp_over(x: f64) -> f64 {
    if x < 0.5 {
        let mut sum = 1.0;
        for n in (2..=19).rev() {
            sum = 1.0 - x / n as f64 * sum;
        }
        sum
    } else {
        (1.0 - pantometry_core::math::exp(-x)) / x
    }
}

/// The Bateman function: the concentration a single first-order oral dose produces in one
/// compartment,
///
/// ```text
/// C(t) = F·D·kₐ / (V (kₐ − k)) · (e^{−kt} − e^{−kₐt})
/// ```
///
/// with `F` the bioavailability, `kₐ` the absorption rate constant, `k = CL/V` the elimination
/// rate constant and `t` the time since absorption began (subtract a lag before calling). Zero
/// for `t ≤ 0`.
///
/// # Written so the singular limit is a number
///
/// The textbook spelling divides by `kₐ − k`, which is zero when the two rates are equal and
/// catastrophically cancelling when they are close: at `kₐ = k(1 + 10⁻⁸)` the difference of
/// exponentials keeps about eight good digits. This is the same function rearranged with no
/// division by the difference,
///
/// ```text
/// (e^{−kt} − e^{−kₐt}) / (kₐ − k) = t · e^{−s t} · φ(|kₐ − k| t),   s = min(k, kₐ),
/// φ(x) = (1 − e^{−x}) / x,   φ(0) = 1,
/// ```
///
/// with `φ` from its Taylor series near zero. It is symmetric in `k` and `kₐ`, as the function
/// is, and at `kₐ = k` it is `F·D·k·t·e^{−kt}/V` exactly — the limit, rather than a branch
/// written to imitate it. Factoring the *slower* exponential out also keeps a flip-flop curve
/// (`kₐ < k`) from underflowing to zero while its true value is still representable.
///
/// The exponentials are the kernel's, so the bits are the same on every platform.
///
/// # Its domain, which it does not check
///
/// It is the closed form for `kₐ > 0`, `k > 0`, `V > 0` and `F` in `(0, 1]` — the inputs
/// [`CompartmentModel::depot`] accepts — and it validates none of them, because it is a formula
/// rather than a model and has no `Result` to refuse in. Outside that domain it returns whatever
/// the algebra gives: a negative rate grows instead of decaying, and a zero volume divides by
/// zero.
pub fn bateman(
    dose: Mass,
    bioavailability: f64,
    absorption: Frequency,
    elimination: Frequency,
    volume: Volume,
    t: Time,
) -> Concentration {
    let (ka, k, t) = (absorption.to_si(), elimination.to_si(), t.to_si());
    if t <= 0.0 {
        return Concentration::from_si(0.0);
    }
    let slow = ka.min(k);
    let shape = t * pantometry_core::math::exp(-slow * t) * one_minus_exp_over((ka - k).abs() * t);
    Concentration::from_si(bioavailability * dose.to_si() * ka / volume.to_si() * shape)
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
    perfusions: Vec<Perfusion>,
    infusions: Vec<Infusion>,
    depots: Vec<DepotState>,
    doses: Vec<OralDose>,
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
/// different clothes: restore without it and the syringe refills. An oral dose's `given` and
/// `remaining` are the same hazard again — restore the compartments without them and a dose is
/// absorbed twice.
#[derive(Clone)]
struct Saved {
    amounts: Vec<f64>,
    delivered: Vec<f64>,
    doses: Vec<(bool, f64)>,
    eliminated: f64,
    bolused: f64,
}

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
            perfusions: Vec::new(),
            infusions: Vec::new(),
            depots: Vec::new(),
            doses: Vec::new(),
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
    ///
    /// In a perfusion-limited model the liver is this too, and then perfused: the clearance acts
    /// on the concentration *leaving* the compartment, which for a tissue declared by
    /// [`CompartmentModel::perfuse`] is its venous concentration `C_t/K_p` rather than `C_t`.
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
            partition: 1.0,
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
    ///
    /// A tissue whose equilibrium with the blood is not equal concentrations is
    /// [`CompartmentModel::perfuse`], which carries `Q·(C_b − C_t/K_p)`; this is that with
    /// `K_p = 1`, to the bit.
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

    /// Perfuse a tissue from a blood pool: a blood flow `Q` through the tissue, which holds drug
    /// at `K_p` times the blood's concentration once the two have equilibrated.
    ///
    /// Carries `Q·(C_b − C_t/K_p)` kilograms a second from `blood` into `tissue`. Blood arrives
    /// at the blood pool's concentration `C_b`, and leaves the tissue at its **venous**
    /// concentration `C_t/K_p` — the perfusion-limited, well-stirred tissue, in which the blood
    /// leaving is in equilibrium with the tissue it has just passed through. The flux vanishes
    /// at `C_t = K_p·C_b`, which is what `K_p` means. The micro rate constants are `Q/V_b` out
    /// of the blood and `Q/(K_p·V_t)` out of the tissue, and that second one is the one
    /// [`Domain::max_stable_dt`] must see: a small tissue that partitions weakly empties fast.
    ///
    /// **Not symmetric**, which is the whole difference from [`CompartmentModel::link`]: the
    /// order of the arguments is the physics. `K_p` belongs to the second, and swapping them
    /// is a different model that conserves every milligram.
    ///
    /// **Elimination from a perfused tissue acts on its venous concentration.** A tissue made
    /// with [`CompartmentModel::eliminating`] and then perfused clears `CL·C_t/K_p`: the
    /// well-stirred liver, in which the enzymes see the blood leaving the organ. At steady state
    /// that extracts `E = CL/(Q + CL)` of the drug arriving, and the organ clears
    /// `Q·E = Q·CL/(Q + CL)` of blood — see [`CompartmentModel::extraction_ratio`]. `CL` is the
    /// intrinsic clearance of what the enzymes can reach; protein binding is not modelled, so a
    /// caller with an unbound fraction `f_u` and an unbound intrinsic clearance passes
    /// `f_u·CL_int` and gets the textbook `E = f_u·CL_int/(Q + f_u·CL_int)`.
    ///
    /// **Tissues in parallel off one blood pool.** A tissue is perfused once, and a blood pool is
    /// not itself a perfused tissue: so a tissue cannot also be another's blood, and nothing is
    /// in series. That is the lumped-blood model — see the module docs for what it implies.
    /// Symmetric links to and from a perfused tissue are allowed and act on its concentration
    /// `C_t`, as a link always has.
    ///
    /// `K_p = 1` is the symmetric link of the same `Q` **to the bit**: `x/1.0` is `x` exactly,
    /// and the flux, the elimination, the step limit and the rate constants are then computed in
    /// the same operations in the same order.
    ///
    /// Refuses a blood flow that is negative or not finite, a `K_p` that is not positive and
    /// finite, a tissue perfused from itself, a tissue already perfused, a tissue that is some
    /// other tissue's blood, a blood pool that is itself a perfused tissue, and a handle from
    /// another model. Nothing is changed by a refusal.
    pub fn perfuse(
        &mut self,
        blood: Compartment,
        tissue: Compartment,
        q: Clearance,
        partition: f64,
    ) -> Result<(), Violation> {
        let (b, t) = (self.resolve(blood)?, self.resolve(tissue)?);
        let site = format!(
            "{}/{}-{}",
            self.name, self.compartments[b].label, self.compartments[t].label
        );
        if b == t {
            return Err(Violation::at(
                site,
                "a tissue cannot be perfused from itself",
                0.0,
            ));
        }
        let w = q.to_si();
        if !w.is_finite() || w < 0.0 {
            return Err(Violation::at(
                site,
                "a blood flow must be finite and not negative",
                w,
            ));
        }
        // `!is_finite` first, for the reason `link` gives: NaN is refused by the branch that
        // reads as refusing it.
        if !partition.is_finite() || partition <= 0.0 {
            return Err(Violation::at(
                site,
                "a partition coefficient must be finite and positive",
                partition,
            ));
        }
        if self.perfusions.iter().any(|p| p.tissue == t) {
            return Err(Violation::at(
                site,
                "a tissue already perfused cannot be perfused again",
                t as f64,
            ));
        }
        if self.perfusions.iter().any(|p| p.blood == t) {
            return Err(Violation::at(
                site,
                "a blood pool cannot be perfused as a tissue",
                t as f64,
            ));
        }
        if self.perfusions.iter().any(|p| p.tissue == b) {
            return Err(Violation::at(
                site,
                "a perfused tissue cannot be a blood pool",
                b as f64,
            ));
        }
        self.compartments[t].partition = partition;
        self.perfusions.push(Perfusion {
            blood: b,
            tissue: t,
            q: w,
        });
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

    /// A first-order absorption route into a compartment: the gut for an oral dose, or the
    /// injection site for an intramuscular or subcutaneous one.
    ///
    /// `absorption` is `kₐ`, a rate constant on the *amount* in the depot — not a clearance,
    /// which is why it is a [`Frequency`] where everything else here is a [`Clearance`].
    /// `bioavailability` is `F`, the fraction of each dose that reaches the compartment at all;
    /// the rest is lost before absorption, to first-pass metabolism or to the gut. `lag` is the
    /// time from taking a dose to the start of its absorption, the stomach's emptying time for a
    /// tablet. Dose it with [`CompartmentModel::oral`].
    ///
    /// Refuses a `kₐ` that is not positive and finite — zero would be a depot that never empties,
    /// whose dose the books would hold forever while every curve stayed flat — an `F` outside
    /// `(0, 1]` or not a number, a negative or non-finite lag, and a compartment from another
    /// model. Several depots may feed one compartment, which is how a regimen of two
    /// formulations is written.
    pub fn depot(
        &mut self,
        label: impl Into<String>,
        into: Compartment,
        absorption: Frequency,
        bioavailability: f64,
        lag: Time,
    ) -> Result<Depot, Violation> {
        let label = label.into();
        let i = self.resolve(into)?;
        let site = format!("{}/{}", self.name, label);
        let ka = absorption.to_si();
        if !ka.is_finite() || ka <= 0.0 {
            return Err(Violation::at(
                site,
                "an absorption rate constant must be finite and positive",
                ka,
            ));
        }
        // NaN named first, so it is refused by the branch that reads as refusing it rather than
        // by two comparisons that both happen to be false.
        let f = bioavailability;
        if f.is_nan() || f <= 0.0 || f > 1.0 {
            return Err(Violation::at(
                site,
                "a bioavailability must be in (0, 1]",
                f,
            ));
        }
        let lag = lag.to_si();
        if !lag.is_finite() || lag < 0.0 {
            return Err(Violation::at(
                site,
                "an absorption lag must be finite and not negative",
                lag,
            ));
        }
        let index = self.depots.len() as u32;
        self.depots.push(DepotState {
            label,
            into: i,
            rate: ka,
            bioavailability: f,
            lag,
        });
        Ok(Depot {
            index,
            model: self.id,
        })
    }

    /// An oral dose — or any dose by a depot's route — taken at `at` on the simulation clock.
    ///
    /// Scheduled rather than given now, the way [`CompartmentModel::infuse`] is, because a lag
    /// has to be measured from somewhere and the model's only clock is the one `step` is handed.
    /// From the moment it is scheduled the whole dose is on the books as
    /// [`pending_dose`](CompartmentModel::pending_dose). At the first step whose interval reaches
    /// `at` it is taken: `F·dose` enters the depot, `(1 − F)·dose` is booked as
    /// [lost](CompartmentModel::presystemic_loss), and absorption begins at `at + lag`, exactly
    /// there even part-way through a step. A dose scheduled before the clock a model is next
    /// stepped from is taken at that step and starts absorbing from its start; nothing that
    /// would have been absorbed earlier arrives all at once.
    ///
    /// Any number of doses may be scheduled, at any times and in any order. The model is linear,
    /// so they superpose: two doses give the sum of the two curves each would give alone.
    ///
    /// Refuses a negative or non-finite dose, a non-finite time, and a depot from another model.
    /// A zero dose is allowed and does nothing, as for [`CompartmentModel::bolus`].
    pub fn oral(&mut self, depot: Depot, dose: Mass, at: Time) -> Result<(), Violation> {
        let d = self.resolve_depot(depot)?;
        let site = format!("{}/{}", self.name, self.depots[d].label);
        let m = dose.to_si();
        if !m.is_finite() || m < 0.0 {
            return Err(Violation::at(
                site,
                "a dose must be finite and not negative",
                m,
            ));
        }
        let at = at.to_si();
        if !at.is_finite() {
            return Err(Violation::at(
                site,
                "a dose must be taken at a finite time",
                at,
            ));
        }
        self.doses.push(OralDose {
            depot: d,
            at,
            dose: m,
            bioavailable: self.depots[d].bioavailability * m,
            given: false,
            remaining: 0.0,
        });
        Ok(())
    }

    fn resolve_depot(&self, d: Depot) -> Result<usize, Violation> {
        if d.model != self.id {
            return Err(Violation::at(
                self.name.clone(),
                "a depot handle from a different model",
                d.index as f64,
            ));
        }
        let i = d.index as usize;
        if i >= self.depots.len() {
            return Err(Violation::at(
                self.name.clone(),
                "a depot handle from a later state of this model",
                d.index as f64,
            ));
        }
        Ok(i)
    }

    /// A depot by label, for a caller holding names off a file — the counterpart of
    /// [`CompartmentModel::compartment_named`].
    pub fn depot_named(&self, label: &str) -> Option<Depot> {
        self.depots
            .iter()
            .position(|d| d.label == label)
            .map(|i| Depot {
                index: i as u32,
                model: self.id,
            })
    }

    /// Drug in one depot right now, from every dose given into it: absorbed by nobody yet.
    ///
    /// Like every depot accessor, this **panics** on a handle from another model, rather than
    /// returning zero or reading whichever depot sits at that position here — see
    /// [`Depot`].
    pub fn depot_amount(&self, d: Depot) -> Mass {
        let i = self.depot_index(d);
        Mass::from_si(
            self.doses
                .iter()
                .filter(|x| x.depot == i)
                .map(|x| x.remaining)
                .fold(0.0, |a, x| a + x),
        )
    }

    /// Drug in every depot right now. Held by the model and not in the body — see the module
    /// docs for why it is counted and why it is not part of [`CompartmentModel::body_burden`].
    pub fn depot_mass(&self) -> Mass {
        Mass::from_si(
            self.doses
                .iter()
                .map(|x| x.remaining)
                .fold(0.0, |a, x| a + x),
        )
    }

    /// Drug lost before absorption over the run: `(1 − F)·dose` for every dose taken so far.
    ///
    /// Summed from the doses rather than tallied as it goes, for the reason
    /// [`administered_mass`](CompartmentModel::administered_mass) is: then it is always the same
    /// number each dose's own subtraction gave.
    pub fn presystemic_loss(&self) -> Mass {
        Mass::from_si(
            self.doses
                .iter()
                .filter(|x| x.given)
                .map(|x| x.dose - x.bioavailable)
                .fold(0.0, |a, x| a + x),
        )
    }

    /// A depot's absorption rate constant, `kₐ`.
    pub fn absorption_rate(&self, d: Depot) -> Frequency {
        Frequency::from_si(self.depots[self.depot_index(d)].rate)
    }

    /// A depot's bioavailability, `F`.
    pub fn bioavailability(&self, d: Depot) -> f64 {
        self.depots[self.depot_index(d)].bioavailability
    }

    /// A depot's absorption lag.
    pub fn absorption_lag(&self, d: Depot) -> Time {
        Time::from_si(self.depots[self.depot_index(d)].lag)
    }

    /// The compartment a depot empties into.
    pub fn absorbed_into(&self, d: Depot) -> Compartment {
        Compartment {
            index: self.depots[self.depot_index(d)].into as u32,
            model: self.id,
        }
    }

    /// A depot's position, or a panic naming why there is none.
    ///
    /// The readers panic where the methods that change the model return a [`Violation`], and
    /// the difference is deliberate: a reader has no `Result` to put a refusal in, and the two
    /// other answers it could give — zero, or another model's depot at the same position — are
    /// both a number that looks like a measurement.
    fn depot_index(&self, d: Depot) -> usize {
        match self.resolve_depot(d) {
            Ok(i) => i,
            Err(e) => panic!("{e}"),
        }
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
    ///
    /// For a perfused tissue this is the *intrinsic* clearance, acting on the venous
    /// concentration; what the organ clears of blood is `Q·E`, with `E` from
    /// [`CompartmentModel::extraction_ratio`].
    pub fn clearance(&self, c: Compartment) -> Clearance {
        Clearance::from_si(self.compartments[c.index as usize].clearance)
    }

    /// A compartment's tissue:blood partition coefficient `K_p`. Exactly 1 for a compartment
    /// that is not a perfused tissue.
    pub fn partition(&self, c: Compartment) -> f64 {
        self.compartments[c.index as usize].partition
    }

    /// The concentration *leaving* a compartment: `C_t/K_p` for a perfused tissue — the venous
    /// blood in equilibrium with it — and the compartment's own concentration, to the bit, for
    /// anything else. The concentration its clearance acts on.
    pub fn venous_concentration(&self, c: Compartment) -> Concentration {
        let i = c.index as usize;
        Concentration::from_si(self.conc(i) / self.compartments[i].partition)
    }

    /// The blood flow through a perfused tissue. Zero for a compartment that is not perfused,
    /// which is a real answer: no blood flows through it.
    pub fn blood_flow(&self, tissue: Compartment) -> Clearance {
        let t = tissue.index as usize;
        Clearance::from_si(
            self.perfusions
                .iter()
                .find(|p| p.tissue == t)
                .map(|p| p.q)
                .unwrap_or(0.0),
        )
    }

    /// The well-stirred extraction ratio of a perfused tissue, `E = CL/(Q + CL)`: the fraction
    /// of the drug arriving in its blood that it removes in one pass, at steady state.
    ///
    /// `None` for a compartment that is not perfused, rather than zero, because zero is the
    /// answer for a perfused tissue that does not eliminate and the two must not be confused.
    /// Zero, too, for one with neither flow nor clearance, where `0/0` has no better value.
    pub fn extraction_ratio(&self, tissue: Compartment) -> Option<f64> {
        let t = tissue.index as usize;
        let q = self.perfusions.iter().find(|p| p.tissue == t)?.q;
        let cl = self.compartments[t].clearance;
        Some(if q + cl > 0.0 { cl / (q + cl) } else { 0.0 })
    }

    /// The perfusion between two compartments, either way round, if there is one.
    fn perfusion_between(&self, i: usize, j: usize) -> Option<&Perfusion> {
        self.perfusions
            .iter()
            .find(|p| (p.blood == i && p.tissue == j) || (p.blood == j && p.tissue == i))
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
    ///
    /// Not including what is still in a depot: drug in the gut lumen has not reached the
    /// systemic circulation, and that is what a body burden is the burden of.
    pub fn body_burden(&self) -> Mass {
        Mass::from_si(self.compartments.iter().map(|c| c.amount).sum())
    }

    /// Drug removed from the body by clearance over the run.
    ///
    /// Still on this domain's books — see the module docs. `body_burden + eliminated_mass +
    /// depot_mass + presystemic_loss + pending_dose` is the dose, exactly.
    pub fn eliminated_mass(&self) -> Mass {
        Mass::from_si(self.eliminated)
    }

    /// Drug given over the run, by bolus, by delivered infusion, and by oral doses taken.
    ///
    /// Summed from the infusions and the doses rather than tallied as it goes, so it agrees
    /// exactly with `scheduled − pending` instead of to within an accumulation of roundings. An
    /// oral dose counts whole, the part that will be lost before absorption included: it was
    /// administered, and bioavailability is a property of the route rather than of the dose.
    pub fn administered_mass(&self) -> Mass {
        Mass::from_si(
            self.bolused
                + self.infusions.iter().map(|i| i.delivered).sum::<f64>()
                + self
                    .doses
                    .iter()
                    .filter(|d| d.given)
                    .map(|d| d.dose)
                    .sum::<f64>(),
        )
    }

    /// Scheduled drug not given yet — what is still in the syringe, and the oral doses whose time
    /// has not come.
    pub fn pending_dose(&self) -> Mass {
        Mass::from_si(self.pending_infusion() + self.pending_oral())
    }

    fn pending_infusion(&self) -> f64 {
        self.infusions
            .iter()
            .map(|i| (i.total - i.delivered).max(0.0))
            .sum()
    }

    fn pending_oral(&self) -> f64 {
        // From `+0.0` and not by `sum`, which starts at `-0.0`: a filter that keeps nothing
        // printed `still to take -0.000000000e0` in the CSV of the first scene with a tablet.
        self.doses
            .iter()
            .filter(|d| !d.given)
            .map(|d| d.dose)
            .fold(0.0, |a, x| a + x)
    }

    /// Steady-state volume of distribution, `Σ K_p,i V_i` — `Σ V_i` when nothing is perfused,
    /// and `V_b + Σ K_p V_t` for blood and its tissues.
    ///
    /// The volume that relates the body burden to the plasma concentration once everything has
    /// equilibrated, which is the only time it relates them at all. Early after a bolus the
    /// central volume is the relevant one and this is an overestimate by however much has not
    /// distributed yet.
    ///
    /// A tissue holds `K_p` times the blood's concentration at equilibrium, so it counts `K_p`
    /// times its volume; `1.0·V` is `V` to the bit, so a model with no perfusion reports what it
    /// always did. **Exact when only blood pools eliminate.** A perfused tissue that eliminates
    /// sits below equilibrium at steady state — at `K_p·C_b·(1 − E)` — and is counted here as
    /// if it did not. That is the convention of the prediction formulas `K_p` comes from, and
    /// **not** what plasma data would give: with the test body's liver clearing 60 L/h this
    /// reports 122.2 L, the steady state's `body burden/C_b` is 119.3 L, and `CL·MRT` from a
    /// bolus 117.6 L. Ask for the one that is meant.
    pub fn volume_of_distribution(&self) -> Volume {
        Volume::from_si(
            self.compartments
                .iter()
                .map(|c| c.partition * c.volume)
                .sum(),
        )
    }

    /// Total elimination clearance, `Σ CL_i`.
    ///
    /// At steady state under a constant infusion `Σ CL_i C_i = R`, whatever the distribution
    /// between compartments — so with elimination from one compartment only, its steady-state
    /// concentration is `R/CL` and does not depend on the intercompartmental clearances at all.
    ///
    /// **A sum of what each compartment's enzymes or kidneys do, not of what the blood sees.**
    /// For a perfused tissue `CL_i` is intrinsic and acts on `C_t/K_p`; what it clears of the
    /// blood is `Q·E = Q·CL/(Q + CL)`, from [`CompartmentModel::blood_flow`] and
    /// [`CompartmentModel::extraction_ratio`], and with infusion into the blood the steady state
    /// is `C_b = R/(Σ CL_blood + Σ Q·E)`. The same is true of a linked peripheral compartment that
    /// eliminates, which is why the sentence above says one compartment; this sum is unchanged so
    /// that a perfusion of `K_p = 1` and a link report the same number.
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
    ///
    /// For a perfused tissue it is `CL/(K_p V)`, because the clearance acts on `C_t/K_p`; with
    /// `K_p = 1` that is `CL/V` to the bit.
    pub fn elimination_rate(&self, c: Compartment) -> Frequency {
        let s = &self.compartments[c.index as usize];
        Frequency::from_si(if s.volume > 0.0 {
            s.clearance / (s.partition * s.volume)
        } else {
            0.0
        })
    }

    /// The micro rate constant for drug leaving `from` towards `to`: `Q/V_from` for a link, and
    /// for a perfusion `Q/V_b` out of the blood and `Q/(K_p V_t)` out of the tissue.
    ///
    /// Zero if they are not joined, which is a real answer rather than a missing one. Not
    /// symmetric: the same `Q` gives a different constant in each direction whenever the volumes
    /// differ, and that asymmetry is what makes the bi-exponential bi-exponential. A link and a
    /// perfusion between the same pair add.
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
        let s = &self.compartments[i];
        let v = s.volume;
        let linked = if v > 0.0 { q / v } else { 0.0 };
        // Only added when there is a perfusion, so a model without one computes exactly what it
        // did before there were perfusions.
        Frequency::from_si(match self.perfusion_between(i, j) {
            None => linked,
            Some(p) if v > 0.0 => {
                // Out of the tissue it is the venous concentration that leaves, `A/(K_p V)`.
                let k = if p.tissue == i { s.partition } else { 1.0 };
                linked + p.q / (k * v)
            }
            Some(_) => linked,
        })
    }

    /// Mass crossing the link or perfusion between two compartments right now, positive from `a`
    /// to `b`: `Q·(C_a − C_b)` for a link, and `Q·(C_b − C_t/K_p)` from blood to tissue for a
    /// perfusion.
    ///
    /// Zero if they are not joined. The per-link number the conservation audit cannot see, which
    /// is why it is public: a test asserting it is the only kind that would notice a transposed
    /// index — or a `K_p` on the wrong side.
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
        let linked = q * (self.conc(i) - self.conc(j));
        MassFlow::from_si(match self.perfusion_between(i, j) {
            None => linked,
            Some(p) => {
                let into_tissue = p.q
                    * (self.conc(p.blood)
                        - self.conc(p.tissue) / self.compartments[p.tissue].partition);
                linked
                    + if p.blood == i {
                        into_tissue
                    } else {
                        -into_tissue
                    }
            }
        })
    }

    fn conc(&self, i: usize) -> f64 {
        let s = &self.compartments[i];
        if s.volume > 0.0 {
            s.amount / s.volume
        } else {
            0.0
        }
    }

    /// Everything carrying drug *out* of compartment `i`, as a rate constant on its amount:
    /// `(ΣQ + CL)/V` — and for a perfused tissue `Σ Q_link/V + (Q + CL)/(K_p V)`, because the
    /// perfusion and the clearance both drain it at its venous concentration `A/(K_p V)`. A
    /// blood pool's perfusions drain it at `Q/V_b` each, as a link would.
    ///
    /// The magnitude of the diagonal of the compartment matrix written in amounts, and the only
    /// thing [`Domain::max_stable_dt`] needs. Infinite for a compartment with no volume, which
    /// `step` refuses before this is ever divided by.
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
            // Drained at the compartment's own concentration, and at its venous one. Each sum is
            // empty, so exactly zero, for a compartment no perfusion touches — and `x + 0` is
            // `x`, so a model with no perfusion takes the same operations it always did.
            let as_blood: f64 = self
                .perfusions
                .iter()
                .filter(|p| p.blood == i)
                .map(|p| p.q)
                .sum();
            let as_tissue: f64 = self
                .perfusions
                .iter()
                .filter(|p| p.tissue == i)
                .map(|p| p.q)
                .sum();
            let (own, venous) = (out + as_blood, s.clearance + as_tissue);
            if s.partition == 1.0 {
                (own + venous) / s.volume
            } else {
                own / s.volume + venous / (s.partition * s.volume)
            }
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
    /// only into `eliminated` or the pre-systemic loss, both of which the ledger reports. There is no bus channel, so there is
    /// nothing for the difference to be attributed to and the audit is comparing the books
    /// against themselves one step apart. That is the strictest form the check has and this
    /// domain is the cheapest possible case of it.
    fn books_balance(&self) -> bool {
        true
    }

    /// `min_i 1/r_i` — the reciprocal of the fastest compartment's turnover, with
    /// `r_i = (Σ_j Q_ij + CL_i)/V_i` for a compartment and
    /// `r_t = Σ Q_link/V_t + (Q + CL)/(K_p V_t)` for a perfused tissue.
    ///
    /// # Where that comes from
    ///
    /// Write the system in amounts, `dA/dt = M A`. Every term draining compartment `i` is
    /// proportional to its own amount, so `M_ii = −r_i`; and `M_ji ≥ 0` for `j ≠ i` is the rate
    /// at which `A_i` feeds `j` — `Q/V_i` through a link, `Q/V_b` from a blood pool into a
    /// tissue, `Q/(K_p V_t)` from a tissue back to its blood. Column `i` of `M` therefore sums
    /// to `−CL_i/(K_p,i V_i) ≤ 0`, which is mass conservation said as a matrix.
    ///
    /// Explicit Euler is stable when `|1 + λ h| ≤ 1` for every eigenvalue. Gershgorin, by
    /// columns, puts every eigenvalue in a disc centred on `−r_i` whose radius — the rest of the
    /// column — is at most `r_i`. For `h r_i ≤ 1`, `1 + hz` over that disc lies in the disc of
    /// centre `1 − h r_i` and radius `h r_i`, which is inside the unit circle because the two add
    /// to 1. So `h ≤ 1/max_i r_i` is stable, and it needs no claim that the eigenvalues are
    /// real. They are for a network of symmetric links, which is similar to a symmetric matrix,
    /// and for blood and its tissues, which satisfy detailed balance with weights `V_b` and
    /// `K_p V_t`; a link added beside a perfusion can make a cycle that breaks that, and the
    /// disc argument does not care.
    ///
    /// The concentration form this used to be derived in does **not** survive perfusion: a
    /// tissue's concentration row has diagonal `Q/(K_p V_t)` and off-diagonal `Q/V_t`, so its
    /// disc is wider than its centre whenever `K_p > 1`; and the blood's row, with off-diagonals
    /// `Q_t/(K_p V_b)`, is too whenever `Σ Q_t/K_p > Σ Q_t + CL` — a small `K_p` breaks it from
    /// the other side. In amounts the columns are what
    /// conservation bounds, which is why the derivation moved there. For a network of links the
    /// two forms give the same number.
    ///
    /// # And it is the positivity limit as well, which is the reason not to loosen it
    ///
    /// The same expression falls out of a different requirement. A compartment's update is
    /// `A_i ← A_i(1 − h r_i) + h Σ_{j≠i} M_ij A_j + dose`, and every term after the first is
    /// non-negative, so the amount stays non-negative exactly while `h r_i ≤ 1`. A negative mass
    /// of drug is not a large error, it is a meaningless one — and unlike a divergence it does
    /// not announce itself, because the system can still be stable past this — up to twice
    /// it for one compartment, `2/α` for the classical pair. Not in general: Gershgorin gives
    /// `|λ| ≤ 2 max r`, so a doubled step can diverge, and once a cycle makes the eigenvalues
    /// complex `2/|λ|` is not the criterion at all.
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
    /// Infinite for a model with no clearances, no links and no perfusions: nothing is moving, so nothing limits
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

        // Oral doses: taken when the clock reaches them, absorbed **exactly** — see the module
        // docs. A depot is fed by nothing between doses, so over the part of this step after its
        // lag it keeps `e^{−kₐs}` of what it held, and the difference arrives in the compartment
        // through `delta`, beside the snapshot, as an infusion's slice does.
        //
        // `s` is `h` itself whenever absorption was already under way when the step began,
        // rather than `until − start`: the two differ in the last bit whenever `now` is not
        // round, and then two runs that should be one shifted by a lag would not be the same
        // bits. Only the step in which absorption *starts* takes the partial span.
        //
        // `remaining − kept` rather than `remaining·(1 − e)`: when `e ≥ ½` Sterbenz makes the
        // subtraction exact, so the depot and the compartment exchange one number rather than
        // two roundings of it.
        for d in self.doses.iter_mut() {
            let route = &self.depots[d.depot];
            if !d.given && d.at < until {
                d.given = true;
                d.remaining = d.bioavailable;
            }
            let start = d.at + route.lag;
            if d.given && d.remaining > 0.0 && until > start {
                let span = if start <= now { h } else { until - start };
                let kept = d.remaining * pantometry_core::math::exp(-route.rate * span);
                delta[route.into] += d.remaining - kept;
                d.remaining = kept;
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

        // What leaves each compartment is at its venous concentration, `C/K_p`: the tissue's
        // blood for a perfused tissue, and the concentration itself, to the bit, for everything
        // else — `x/1.0` is `x`.
        let venous: Vec<f64> = self
            .compartments
            .iter()
            .zip(&conc)
            .map(|(c, x)| x / c.partition)
            .collect();

        // Perfusions, once each and applied twice, for the links' reason. Arterial blood in at
        // the blood pool's concentration, venous blood out at the tissue's `C_t/K_p`. A blood
        // pool is never a perfused tissue — `perfuse` refuses it — so its own concentration is
        // its arterial one.
        for p in &self.perfusions {
            let q = p.q * (conc[p.blood] - venous[p.tissue]) * h;
            delta[p.blood] -= q;
            delta[p.tissue] += q;
        }

        // Elimination, on the venous concentration. Accumulated with exactly the value
        // subtracted, so what left the body and what the books gained are the same number rather
        // than two roundings of one.
        let mut cleared = 0.0;
        for (i, c) in self.compartments.iter().enumerate() {
            let out = c.clearance * venous[i] * h;
            delta[i] -= out;
            cleared += out;
        }
        self.eliminated += cleared;

        for (c, d) in self.compartments.iter_mut().zip(&delta) {
            c.amount += d;
        }
        Ok(())
    }

    /// Every compartment's holding, plus what has been cleared, what is still in the syringe or
    /// the bottle, what is still in a depot, and what was lost before absorption.
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
        // One entry per oral dose still absorbing, for the scale's reason above, and the loss as
        // one total, for `eliminated`'s. Without the depot the books would drop by `F·D` the
        // moment a dose is taken and recover it as it is absorbed — a correct model reported as
        // leaking and then as inventing drug, which is what the depot is held for.
        for d in &self.doses {
            ledger.add(quantity::MASS, d.remaining);
        }
        ledger.add(quantity::MASS, self.presystemic_loss().to_si());
        ledger
    }

    /// Every amount, every infusion's progress, every oral dose's, and both running totals.
    ///
    /// All of it, because `ledger` reads all of it. An infusion's `delivered` is the one easiest
    /// to forget and the one whose loss is silent: restore without it and the syringe refills,
    /// the books gain a dose, and the sweep that was meant to be re-run from the same state is
    /// re-run from a richer one.
    ///
    /// **What it does not rewind is the schedule.** An infusion or an oral dose *scheduled* after
    /// the checkpoint is not removed by [`restore`](Domain::restore), and keeps whatever progress
    /// it has made: it was not in the saved state, so there is nothing to put back. Schedule a
    /// regimen before checkpointing it, which is what an iterative coupling does anyway.
    fn checkpoint(&mut self) {
        self.saved = Some(Saved {
            amounts: self.compartments.iter().map(|c| c.amount).collect(),
            delivered: self.infusions.iter().map(|i| i.delivered).collect(),
            doses: self.doses.iter().map(|d| (d.given, d.remaining)).collect(),
            eliminated: self.eliminated,
            bolused: self.bolused,
        });
    }

    fn restore(&mut self) {
        if let Some(saved) = self.saved.clone() {
            for (c, a) in self.compartments.iter_mut().zip(saved.amounts) {
                c.amount = a;
            }
            for (inf, d) in self.infusions.iter_mut().zip(saved.delivered) {
                inf.delivered = d;
            }
            for (dose, (given, remaining)) in self.doses.iter_mut().zip(saved.doses) {
                dose.given = given;
                dose.remaining = remaining;
            }
            self.eliminated = saved.eliminated;
            self.bolused = saved.bolused;
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
                self.pending_infusion() * 1e6,
                "mg",
            ));
        }
        // The oral half of the books, present when there is an oral dose and absent otherwise,
        // the rule the syringe follows. All three, so the table sums to the dose by itself.
        if !self.doses.is_empty() {
            out.push(Reading::new(
                &self.name,
                "still to take",
                self.pending_oral() * 1e6,
                "mg",
            ));
            out.push(Reading::new(
                &self.name,
                "in the depot",
                self.depot_mass().to_si() * 1e6,
                "mg",
            ));
            out.push(Reading::new(
                &self.name,
                "lost before absorption",
                self.presystemic_loss().to_si() * 1e6,
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
            "the central one should have lost the link and the clearance together: {lost:e} \
             against {:e}",
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
    /// announce itself: the scheme is *stable* to `2/α`, 1.74 times this step, so a model run between the two
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

    // ---- Oral doses ------------------------------------------------------------------------

    /// 100 mg, the oral dose every check below gives unless it says otherwise.
    const ORAL: f64 = 1e-4;

    fn per_hour(x: f64) -> Frequency {
        Frequency::from_si(x / HOUR)
    }

    /// The one-compartment model above — 10 L, 5 L/h, `k = 0.5/h` — with a depot into it.
    fn oral_one(ka_per_h: f64, f: f64, lag_s: f64) -> (CompartmentModel, Compartment, Depot) {
        let (mut m, plasma) = one_compartment();
        let gut = m
            .depot("gut", plasma, per_hour(ka_per_h), f, Time::s(lag_s))
            .unwrap();
        (m, plasma, gut)
    }

    /// The Bateman function as every textbook writes it, in SI, for a dose of [`ORAL`].
    ///
    /// The reference the stepped curves are held to, and deliberately *not* [`bateman`], which is
    /// code under test and is checked against this separately where the two rates are far apart
    /// enough for the division to be harmless.
    fn bateman_textbook(f: f64, ka: f64, k: f64, v: f64, t: f64) -> f64 {
        if t <= 0.0 {
            return 0.0;
        }
        f * ORAL * ka / (v * (ka - k)) * ((-k * t).exp() - (-ka * t).exp())
    }

    /// The concentration at every point of the grid `t = n h`, `n = 0..=steps`.
    ///
    /// The clock is `n·h` rather than a running sum, and every step size below is a whole number
    /// of seconds, so every time is an exact integer and two runs meant to agree are handed the
    /// same clock.
    fn trace(m: &mut CompartmentModel, c: Compartment, h: f64, steps: usize) -> Vec<f64> {
        let mut bus = Exchange::new();
        let mut out = vec![m.concentration(c).to_si()];
        for n in 0..steps {
            m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
            out.push(m.concentration(c).to_si());
        }
        out
    }

    /// The order of convergence between successive halvings: `log₂(e_h / e_{h/2})`.
    fn orders(errors: &[f64]) -> Vec<f64> {
        errors
            .windows(2)
            .map(|w| (w[0] / w[1]).ln() / 2.0f64.ln())
            .collect()
    }

    /// **The exact solution of the recurrence the model steps**, one compartment, one dose of
    /// [`ORAL`] at `t = 0`, on the grid of [`trace`]: element `n` is the concentration after `n`
    /// steps.
    ///
    /// # The derivation
    ///
    /// Step `n` runs from `t = nh`. The depot holds `F·D·Eⁿ` at its start, `E = e^{−kₐh}`, and
    /// hands the compartment `F·D·Eⁿ(1 − E)`; the compartment loses `CL·Cₙ·h = kh·Aₙ` from the
    /// snapshot taken before anything arrives. So with `r = 1 − kh`,
    ///
    /// ```text
    /// A₀ = 0,   Aₙ₊₁ = r Aₙ + F·D(1 − E)Eⁿ
    /// Aₙ = F·D(1 − E) Σ_{j<n} Eʲ r^{n−1−j} = F·D(1 − E)(Eⁿ − rⁿ)/(E − r)
    /// ```
    ///
    /// and `F·D(1 − E)·n·r^{n−1}` in the limit `E = r`. The first two terms check the indexing
    /// against the trace by hand: `A₁ = F·D(1 − E)`, `A₂ = F·D(1 − E)(r + E)`.
    ///
    /// # Written so `kₐ = k` does not cancel
    ///
    /// At `kₐ = k`, `E` and `r` are *not* equal — `E − r = (kh)²/2 + …`, 3e-6 at `h = 18 s` — so
    /// `(Eⁿ − rⁿ)/(E − r)` would keep only ten digits. Instead, with `q = E/r`,
    /// `Σ qʲ = expm1(nL)/expm1(L)` and `L = ln q = (k − kₐ)h + Σ_{m≥2}(kh)ᵐ/m`, the series being
    /// `−kh − ln(1 − kh)` without the subtraction. `r^{n−1}` is built by repeated multiplication.
    /// `E` is the kernel's `exp`, the one the model uses, because a reference built on a different
    /// rounding of `E` would differ from the model by `ε/(1 − E)` of itself for a reason that is
    /// not a defect. `expm1` is the platform's; its error is an ulp and is inside the bound.
    fn stepped_exact(f: f64, ka: f64, k: f64, v: f64, h: f64, steps: usize) -> Vec<f64> {
        let e = pantometry_core::math::exp(-ka * h);
        let (kh, r) = (k * h, 1.0 - k * h);
        let (mut series, mut power) = (0.0, kh);
        for m in 2..=40 {
            power *= kh;
            series += power / m as f64;
        }
        let l = (k - ka) * h + series;
        let amplitude = f * ORAL * (1.0 - e) / v;
        let mut out = vec![0.0];
        let mut r_power = 1.0; // r^{n−1}
        for n in 1..=steps {
            let geometric = if l == 0.0 {
                n as f64
            } else {
                (n as f64 * l).exp_m1() / l.exp_m1()
            };
            out.push(amplitude * r_power * geometric);
            r_power *= r;
        }
        out
    }

    /// The worst disagreement between a stepped curve and [`stepped_exact`], in units of a
    /// rounding bound, asserted to be under one. The bound at step `n` is
    ///
    /// ```text
    /// (3n + 40) ε Cₙ  +  2ε (F·D/V) Σ_{j<n} Eʲ
    /// ```
    ///
    /// The second term is the depot's, and is the larger one at small `n` and small `kₐh`. The
    /// step computes `kept = remaining·E` and hands over `remaining − kept`: the product rounds by
    /// up to `ε·remaining/2`, and the subtraction is exact, so the *absorbed* amount carries that
    /// whole error — `ε/(2(1 − E))` of itself, some 130 ulps at `kₐh = 3.75e-3`. Step `j` injects
    /// at most that much, `ε·F·D·Eʲ/(2V)` in concentration, propagated by powers of `r ≤ 1`,
    /// which sums to under the term as written; the factor of four over that covers the rounding
    /// of `E` itself. The depot does not leak by it — what one side lost the other gained, to the
    /// bit. Found by this test, which failed at step 2 for `h = 9 s` at 1.3e-14 of the value
    /// before the term was derived.
    ///
    /// The first term, term by term: the model's step rounds about twice per step relative to the
    /// amount it carries (the clearance product, the addition), `2nε`; the reference's `r^{n−1}`
    /// is `n` multiplications at half an ulp each, `nε/2`, rounded up to `nε`; `expm1(nL)` carries
    /// `|nL| ≤ 12` times `L`'s rounding, and `L` itself inherits `E`'s half-ulp, which with the
    /// remaining handful of products is the `40ε`.
    fn on_the_recurrence(c: &[f64], f: f64, ka: f64, k: f64, v: f64, h: f64) -> f64 {
        let exact = stepped_exact(f, ka, k, v, h, c.len() - 1);
        let e = pantometry_core::math::exp(-ka * h);
        // Σ_{j<n} Eʲ, kept beside the loop: see the second term of the bound.
        let (mut absorbed_sum, mut e_power) = (1.0, e);
        let (mut worst, mut worst_absolute, mut peak) = (0.0f64, 0.0f64, 0.0f64);
        for (n, (&got, &want)) in c.iter().zip(&exact).enumerate() {
            if n == 0 {
                assert_eq!(got, 0.0);
                continue;
            }
            let bound = (3.0 * n as f64 + 40.0) * f64::EPSILON * want
                + 2.0 * f64::EPSILON * f * ORAL / v * absorbed_sum;
            absorbed_sum += e_power;
            e_power *= e;
            let ratio = (got - want).abs() / bound;
            assert!(
                ratio < 1.0,
                "kₐ = {}/h, h = {h} s, step {n}: {got:e} against the recurrence's {want:e}, {:.2e} of the value",
                ka * HOUR,
                (got / want - 1.0).abs()
            );
            worst = worst.max(ratio);
            worst_absolute = worst_absolute.max((got - want).abs());
            peak = peak.max(want);
        }
        println!(
            "    worst {:.2e} of the peak concentration",
            worst_absolute / peak
        );
        worst
    }

    /// **The model steps exactly the recurrence its scheme defines, to rounding.**
    ///
    /// The convergence tests below say the stepped curve approaches the Bateman function at first
    /// order, and an order band cannot see a wrong error *constant*: an Euler depot, or absorbed
    /// drug eliminated in the step it arrives, still converges at order one. This pins the
    /// constant. [`stepped_exact`] is the closed-form solution of the discrete map — derived
    /// there — and the model must agree with it at every point to `(3n + 40)ε`, at four steps
    /// and three absorption rates: faster than elimination, equal to it, and slower (flip-flop).
    /// Measured: the worst point of each of the twelve runs is 0.05 to 0.18 of that bound, and
    /// the worst disagreement anywhere is printed as a fraction of `C_max` too. An Euler depot
    /// is off by about `kₐh/2` relative, and absorption folded into the snapshot by about `kh`,
    /// each ten orders of magnitude over it.
    #[test]
    fn the_stepped_curve_is_the_recurrences_exact_solution() {
        let (f, k, v) = (0.8, 0.5 / HOUR, 0.01);
        for ka_h in [1.5, 0.5, 0.1] {
            for h in [72.0, 36.0, 18.0, 9.0] {
                let (mut m, plasma, gut) = oral_one(ka_h, f, 0.0);
                m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
                let c = trace(&mut m, plasma, h, (12.0 * HOUR / h) as usize);
                let worst = on_the_recurrence(&c, f, ka_h / HOUR, k, v, h);
                println!("  ka = {ka_h}/h, h = {h:>4} s: worst {worst:.3} of the rounding bound");
            }
        }
    }

    /// **Oral closed form 1: the Bateman function, and its peak, at a measured first order.**
    ///
    /// `C(t) = F·D·kₐ/(V(kₐ − k))·(e^{−kt} − e^{−kₐt})`, with `kₐ = 1.5/h`, `k = 0.5/h` and
    /// `F = 0.8`, so `t_max = ln(kₐ/k)/(kₐ − k) = ln 3` hours.
    ///
    /// What the stepped curve *is* is pinned by the recurrence test above, to rounding; this one
    /// says what it converges *on*. Three things converge and each is asserted as a **rate**, at
    /// six steps each half the last, 72 s to 2.25 s: the whole curve's worst error over twelve
    /// hours, relative to `C_max`; the peak time, located as the discrete maximum and refined by
    /// the parabola through it and its two neighbours; and the peak height.
    ///
    /// Measured, at `h` = 72 s down to 2.25 s: the curve's worst error 3.76e-3 down to 1.16e-4
    /// of `C_max`, orders 1.0084, 1.0043, 1.0021, 1.0011, 1.0005; `C_max` orders 1.0091, 1.0046,
    /// 1.0023, 1.0012, 1.0006; `t_max` orders 0.9586, 0.9779, 0.9961, 0.9977, 0.9974, its
    /// discrete peak at 1.09854 h against `ln 3 = 1.09861` at the finest step.
    ///
    /// **The bands are earned per pair, not once.** The departure of a measured order from 1 is
    /// the `O(h²)` term's share of the error, which is proportional to `h`, so the band for the
    /// pair starting at `h` is `slack·(h/72 s)` — it halves as the step does. For the curve and
    /// the peak height `slack = 0.02`, and their departures (0.0084 halving to 0.0005) sit at
    /// 0.40 to 0.46 of it at every pair. `t_max` is located on a grid and refined by a parabola,
    /// so its departures halve only while the `O(h²)` share dominates — 0.041, 0.022 — and then
    /// settle at a floor near 0.003 from the location itself; `slack = 0.08` puts its five
    /// departures at 0.52, 0.55, 0.19, 0.23 and 0.52 of their bands. The single band of 0.05 this
    /// replaced was only 1.2 times the coarsest departure, which is the margin the review flagged.
    /// Each is also required to *shrink* from the coarsest pair to the finest, which a constant
    /// bias in the order would not, and the error is required to be real at the finest step,
    /// which a convergence test that has stopped measuring would fail.
    #[test]
    fn an_oral_dose_follows_the_bateman_function_at_first_order() {
        let (ka_h, f) = (1.5, 0.8);
        let (ka, k, v) = (ka_h / HOUR, 0.5 / HOUR, 0.01);
        let t_max = (ka / k).ln() / (ka - k);
        let c_max = bateman_textbook(f, ka, k, v, t_max);
        assert!(
            (t_max / HOUR - 3.0f64.ln()).abs() < 1e-12,
            "{}",
            t_max / HOUR
        );
        // The closed form's own peak, independent of any code under test: its derivative
        // vanishes there, `k e^{−k t} = kₐ e^{−kₐ t}`.
        assert!(((-k * t_max).exp() * k / ((-ka * t_max).exp() * ka) - 1.0).abs() < 1e-12);

        let (mut curve, mut peak_t, mut peak_c) = (Vec::new(), Vec::new(), Vec::new());
        for h in [72.0, 36.0, 18.0, 9.0, 4.5, 2.25] {
            let (mut m, plasma, gut) = oral_one(ka_h, f, 0.0);
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let c = trace(&mut m, plasma, h, (12.0 * HOUR / h) as usize);
            let worst = c
                .iter()
                .enumerate()
                .map(|(n, &got)| (got - bateman_textbook(f, ka, k, v, n as f64 * h)).abs())
                .fold(0.0, f64::max);
            curve.push(worst / c_max);

            let n = (1..c.len() - 1)
                .max_by(|&a, &b| c[a].total_cmp(&c[b]))
                .unwrap();
            let (y0, y1, y2) = (c[n - 1], c[n], c[n + 1]);
            let off = 0.5 * (y0 - y2) / (y0 - 2.0 * y1 + y2);
            assert!(
                off.abs() <= 0.5,
                "the parabola's vertex lies between the neighbours"
            );
            let (t_peak, c_peak) = ((n as f64 + off) * h, y1 - 0.25 * (y0 - y2) * off);
            peak_t.push((t_peak - t_max).abs() / t_max);
            peak_c.push((c_peak - c_max).abs() / c_max);
            println!(
                "  h = {h:>4} s: curve {:.3e}, t_max {:.5} h ({:.3e}), C_max {:.6e} ({:.3e})",
                worst / c_max,
                t_peak / HOUR,
                (t_peak - t_max).abs() / t_max,
                c_peak,
                (c_peak - c_max).abs() / c_max
            );
        }
        for (what, errors, slack) in [
            ("curve", &curve, 0.02),
            ("C_max", &peak_c, 0.02),
            ("t_max", &peak_t, 0.08),
        ] {
            let measured = orders(errors);
            println!("  {what}: orders {measured:.4?}");
            // The band shrinks with the step, because the departure it bounds is the `O(h²)`
            // share and is proportional to `h`: `slack·(h/72 s)` for the pair starting at `h`.
            let mut h = 72.0;
            for p in &measured {
                let band = slack * h / 72.0;
                h /= 2.0;
                assert!(
                    (p - 1.0).abs() < band,
                    "{what} should converge at first order and measured {p:.4}: {errors:?}"
                );
            }
            // The departure from 1 is the `O(h²)` term's share, so it must itself shrink as the
            // step does — the finest pair closer to 1 than the coarsest, by at least half.
            let (first, last) = (measured[0] - 1.0, measured[measured.len() - 1] - 1.0);
            assert!(
                last.abs() < 0.5 * first.abs(),
                "{what}: the order should approach 1 as h falls, and went {measured:.4?}"
            );
            assert!(*errors.last().unwrap() > 1e-5, "{what} is above the floor");
        }
    }

    /// **The depot is integrated exactly, to rounding — not to first order.**
    ///
    /// `A_d(t) = F·D·e^{−kₐt}`, and each step multiplies by `e^{−kₐh}` from the kernel's `exp`,
    /// whose error is under 0.52 ulp. Over `n` steps that one rounding repeats `n` times and the
    /// products add half an ulp each, so the bound is `n·1.02ε` plus the reference's own
    /// half-ulp; asserted at `1.1·n·ε`, and measured at 0.28 to 0.50 `n·ε`. Explicit Euler on the
    /// depot would be `(1 − kₐh)ⁿ`, low by a factor `e^{−t kₐ² h/2}` — 24% at twelve hours and
    /// `h = 72 s` — and the test is shaped to tell the two apart by twelve orders of magnitude.
    /// The error is also the same at every step size, so it is not a convergence at all: halving
    /// `h` doubles `n` and so doubles the bound, and the measurement stays under it.
    #[test]
    fn the_depot_is_exact_to_rounding_at_any_step() {
        let (ka_h, f) = (1.5, 0.8);
        let ka = ka_h / HOUR;
        for h in [72.0, 9.0, 1800.0] {
            let (mut m, _, gut) = oral_one(ka_h, f, 0.0);
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let mut bus = Exchange::new();
            let mut worst: f64 = 0.0;
            let steps = (12.0 * HOUR / h) as usize;
            for n in 0..steps {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                let t = (n + 1) as f64 * h;
                let exact = f * ORAL * (-ka * t).exp();
                let rel = (m.depot_amount(gut).to_si() / exact - 1.0).abs();
                let bound = 1.1 * (n + 1) as f64 * f64::EPSILON;
                assert!(
                    rel <= bound,
                    "h = {h}: after {} steps the depot is off by {rel:e}, over {bound:e}",
                    n + 1
                );
                worst = worst.max(rel / ((n + 1) as f64 * f64::EPSILON));
            }
            println!("  h = {h:>6} s: worst depot error {worst:.3} n eps");
        }
        // And a step as long as twelve absorption half-lives is still exact, and still stable:
        // the depot adds nothing to the step limit, which stays the compartment's own `1/k`.
        let (m, plasma, _) = oral_one(100.0, f, 0.0);
        assert_eq!(
            m.max_stable_dt(Time::ZERO).to_si(),
            1.0 / m.elimination_rate(plasma).to_si()
        );
    }

    /// **Oral closed form 2: `kₐ = k`, the singular limit, in the helper and in the model.**
    ///
    /// The textbook Bateman function is `0/0` at `kₐ = k`; its limit is `F·D·k·t·e^{−kt}/V`.
    ///
    /// **The helper.** [`bateman`] at equality must be that limit to rounding, and near it the
    /// relative change must be what the expansion says: at `kₐ = k(1 + δ)`,
    /// `ln C − ln C₀ = δ(1 − kt/2) + O(δ²)`. At `δ = ±1e-8` the `O(δ²)` is `1e-16`, so the
    /// assertion is that the helper is within `8ε` of the first-order expansion — where the
    /// textbook spelling, measured here too, is off by around `1e-8`, because its difference of
    /// two exponentials has kept eight digits. At `kt = 2` the linear term vanishes and the
    /// helper must not move at all beyond rounding. Far from equality it must agree with the
    /// textbook form, which is accurate there, to the textbook's own rounding.
    ///
    /// **The model.** It divides by nothing — the depot is an exponential and the compartment is
    /// Euler — so `kₐ = k` is an ordinary case for it. It must be the exact solution of its own
    /// recurrence at every point, to rounding, and converge on the limit at first order like any
    /// other case: measured at orders 1.0053 and 1.0027.
    #[test]
    fn equal_absorption_and_elimination_is_the_limit_not_a_division_by_zero() {
        let (f, v) = (0.8, Volume::litres(10.0));
        let k = 0.5 / HOUR;
        let dose = Mass::from_si(ORAL);
        let limit = |t: f64| f * ORAL * k * t * (-k * t).exp() / v.to_si();
        let near_rel = |ka: Frequency, t: Time, at: f64| {
            bateman(dose, f, ka, Frequency::from_si(k), v, t).to_si() / at - 1.0
        };
        for kt in [1e-3, 0.5, 2.0, 5.0, 40.0] {
            let t = Time::s(kt / k);
            let at = bateman(dose, f, Frequency::from_si(k), Frequency::from_si(k), v, t).to_si();
            let want = limit(t.to_si());
            assert!(
                (at / want - 1.0).abs() < 8.0 * f64::EPSILON,
                "at kt = {kt}: {at:e} against the limit {want:e}"
            );
            for asked in [1e-8, -1e-8] {
                let ka = Frequency::from_si(k * (1.0 + asked));
                // The separation actually realised, which is not the one asked for: `k(1 + δ)`
                // rounds to an ulp of `k`, so `kₐ − k` is `kδ` only to `ε/δ ≈ 2e-8` of itself.
                // The subtraction is exact (Sterbenz), so this is the δ the helper is handed.
                let delta = (ka.to_si() - k) / k;
                // ln(C/C₀) = ln(1 + δ) + ln φ(δkt) = δ − δ²/2 − δkt/2 + δ²(kt)²/24 + O(δ³), the
                // same for either sign; and e^L − 1 = L + L²/2. At kt = 40 the second-order
                // terms are 2.5e-14, a hundred ulps of the answer, so they are kept.
                let a = delta * (1.0 - kt / 2.0);
                let expected = a + a * a / 2.0 + delta * delta * (kt * kt / 24.0 - 0.5);
                let rel = near_rel(ka, t, at);
                // Rounding: a few ulps for the products, and `kt·ε` because an exponential of
                // an argument of size `kt` inherits that argument's rounding times `kt`.
                assert!(
                    (rel - expected).abs() < (8.0 + kt) * f64::EPSILON,
                    "kt = {kt}, delta = {delta:e}: moved {rel:e}, the expansion says {expected:e}"
                );
                let textbook = bateman_textbook(f, ka.to_si(), k, v.to_si(), t.to_si());
                if kt == 0.5 {
                    println!(
                        "  kt = 0.5, delta = {delta:e}: helper off the expansion by {:.1e}, \
                         textbook by {:.1e}",
                        (rel - expected).abs(),
                        (textbook / at - 1.0 - expected).abs()
                    );
                }
            }
        }
        // Far from equality the textbook form is the closed form to rounding, and the helper is
        // it — both orderings, flip-flop included.
        for ka_h in [0.05, 0.2, 1.5, 6.0] {
            for t_h in [0.1, 1.0, 3.0, 10.0] {
                let got = bateman(
                    dose,
                    f,
                    per_hour(ka_h),
                    Frequency::from_si(k),
                    v,
                    hours(t_h),
                )
                .to_si();
                let want = bateman_textbook(f, ka_h / HOUR, k, v.to_si(), t_h * HOUR);
                assert!(
                    (got / want - 1.0).abs() < 1e-13,
                    "ka = {ka_h}/h at {t_h} h: {got:e} against {want:e}"
                );
            }
        }
        assert_eq!(
            bateman(dose, f, per_hour(1.0), per_hour(0.5), v, Time::ZERO).to_si(),
            0.0
        );

        // The model at kₐ = k, converging on the limit at first order.
        let c_peak = limit(1.0 / k);
        let mut errors = Vec::new();
        for h in [72.0, 36.0, 18.0] {
            let (mut m, plasma, gut) = oral_one(0.5, f, 0.0);
            m.oral(gut, dose, Time::ZERO).unwrap();
            let c = trace(&mut m, plasma, h, (12.0 * HOUR / h) as usize);
            let worst = c
                .iter()
                .enumerate()
                .map(|(n, &got)| (got - limit(n as f64 * h)).abs())
                .fold(0.0, f64::max);
            assert!(c.iter().all(|x| x.is_finite()));
            // Not only an order: the curve is the discrete map's own solution at kₐ = k, where
            // `E ≠ r` by `(kh)²/2` and the reference must not cancel — see `stepped_exact`.
            let pinned = on_the_recurrence(&c, f, k, k, v.to_si(), h);
            println!("  ka = k, h = {h} s: {pinned:.3} of the rounding bound");
            errors.push(worst / c_peak);
        }
        let measured = orders(&errors);
        println!("  ka = k: errors {errors:?}, orders {measured:.4?}");
        for p in measured {
            assert!((p - 1.0).abs() < 0.02, "first order at ka = k: {p:.4}");
        }
    }

    /// **Oral closed form 3: flip-flop. When absorption is slower, the terminal slope is `−kₐ`.**
    ///
    /// With `kₐ = 0.1/h < k = 0.5/h` the curve's tail is absorption-limited: `C ∝ e^{−kₐt}` once
    /// `e^{−(k − kₐ)t}` has gone. Measured as `Δln C/Δt` between 80 and 120 hours (below).
    ///
    /// **And it is exact, not first order**, which is a consequence worth stating. The depot is
    /// integrated exactly, so the compartment is a discrete linear map driven by `E = e^{−kₐh}`;
    /// its particular solution decays by exactly `E` a step and its homogeneous one by `1 − kh`,
    /// which is gone. So the measured slope is `ln E/h = −kₐ` to rounding — where an Euler depot
    /// would read `ln(1 − kₐh)/h`, off by `kₐh/2 = 5e-4`. The ordinary case, `kₐ = 1.5/h`, is the
    /// contrast: its terminal slope is the compartment's `−k`, and that *is* Euler's,
    /// `ln(1 − kh)/h`, high in magnitude by `kh/2 = 2.5e-3`.
    ///
    /// The window is 80 to 120 hours, where the faster exponential is `e^{−32}` of the slower.
    /// The bound is rounding: the 4 000 steps between the two samples each carry about `ε` of
    /// relative error into `C`, which a difference of logarithms over `Δt` turns into
    /// `4000ε/(kₐΔt) = 2.2e-13` of the slope. Measured: 3.4e-14.
    #[test]
    fn a_slow_absorption_sets_the_terminal_slope() {
        let h = 36.0;
        let slope = |ka_h: f64| {
            let (mut m, plasma, gut) = oral_one(ka_h, 0.8, 0.0);
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let c = trace(&mut m, plasma, h, (120.0 * HOUR / h) as usize);
            let (n1, n2) = ((80.0 * HOUR / h) as usize, (120.0 * HOUR / h) as usize);
            (c[n2].ln() - c[n1].ln()) / ((n2 - n1) as f64 * h)
        };
        let ka = 0.1 / HOUR;
        let k = 0.5 / HOUR;
        let bound = 4000.0 * f64::EPSILON / (ka * 40.0 * HOUR);
        let flip = slope(0.1);
        println!("  flip-flop slope / -ka - 1 = {:.3e}", flip / -ka - 1.0);
        assert!(
            (flip / -ka - 1.0).abs() < bound,
            "the terminal slope should be -ka = {:e}/s, got {flip:e}",
            -ka
        );
        assert!((flip / -k - 1.0).abs() > 0.5, "and it is not -k");

        let ordinary = slope(1.5);
        let euler = (1.0 - k * h).ln() / h;
        println!(
            "  ordinary slope / -k - 1 = {:.3e}, Euler's own {:.3e}, the two {:.2e} apart",
            ordinary / -k - 1.0,
            euler / -k - 1.0,
            ordinary / euler - 1.0
        );
        // The same rounding bound, now relative to `k` rather than `kₐ`, which is smaller.
        assert!(
            (ordinary / euler - 1.0).abs() < bound,
            "{ordinary:e} against {euler:e}"
        );
        assert!((ordinary / -k - 1.0).abs() < 1.1 * k * h / 2.0);
    }

    /// **Oral closed form 4: `AUC₀^∞ = F·D/CL`, whatever `kₐ` is — and what that pins.**
    ///
    /// Integrated numerically: the left Riemann sum `Σ Cₙh` of the stepped curve out to 240 hours,
    /// where the slowest tail here (`kₐ = 0.2/h`) is `e^{−48}` of its start. Three absorption
    /// rates, flip-flop included, at two step sizes.
    ///
    /// **This is the mass-balance identity, and it is asserted as one.** The Euler step removes
    /// exactly `CL·Cₙ·h` each step, so `CL·Σ Cₙh` *is* the eliminated mass, and the books make
    /// that `F·D` less what is still in the body and the depot. So the area is asserted at `n·ε`
    /// with `n` the step count, and does **not** halve with the step: it was never a truncation
    /// error. What it pins is therefore what the identity contains — `F`, and that the samples
    /// `Cₙ` are the snapshot elimination acts on (absorbed drug folded into the snapshot first
    /// breaks it). It says **nothing about timing**: an Euler depot, which absorbs at the wrong
    /// rate, delivers the same total and passes it at 4.1e-15. Timing is what the recurrence and
    /// convergence tests are for. Measured: 2.7e-15 to 1.1e-14 at `h = 72 s` and
    /// 1.0e-14 to 2.6e-14 at 36 s — growing with the step count, as rounding does, and two orders
    /// under `n·ε`. `F` applied twice reads `F²D/CL`, 20% low.
    #[test]
    fn the_area_under_the_curve_is_the_bioavailable_dose_over_the_clearance() {
        let f = 0.8;
        let cl = Clearance::l_per_h(5.0).to_si();
        let want = f * ORAL / cl;
        for h in [72.0, 36.0] {
            let steps = (240.0 * HOUR / h) as usize;
            for ka_h in [0.2, 1.5, 6.0] {
                let (mut m, plasma, gut) = oral_one(ka_h, f, 0.0);
                m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
                let c = trace(&mut m, plasma, h, steps);
                let auc: f64 = c[..steps].iter().map(|x| x * h).sum();
                let rel = (auc / want - 1.0).abs();
                println!(
                    "  h = {h} s, ka = {ka_h}/h: AUC off by {rel:.2e} ({:.2} n eps)",
                    rel / (steps as f64 * f64::EPSILON)
                );
                assert!(
                    rel < steps as f64 * f64::EPSILON,
                    "h = {h}, ka = {ka_h}/h: AUC {auc:e} against F D/CL = {want:e}"
                );
            }
        }
    }

    /// **A lag shifts the curve by exactly `t_lag` — to the bit when the grid allows it.**
    ///
    /// On the grid: a lag of 1800 s at `h = 36 s` is fifty steps, and every time involved is a
    /// whole number of seconds, so the lagged run must be the unlagged one shifted by fifty steps
    /// **bit for bit** — every concentration, the depot and the cleared mass — and exactly zero
    /// before. The step that starts the absorption takes `h` itself as its span, not
    /// `until − start`, which is what lets this be an `assert_eq`.
    ///
    /// Off the grid: a lag of 1813 s starts absorption half-way through a step, and **the depot
    /// does not care** — it is `F·D·e^{−kₐ(t − t_lag)}` at every grid point after the lag, to
    /// rounding, because the first span is `until − start` exactly. Where the grid *does* matter is
    /// the compartment: its Euler step sees the onset as a slice inside one step, and the curve
    /// against the shifted Bateman function is first order, measured at a ratio near 2 between
    /// `h = 36` and `18 s` — the same order the unlagged curve converges at.
    ///
    /// And a dose *taken* off the grid, with no lag, is the same run as a lag ending there: see
    /// the last block, which is what notices a dose taken one step late.
    #[test]
    fn a_lag_shifts_the_curve_and_nothing_else() {
        let (ka_h, f, h) = (1.5, 0.8, 36.0);
        let steps = (12.0 * HOUR / h) as usize;
        let lagged_steps = 50usize;
        let lag = lagged_steps as f64 * h;
        let run_with = |lag: f64| {
            let (mut m, plasma, gut) = oral_one(ka_h, f, lag);
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let mut bus = Exchange::new();
            let mut out = vec![(0.0, m.depot_mass().to_si(), 0.0)];
            for n in 0..steps {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                out.push((
                    m.amount(plasma).to_si(),
                    m.depot_mass().to_si(),
                    m.eliminated_mass().to_si(),
                ));
            }
            out
        };
        let (plain, late) = (run_with(0.0), run_with(lag));
        // Index 0 is before the first step, when the dose has not been taken: the depot is empty
        // and the dose is pending. From the first step it is in the depot, untouched, until the
        // lag is over.
        for (n, row) in late.iter().enumerate().take(lagged_steps + 1).skip(1) {
            assert_eq!(row.0, 0.0, "nothing absorbed before the lag, step {n}");
            assert_eq!(row.1, f * ORAL, "the depot is full until the lag, step {n}");
        }
        for n in lagged_steps + 1..=steps {
            assert_eq!(
                late[n],
                plain[n - lagged_steps],
                "step {n}: the lagged run is the plain one {lagged_steps} steps later"
            );
        }

        // Off the grid.
        let ka = ka_h / HOUR;
        let (k, v) = (0.5 / HOUR, 0.01);
        let odd = lag + 13.0;
        let mut errors = Vec::new();
        for h in [36.0, 18.0] {
            let (mut m, plasma, gut) = oral_one(ka_h, f, odd);
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let mut bus = Exchange::new();
            let (mut worst, mut depot_worst) = (0.0f64, 0.0f64);
            let steps = (12.0 * HOUR / h) as usize;
            for n in 0..steps {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                let t = (n + 1) as f64 * h;
                let exact = bateman_textbook(f, ka, k, v, t - odd);
                worst = worst.max((m.concentration(plasma).to_si() - exact).abs());
                let depot_exact = if t > odd {
                    f * ORAL * (-ka * (t - odd)).exp()
                } else {
                    f * ORAL
                };
                depot_worst =
                    depot_worst.max((m.depot_amount(gut).to_si() / depot_exact - 1.0).abs());
            }
            println!("  lag {odd} s, h = {h} s: curve {worst:.3e}, depot {depot_worst:.2e}");
            assert!(
                depot_worst < 1.1 * steps as f64 * f64::EPSILON,
                "the depot is exact across a lag that is not on the grid: {depot_worst:e}"
            );
            errors.push(worst);
        }
        // **A dose taken part-way through a step**, rather than a lag ending there: no lag, and
        // the dose at 1 813 s. The step from 1 800 s takes it and absorbs from 1 813 s on, so the
        // depot must be `F·D·e^{−kₐ(t − 1813)}` at every grid point after, and the whole run must
        // be the `lag = 1813, at = 0` run **to the bit** in everything that is not the schedule:
        // the two have the same absorption start, `at + lag`, and that is all the step reads.
        // Before 1 800 s the late dose is pending rather than in the depot, which is the one
        // difference, and is asserted rather than skipped.
        let h = 36.0;
        let steps = (12.0 * HOUR / h) as usize;
        let observe = |lag: f64, at: f64| {
            let (mut m, plasma, gut) = oral_one(ka_h, f, lag);
            m.oral(gut, Mass::from_si(ORAL), Time::s(at)).unwrap();
            let mut bus = Exchange::new();
            let mut out = Vec::new();
            for n in 0..steps {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                out.push((
                    m.amount(plasma).to_si(),
                    m.eliminated_mass().to_si(),
                    m.depot_amount(gut).to_si(),
                    m.pending_dose().to_si(),
                ));
            }
            out
        };
        let (lagged, taken_late) = (observe(odd, 0.0), observe(0.0, odd));
        let mut depot_worst: f64 = 0.0;
        for n in 0..steps {
            let t = (n + 1) as f64 * h;
            let (a, b) = (lagged[n], taken_late[n]);
            assert_eq!((a.0, a.1), (b.0, b.1), "step {n}: the body and the books");
            if t <= 1800.0 {
                assert_eq!((b.2, b.3), (0.0, ORAL), "step {n}: not taken yet");
            } else {
                assert_eq!(a.2, b.2, "step {n}: the depot");
                let exact = f * ORAL * (-ka * (t - odd)).exp();
                depot_worst = depot_worst.max((b.2 / exact - 1.0).abs());
            }
        }
        println!("  dose taken at {odd} s: depot {depot_worst:.2e}");
        assert!(depot_worst < 1.1 * steps as f64 * f64::EPSILON);

        let p = orders(&errors)[0];
        println!("  off-grid lag: order {p:.4}");
        // Measured 0.9997: 8.59e-6 then 4.30e-6 mg/L. The bound is the unlagged curve's 0.02.
        assert!(
            (p - 1.0).abs() < 0.02,
            "first order with the lag off the grid: {p:.4}"
        );
    }

    /// **Mass balance with every route at once, at every step, to machine precision.**
    ///
    /// Two eliminating compartments and a link; a bolus; an infusion on neither a step boundary
    /// nor an hour; two depots — one into each compartment, one with a lag and `F = 0.62`, the
    /// other with `F = 1` — and three oral doses, one at zero, one exactly on a step boundary and
    /// one between two. At every step,
    ///
    /// ```text
    /// committed = in compartments + eliminated + in depots + lost + pending
    /// ```
    ///
    /// read twice: from the accessors, and from [`Domain::ledger`] — which is the one the
    /// scheduler's audit reads, so a depot left out of it is caught here and not only by a
    /// scheduler. The bound is the existing test's `1e-12`, and its reasoning carries over: about
    /// `ε` relative per step, 14 400 steps, errors of both signs.
    ///
    /// Then the loss is checked by dose rather than by total: exactly `(1 − F)·D` per dose, and
    /// zero for the route with `F = 1`. And the readings sum to the dose on their own.
    #[test]
    fn mass_balance_holds_through_the_depot_and_the_loss() {
        let mut m = CompartmentModel::new("subject");
        let central = m.eliminating("central", Volume::litres(5.0), Clearance::l_per_h(10.0));
        let peripheral = m.eliminating("peripheral", Volume::litres(20.0), Clearance::l_per_h(1.0));
        m.link(central, peripheral, Clearance::l_per_h(30.0))
            .unwrap();
        let gut = m
            .depot("gut", central, per_hour(1.3), 0.62, hours(0.21))
            .unwrap();
        let muscle = m
            .depot("muscle", peripheral, per_hour(0.4), 1.0, Time::ZERO)
            .unwrap();

        let (bolus, infused) = (1e-4, 2.5e-4);
        m.bolus(central, Mass::from_si(bolus)).unwrap();
        m.infuse(peripheral, Mass::from_si(infused), hours(0.37), hours(2.81))
            .unwrap();
        let oral = [(gut, 1e-4, 0.0), (gut, 2e-4, 4.0), (muscle, 3e-5, 1.234567)];
        for (route, d, at) in oral {
            m.oral(route, Mass::from_si(d), hours(at)).unwrap();
        }
        let committed = bolus + infused + oral.iter().map(|x| x.1).sum::<f64>();
        // Before anything has run, the whole oral regimen is pending.
        assert_eq!(m.pending_dose().to_si(), infused + 1e-4 + 2e-4 + 3e-5);

        let h = 6.0;
        let steps = (24.0 * HOUR / h) as usize;
        let mut bus = Exchange::new();
        let (mut worst, mut worst_ledger) = (0.0f64, 0.0f64);
        for n in 0..steps {
            m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
            let books = m.body_burden().to_si()
                + m.eliminated_mass().to_si()
                + m.depot_mass().to_si()
                + m.presystemic_loss().to_si()
                + m.pending_dose().to_si();
            worst = worst.max((books - committed).abs() / committed);
            let ledger = m.ledger().get(quantity::MASS).unwrap();
            worst_ledger = worst_ledger.max((ledger - committed).abs() / committed);
            // The 4-hour dose is pending until the step that reaches it, and taken by it.
            let t_end = (n + 1) as f64 * h;
            if t_end <= 4.0 * HOUR {
                assert!(m.pending_dose().to_si() >= 2e-4, "step {n}");
            } else {
                assert_eq!(m.pending_oral(), 0.0, "step {n}");
            }
        }
        println!("  books {worst:.2e}, ledger {worst_ledger:.2e} over {steps} steps");
        assert!(worst < 1e-12, "the books drifted by {worst:e} of the dose");
        assert!(
            worst_ledger < 1e-12,
            "the ledger drifted by {worst_ledger:e} of the dose"
        );

        // The loss, dose by dose: (1 − 0.62)·(1e-4 + 2e-4) through the gut, nothing through
        // the muscle.
        let lost = m.presystemic_loss().to_si();
        let expected = (1e-4 - 0.62 * 1e-4) + (2e-4 - 0.62 * 2e-4);
        assert_eq!(lost, expected);
        assert!((m.administered_mass().to_si() / committed - 1.0).abs() < 1e-15);
        // Not vacuous: drug went through both depots and most of it has been cleared.
        assert!(m.eliminated_mass().to_si() > 0.5 * committed);
        assert!(m.depot_amount(muscle).to_si() > 0.0);
        assert!(m.depot_amount(gut).to_si() < 1e-3 * 2e-4);

        // The readings carry the oral half of the books and sum to the dose by themselves.
        let r = m.readings();
        let mg = |label: &str| {
            r.iter()
                .find(|x| x.label == label)
                .unwrap_or_else(|| panic!("no {label} reading"))
                .value
        };
        let total = mg("in the body")
            + mg("cleared")
            + mg("still to infuse")
            + mg("still to take")
            + mg("in the depot")
            + mg("lost before absorption");
        assert!(
            (total / (committed * 1e6) - 1.0).abs() < 1e-12,
            "{total} mg"
        );
        assert!((mg("lost before absorption") / (expected * 1e6) - 1.0).abs() < 1e-15);
    }

    /// **No reading is a negative zero.** `Iterator::sum` over `f64` starts at `-0.0`, so a
    /// filter that keeps nothing — no dose taken yet, every dose taken — read `-0` and the CSV of
    /// the first scene with a tablet printed `still to take -0.000000000e0` beside the dose it
    /// had just swallowed. Asked before the dose, at the step it is taken, and of a depot no dose
    /// has gone through.
    #[test]
    fn an_empty_part_of_the_books_reads_plus_zero() {
        let (mut m, c) = one_compartment();
        let gut = m
            .depot("gut", c, Frequency::from_si(1.0 / HOUR), 0.8, Time::s(0.0))
            .unwrap();
        m.oral(gut, Mass::from_si(1e-4), Time::s(0.0)).unwrap();
        let signs = |m: &CompartmentModel| -> Vec<(String, f64)> {
            m.readings()
                .into_iter()
                .filter(|r| r.value == 0.0 && r.value.is_sign_negative())
                .map(|r| (r.label, r.value))
                .collect()
        };
        assert_eq!(signs(&m), [], "before the dose is taken");
        assert!(m.presystemic_loss().to_si().is_sign_positive());
        run(&mut m, Time::s(0.0), Time::s(60.0), 1);
        assert_eq!(signs(&m), [], "the step it is taken");
        assert!(m.pending_dose().to_si().is_sign_positive());
        // And a depot no dose has gone through, which is the empty filter in `depot_amount` and
        // the empty sum in `depot_mass`.
        let (mut n, d) = one_compartment();
        let unused = n
            .depot(
                "unused",
                d,
                Frequency::from_si(1.0 / HOUR),
                1.0,
                Time::s(0.0),
            )
            .unwrap();
        assert!(n.depot_amount(unused).to_si().is_sign_positive());
        assert!(n.depot_mass().to_si().is_sign_positive());
    }

    /// **Superposition: two doses are the sum of two single doses, to rounding.**
    ///
    /// The model is linear and the steps are the same, so the run with both doses and the sum of
    /// the two runs with one each differ only in the order roundings happen. Two compartments, a
    /// depot into the central one, 100 mg at zero and 37 mg at 5.3 hours — off the grid, and
    /// while the first is still absorbing.
    ///
    /// Measured against the largest amount of the run, the difference stays near `ε` per step
    /// and is asserted at `n·ε` for the step count `n`: a bound on independent roundings
    /// accumulating, of which `5.3e-13` at 2 400 steps is the worst case and `2.7e-15` is what
    /// is measured. A model that was not linear — a depot whose rate depended on its content, a
    /// second dose that reset the first — would be off at the size of a dose.
    #[test]
    fn two_oral_doses_superpose() {
        let h = 18.0;
        let steps = (12.0 * HOUR / h) as usize;
        let second_at = hours(5.3);
        let build = |first: bool, second: bool| {
            let (mut m, central, peripheral) = two_compartment();
            let gut = m
                .depot("gut", central, per_hour(1.2), 0.9, hours(0.1))
                .unwrap();
            if first {
                m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            }
            if second {
                m.oral(gut, Mass::from_si(3.7e-5), second_at).unwrap();
            }
            (m, central, peripheral)
        };
        let state = |m: &CompartmentModel, a: Compartment, b: Compartment| {
            [
                m.amount(a).to_si(),
                m.amount(b).to_si(),
                m.depot_mass().to_si(),
                m.eliminated_mass().to_si(),
            ]
        };
        let (mut one, c1, p1) = build(true, false);
        let (mut two, c2, p2) = build(false, true);
        let (mut both, c, p) = build(true, true);
        let mut bus = Exchange::new();
        let mut worst: f64 = 0.0;
        for n in 0..steps {
            let (t, dt) = (Time::s(n as f64 * h), Time::s(h));
            one.step(t, dt, &mut bus).unwrap();
            two.step(t, dt, &mut bus).unwrap();
            both.step(t, dt, &mut bus).unwrap();
            let (a, b, s) = (state(&one, c1, p1), state(&two, c2, p2), state(&both, c, p));
            for i in 0..4 {
                worst = worst.max((s[i] - (a[i] + b[i])).abs() / ORAL);
            }
        }
        println!("  superposition: worst {worst:.2e} of the dose over {steps} steps");
        assert!(
            worst < steps as f64 * f64::EPSILON,
            "two doses differ from the sum of one and one by {worst:e} of a dose"
        );
        // And the second dose did something: the depot was refilled part-way through.
        assert!(two.eliminated_mass().to_si() > 0.5 * 0.9 * 3.7e-5);
    }

    /// **Oral closed form 5: two compartments, oral input into the central one — the
    /// tri-exponential.**
    ///
    /// # The derivation
    ///
    /// With `A_d = F·D·e^{−kₐt}` feeding the central compartment,
    ///
    /// ```text
    /// A₁' = −(k10 + k12) A₁ + k21 A₂ + kₐ A_d
    /// A₂' =  k12 A₁ − k21 A₂
    /// ```
    ///
    /// Laplace-transformed from rest, the second gives `Ā₂ = k12 Ā₁/(s + k21)`, and substituting
    /// into the first, `Ā₁ [(s + k10 + k12)(s + k21) − k12 k21] = kₐ F D (s + k21)/(s + kₐ)`. The
    /// bracket is `s² + (k10 + k12 + k21)s + k10 k21 = (s + α)(s + β)`, the bolus case's
    /// characteristic polynomial, so
    ///
    /// ```text
    /// Ā₁ = F D kₐ (s + k21) / ((s + kₐ)(s + α)(s + β)),    Ā₂ = F D kₐ k12 / ((s + kₐ)(s + α)(s + β))
    /// ```
    ///
    /// and the residues at `−α`, `−β`, `−kₐ` give the three terms:
    ///
    /// ```text
    /// A₁(t) = F D kₐ [ (k21 − α)/((kₐ − α)(β − α)) e^{−αt} + (k21 − β)/((kₐ − β)(α − β)) e^{−βt}
    ///                 + (k21 − kₐ)/((α − kₐ)(β − kₐ)) e^{−kₐt} ]
    /// A₂(t) = F D kₐ k12 [ e^{−αt}/((kₐ − α)(β − α)) + e^{−βt}/((kₐ − β)(α − β))
    ///                    + e^{−kₐt}/((α − kₐ)(β − kₐ)) ]
    /// ```
    ///
    /// Checked before any code is: each set of coefficients must give `A(0) = 0`, `A₁'(0) = kₐFD`,
    /// `A₂'(0) = 0` and `A₂''(0) = k12 kₐ F D`, which the equations above say at `t = 0` and which
    /// a sign or an index wrong in a residue would break.
    ///
    /// # What it catches
    ///
    /// A depot feeding the peripheral compartment conserves every milligram and is wrong in both
    /// curves; a single-compartment test cannot see it because there is only one place to feed.
    /// The tolerance is a rate again: the worst error over twelve hours, relative to the dose,
    /// must halve with the step, for each compartment separately, at `h` = 18, 9 and 4.5 s.
    ///
    /// Measured: central 8.18e-4, 4.05e-4, 2.01e-4 of the dose, orders 1.0153 and 1.0078;
    /// peripheral 5.81e-4, 2.87e-4, 1.42e-4, orders 1.0178 and 1.0088 — a larger `O(h²)` share
    /// than one compartment's because the fast root `α = 9.17/h` makes `αh = 4.6e-2` at the
    /// coarse step. The band is `0.035·(h/18 s)`, halving with the pair as in the one-compartment
    /// test; the four departures sit at 0.44 to 0.51 of their bands.
    #[test]
    fn two_compartments_with_oral_input_follow_the_tri_exponential() {
        let (k10, k12, k21) = (2.0 / HOUR, 6.0 / HOUR, 1.5 / HOUR);
        let (ka, f) = (1.2 / HOUR, 0.9);
        let sum = k10 + k12 + k21;
        let disc = (sum * sum - 4.0 * k10 * k21).sqrt();
        let (alpha, beta) = (0.5 * (sum + disc), 0.5 * (sum - disc));
        let fd = f * ORAL;
        let rates = [alpha, beta, ka];
        let central = [
            (k21 - alpha) / ((ka - alpha) * (beta - alpha)),
            (k21 - beta) / ((ka - beta) * (alpha - beta)),
            (k21 - ka) / ((alpha - ka) * (beta - ka)),
        ];
        let peripheral = [
            1.0 / ((ka - alpha) * (beta - alpha)),
            1.0 / ((ka - beta) * (alpha - beta)),
            1.0 / ((alpha - ka) * (beta - ka)),
        ];
        // The derivation's own checks, in units that make each one a pure number.
        let moment = |c: &[f64; 3], p: i32| -> f64 {
            let mut s = 0.0;
            for i in 0..3 {
                let mut lam = 1.0;
                for _ in 0..p {
                    lam *= -rates[i];
                }
                s += c[i] * lam;
            }
            s
        };
        // Made dimensionless before comparing to 1e-12: a central coefficient is a time and a
        // peripheral one a time squared, so the `p`-th moment of the central set is in
        // `s^{1−p}` and of the peripheral set in `s^{2−p}`, and each is divided by the matching
        // power of an hour — the scale the rates are O(1) in. Each moment is then a sum of terms
        // of order one and 1e-12 is about 5 000 ulps of them.
        let central_moment = |p: i32| moment(&central, p) / [HOUR, 1.0][p as usize];
        let peripheral_moment =
            |p: i32| moment(&peripheral, p) / [HOUR * HOUR, HOUR, 1.0][p as usize];
        assert!(
            central_moment(0).abs() < 1e-12,
            "A1(0) = 0: {:e}",
            central_moment(0)
        );
        assert!((central_moment(1) - 1.0).abs() < 1e-12, "A1'(0) = ka F D");
        assert!(
            peripheral_moment(0).abs() < 1e-12,
            "A2(0) = 0: {:e}",
            peripheral_moment(0)
        );
        assert!(
            peripheral_moment(1).abs() < 1e-12,
            "A2'(0) = 0: {:e}",
            peripheral_moment(1)
        );
        assert!(
            (peripheral_moment(2) - 1.0).abs() < 1e-12,
            "A2''(0) = k12 ka F D"
        );

        let exact = |c: &[f64; 3], gain: f64, t: f64| {
            fd * ka * gain * (0..3).map(|i| c[i] * (-rates[i] * t).exp()).sum::<f64>()
        };
        let (mut errors_1, mut errors_2) = (Vec::new(), Vec::new());
        for h in [18.0, 9.0, 4.5] {
            let (mut m, c1, c2) = two_compartment();
            let gut = m
                .depot("gut", c1, Frequency::from_si(ka), f, Time::ZERO)
                .unwrap();
            m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let mut bus = Exchange::new();
            let (mut worst_1, mut worst_2) = (0.0f64, 0.0f64);
            for n in 0..(12.0 * HOUR / h) as usize {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                let t = (n + 1) as f64 * h;
                worst_1 = worst_1.max((m.amount(c1).to_si() - exact(&central, 1.0, t)).abs());
                worst_2 = worst_2.max((m.amount(c2).to_si() - exact(&peripheral, k12, t)).abs());
            }
            errors_1.push(worst_1 / ORAL);
            errors_2.push(worst_2 / ORAL);
        }
        // Each compartment on its own, because a depot feeding the wrong one is wrong in both and
        // a single maximum would let the larger hide which. Three resolutions, so each order is
        // two ratios and the band rests on more than one; and the departures must shrink.
        for (which, errors) in [("central", &errors_1), ("peripheral", &errors_2)] {
            let measured = orders(errors);
            println!("  tri-exponential, {which}: errors {errors:?}, orders {measured:.4?}");
            assert!(
                errors[0] < 0.01,
                "{which} near the closed form at all: {:e}",
                errors[0]
            );
            // `0.035·(h/18 s)` for the pair starting at `h`, the band shrinking with the
            // `O(h²)` share it bounds, as in the one-compartment test.
            let mut h = 18.0;
            for p in &measured {
                let band = 0.035 * h / 18.0;
                h /= 2.0;
                assert!(
                    (p - 1.0).abs() < band,
                    "{which}: first order, measured {p:.4}"
                );
            }
            assert!(
                (measured[1] - 1.0).abs() < (measured[0] - 1.0).abs(),
                "{which}: the order should approach 1 as h falls: {measured:.4?}"
            );
        }
    }

    /// **The same bits however the run is cut into calls.**
    ///
    /// Three cuts of one regimen — two compartments, a bolus, a lagged depot with two doses — and
    /// all must agree to the bit: eighteen scheduler windows of 400 s against nine of 800, which
    /// the multirate schedule splits into the same 400 s substeps (asserted, or the comparison is
    /// of two different runs and means nothing); and a run rewound by `restore` part-way through
    /// absorption and run again, which must land where the uninterrupted one did. The scheduler
    /// also audits the ledger every window at `1e-12`, so a depot missing from the books fails
    /// here as well as in the mass-balance test.
    #[test]
    fn an_oral_run_is_the_same_bits_however_it_is_cut() {
        let build = || {
            let (mut m, central, _) = two_compartment();
            m.bolus(central, Mass::from_si(5e-5)).unwrap();
            let gut = m
                .depot("gut", central, per_hour(1.7), 0.75, Time::s(500.0))
                .unwrap();
            m.oral(gut, Mass::from_si(ORAL), Time::s(130.0)).unwrap();
            m.oral(gut, Mass::from_si(6e-5), Time::s(3333.0)).unwrap();
            m
        };
        let state = |m: &CompartmentModel| -> Vec<f64> {
            let mut v: Vec<f64> = m.handles().map(|(c, _)| m.amount(c).to_si()).collect();
            v.push(m.depot_mass().to_si());
            v.push(m.eliminated_mass().to_si());
            v.push(m.presystemic_loss().to_si());
            v.push(m.pending_dose().to_si());
            v
        };
        let cut = |window: f64, windows: usize| {
            let mut sim = Simulation::new(Schedule::Multirate)
                .conservation_tolerance(1e-12)
                .with(build());
            let mut substeps = 0;
            for _ in 0..windows {
                let report = sim.advance(Time::s(window)).unwrap();
                substeps = report.substeps[0].1;
            }
            let m: &CompartmentModel = sim.domain_as("subject").unwrap();
            (state(m), substeps)
        };
        let (fine, n_fine) = cut(400.0, 18);
        let (coarse, n_coarse) = cut(800.0, 9);
        assert_eq!((n_fine, n_coarse), (1, 2), "both cuts take 400 s substeps");
        assert_eq!(fine, coarse, "two cuts of one run");
        assert!(
            fine[2] > 0.0 && fine[4] > 0.0,
            "the run is not vacuous: {fine:?}"
        );

        // Rewound part-way through absorption and run again.
        let mut straight = build();
        run(&mut straight, Time::ZERO, Time::s(400.0), 18);
        let mut rewound = build();
        run(&mut rewound, Time::ZERO, Time::s(400.0), 5);
        rewound.checkpoint();
        let mid = state(&rewound);
        assert!(
            mid[2] > 0.0,
            "the checkpoint is taken while the depot holds drug"
        );
        run(&mut rewound, Time::s(2000.0), Time::s(400.0), 13);
        rewound.restore();
        assert_eq!(state(&rewound), mid, "restore puts the depot back too");
        run(&mut rewound, Time::s(2000.0), Time::s(400.0), 13);
        assert_eq!(state(&rewound), state(&straight));
        assert_eq!(
            state(&straight),
            fine,
            "and the scheduler ran the same steps"
        );
    }

    /// **A bad route or a bad dose is refused where it is declared, and named.**
    ///
    /// Every `Violation` the oral half can raise: a `kₐ` that is zero, negative, infinite or NaN;
    /// an `F` of zero, below it, above one or NaN; a lag that is negative, infinite or NaN; a
    /// compartment or a depot from another model; a dose that is negative or NaN; a time that is
    /// not finite. And the edges that are allowed: `F = 1`, a lag of zero, and a zero dose that
    /// does nothing.
    #[test]
    fn a_bad_route_or_dose_is_refused_and_named() {
        let (mut m, plasma) = one_compartment();
        let ok = |m: &mut CompartmentModel, ka: f64, f: f64, lag: f64| {
            m.depot("gut", plasma, per_hour(ka), f, Time::s(lag))
        };
        for ka in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            let err = ok(&mut m, ka, 0.5, 0.0).expect_err("a bad absorption rate");
            assert!(err.quantity.contains("absorption rate"), "{err}");
            assert!(err.site.contains("gut"), "{err}");
        }
        for f in [0.0, -0.1, 1.000_000_1, f64::NAN] {
            let err = ok(&mut m, 1.0, f, 0.0).expect_err("a bad bioavailability");
            assert!(err.quantity.contains("bioavailability"), "{err}");
        }
        for lag in [-1.0, f64::INFINITY, f64::NAN] {
            let err = ok(&mut m, 1.0, 0.5, lag).expect_err("a bad lag");
            assert!(err.quantity.contains("lag"), "{err}");
        }
        let mut other = CompartmentModel::new("another subject");
        let elsewhere = other.eliminating("plasma", Volume::litres(10.0), Clearance::l_per_h(5.0));
        let err = m
            .depot("gut", elsewhere, per_hour(1.0), 0.5, Time::ZERO)
            .expect_err("a compartment from another model");
        assert!(err.quantity.contains("different model"), "{err}");
        let foreign = other
            .depot("gut", elsewhere, per_hour(1.0), 0.5, Time::ZERO)
            .unwrap();
        let err = m
            .oral(foreign, Mass::from_si(1e-4), Time::ZERO)
            .expect_err("a depot from another model");
        assert!(err.quantity.contains("different model"), "{err}");
        // Nothing refused above was half-built.
        assert!(m.depot_named("gut").is_none());

        let gut = ok(&mut m, 1.0, 1.0, 0.0).expect("F = 1 and no lag are allowed");
        assert_eq!(m.depot_named("gut"), Some(gut));
        assert_eq!(m.absorbed_into(gut), plasma);
        assert_eq!(m.bioavailability(gut), 1.0);
        assert_eq!(m.absorption_lag(gut).to_si(), 0.0);
        assert!((m.absorption_rate(gut).to_si() * HOUR - 1.0).abs() < 1e-15);
        for dose in [-1e-6, f64::NAN, f64::INFINITY] {
            let err = m
                .oral(gut, Mass::from_si(dose), Time::ZERO)
                .expect_err("a bad dose");
            assert!(err.quantity.contains("dose"), "{err}");
        }
        for at in [f64::INFINITY, f64::NAN] {
            let err = m
                .oral(gut, Mass::from_si(1e-4), Time::from_si(at))
                .expect_err("a dose at no time");
            assert!(err.quantity.contains("finite time"), "{err}");
        }
        assert_eq!(
            m.pending_dose().to_si(),
            0.0,
            "nothing refused was scheduled"
        );
        m.oral(gut, Mass::ZERO, Time::ZERO).unwrap();
        run(&mut m, Time::ZERO, Time::s(60.0), 10);
        assert_eq!(m.amount(plasma).to_si(), 0.0);
        assert_eq!(m.presystemic_loss().to_si(), 0.0);
    }

    /// **`φ(x) = (1 − e^{−x})/x`, the function [`bateman`] rests on, to a few ulps on both sides
    /// of its switch.**
    ///
    /// The Taylor branch is the one the singular-limit test only reaches near `x = 0`, where its
    /// first two terms are already the answer; cut to twelve terms, through `x¹¹`, instead of
    /// nineteen, every other test passes, and it is wrong at `x = ½` by `½¹²/13! ≈ 4e-14` —
    /// measured, 3.8e-14. So it is checked where its
    /// truncation is largest:
    ///
    /// - at the last double below ½, against `φ(½) = 2(1 − e^{−½}) = 0.786938680574733152792…`,
    ///   hard-coded from mpmath at 40 digits; and the exponential branch at ½ against the same
    ///   number. `φ′(½) ≈ −0.36`, so the half-ulp step between the two arguments moves `φ` by
    ///   `2e-17`, under a tenth of an ulp;
    /// - across `[0, ½]` against the series summed **forwards**, term by term, to forty terms —
    ///   the same mathematics as the nested form but a different evaluation, and with a
    ///   truncation of `½⁴⁰/41!`, nothing; the terms alternate and shrink, so the forward sum's
    ///   rounding is under two ulps of the result;
    /// - and against `−expm1(−x)/x` from the standard library, the platform's `expm1`, whose
    ///   error is an ulp — an independent function, which is why its bound is the loosest here.
    #[test]
    fn phi_is_accurate_on_both_sides_of_its_switch() {
        // The nearest double to 0.78693868057473315279…, 1.3e-18 below it.
        const PHI_HALF: f64 = 0.786_938_680_574_733_2;
        let below = f64::from_bits(0.5f64.to_bits() - 1);
        assert!(below < 0.5);
        let taylor = one_minus_exp_over(below);
        let exponential = one_minus_exp_over(0.5);
        let ulp = f64::EPSILON * PHI_HALF;
        println!(
            "  phi(1/2-) off by {:.2} ulp, phi(1/2) by {:.2} ulp, apart by {:.2} ulp",
            (taylor - PHI_HALF).abs() / ulp,
            (exponential - PHI_HALF).abs() / ulp,
            (taylor - exponential).abs() / ulp
        );
        assert!(
            (taylor - PHI_HALF).abs() <= 2.0 * ulp,
            "the Taylor branch at 1/2-: {taylor:e}"
        );
        assert!(
            (exponential - PHI_HALF).abs() <= 3.0 * ulp,
            "the exponential branch at 1/2: {exponential:e}"
        );
        assert!(
            (taylor - exponential).abs() <= 3.0 * ulp,
            "continuous across the switch: {taylor:e} and {exponential:e}"
        );
        assert_eq!(one_minus_exp_over(0.0), 1.0);

        let forward = |x: f64| {
            let (mut sum, mut term) = (0.0, 1.0);
            for n in 0..40 {
                // term = (−x)ⁿ/(n+1)!
                sum += term;
                term *= -x / (n + 2) as f64;
            }
            sum
        };
        let (mut worst_series, mut worst_expm1) = (0.0f64, 0.0f64);
        for i in 0..=4096 {
            let x = 0.5 * i as f64 / 4096.0;
            let x = if i == 4096 { below } else { x };
            let got = one_minus_exp_over(x);
            let series = (got - forward(x)).abs() / (got * f64::EPSILON);
            worst_series = worst_series.max(series);
            if x > 0.0 {
                let reference = -(-x).exp_m1() / x;
                worst_expm1 = worst_expm1.max((got - reference).abs() / (got * f64::EPSILON));
            }
        }
        println!("  [0, 1/2): against the forward series {worst_series:.2} ulp, against expm1 {worst_expm1:.2} ulp");
        assert!(
            worst_series <= 3.0,
            "{worst_series} ulp from the forward series"
        );
        assert!(worst_expm1 <= 4.0, "{worst_expm1} ulp from -expm1(-x)/x");
    }

    // ---- Perfusion-limited distribution ----------------------------------------------------

    fn l_h(x: f64) -> Clearance {
        Clearance::l_per_h(x)
    }

    /// The four-pool model the steady-state checks use: blood with a renal clearance, a liver
    /// that eliminates, a muscle and a fat that only distribute. **Round illustrative numbers,
    /// not physiology data**: chosen so `E = 0.4` and the arithmetic below is checkable by hand.
    ///
    /// | pool | V | Q | K_p | CL |
    /// | --- | --- | --- | --- | --- |
    /// | blood | 5 L | | | 6 L/h |
    /// | liver | 1.8 L | 90 L/h | 4 | 60 L/h, on `C_t/K_p` |
    /// | muscle | 30 L | 50 L/h | 2 | |
    /// | fat | 10 L | 10 L/h | 5 | |
    struct Body {
        m: CompartmentModel,
        blood: Compartment,
        liver: Compartment,
        muscle: Compartment,
        fat: Compartment,
    }

    fn body(liver_cl: f64) -> Body {
        let mut m = CompartmentModel::new("subject");
        let blood = m.eliminating("blood", Volume::litres(5.0), l_h(6.0));
        let liver = m.eliminating("liver", Volume::litres(1.8), l_h(liver_cl));
        let muscle = m.compartment("muscle", Volume::litres(30.0));
        let fat = m.compartment("fat", Volume::litres(10.0));
        m.perfuse(blood, liver, l_h(90.0), 4.0).unwrap();
        m.perfuse(blood, muscle, l_h(50.0), 2.0).unwrap();
        m.perfuse(blood, fat, l_h(10.0), 5.0).unwrap();
        Body {
            m,
            blood,
            liver,
            muscle,
            fat,
        }
    }

    /// The closed-form steady state of [`body`] under an infusion `R` into the blood:
    /// `[C_b, C_liver, C_muscle, C_fat]`, in kg/m³.
    ///
    /// Blood: `R = CL_r C_b + CL·C_l/K_l`. Liver: `Q_h(C_b − C_l/K_l) = CL·C_l/K_l`, so
    /// `C_l/K_l = C_b·Q_h/(Q_h + CL) = C_b(1 − E)`. A tissue with no clearance has no net flux,
    /// so `C_t = K_p C_b`. Substituting, `R = (CL_r + Q_h E) C_b`.
    fn body_steady_state(rate: f64, liver_cl: f64) -> [f64; 4] {
        let (q_h, cl, cl_r) = (l_h(90.0).to_si(), l_h(liver_cl).to_si(), l_h(6.0).to_si());
        let e = cl / (q_h + cl);
        let c_b = rate / (cl_r + q_h * e);
        [c_b, 4.0 * c_b * (1.0 - e), 2.0 * c_b, 5.0 * c_b]
    }

    /// **PBPK closed form 1: the steady state under a constant infusion — `C_t = K_p C_b` in every
    /// tissue that does not eliminate, `C_b = R/(CL_r + Q_h E)`, and the liver at
    /// `K_p C_b (1 − E)` — reproduced by Euler exactly, at any step.**
    ///
    /// The proof is the one the one-compartment infusion test states. The step is
    /// `A ← A + h(MA + R)`, linear, and its fixed point is `MA* = −R`: the continuous steady
    /// state, with `h` nowhere in it. So the claim is not "close at small `h`" but "the same
    /// point", and it is measured that way: the model is **started at** the closed form — each
    /// pool bolused with its steady-state amount — and run, and must stay there to rounding at
    /// two steps, one of them 0.98 of the stability limit. A transient has nothing to decay
    /// from, so no tolerance is spent on one.
    ///
    /// Then from empty, where the transient is real: the slowest mode must die and leave the same
    /// point at the rounding floor.
    #[test]
    fn under_an_infusion_every_tissue_settles_at_kp_times_the_blood() {
        let rate = 1e-3 / (2000.0 * HOUR); // 1 g over 2000 h: 0.5 mg/h
        let want = body_steady_state(rate, 60.0);
        for h in [60.0, 0.98 * 3600.0 / 31.2] {
            let mut b = body(60.0);
            let pools = [b.blood, b.liver, b.muscle, b.fat];
            for (c, x) in pools.iter().zip(want) {
                let a = x * b.m.volume(*c).to_si();
                b.m.bolus(*c, Mass::from_si(a)).unwrap();
            }
            b.m.infuse(b.blood, Mass::from_si(1e-3), Time::ZERO, hours(2000.0))
                .unwrap();
            let start: Vec<f64> = pools.iter().map(|c| b.m.amount(*c).to_si()).collect();
            let steps = 4000;
            let mut bus = Exchange::new();
            let (mut worst, mut ratio, mut moved) = (0.0f64, 0.0f64, 0usize);
            for n in 0..steps {
                b.m.step(Time::s(n as f64 * h), Time::s(h), &mut bus)
                    .unwrap();
                for (i, (c, x)) in pools.iter().zip(want).enumerate() {
                    let rel = (b.m.concentration(*c).to_si() / x - 1.0).abs();
                    worst = worst.max(rel);
                    ratio = ratio.max(rel / ((n as f64 + 2.0) * f64::EPSILON));
                    moved += usize::from(b.m.amount(*c).to_si() != start[i]);
                }
            }
            println!(
                "  started at the closed form, h = {h:.2} s: worst {:.1} eps, {ratio:.3} of \
                 (n + 2) eps, {moved} of {} amounts off their starting bits",
                worst / f64::EPSILON,
                4 * steps
            );
            assert!(ratio < 1.0, "h = {h}: the fixed point moved by {worst:e}");
        }

        // And the readers say the same closed form. `E = 60/(90 + 60) = 0.4`, and the liver's
        // venous blood is `(1 − E)` of the arterial.
        let b = body(60.0);
        assert_eq!(b.m.extraction_ratio(b.liver), Some(0.4));
        assert_eq!(b.m.extraction_ratio(b.muscle), Some(0.0));
        assert_eq!(b.m.extraction_ratio(b.blood), None);
        assert_eq!(b.m.partition(b.fat), 5.0);
        assert_eq!(b.m.partition(b.blood), 1.0);
        assert_eq!(b.m.blood_flow(b.muscle), l_h(50.0));
        assert_eq!(b.m.blood_flow(b.blood).to_si(), 0.0);

        // From empty.
        let mut b = body(60.0);
        b.m.infuse(b.blood, Mass::from_si(1e-3), Time::ZERO, hours(2000.0))
            .unwrap();
        let h = 60.0;
        let pools = [b.blood, b.liver, b.muscle, b.fat];
        let mut bus = Exchange::new();
        let mut n = 0usize;
        for until in [400.0, 600.0] {
            while (n as f64) * h < until * HOUR {
                b.m.step(Time::s(n as f64 * h), Time::s(h), &mut bus)
                    .unwrap();
                n += 1;
            }
            let rel: Vec<f64> = pools
                .iter()
                .zip(want)
                .map(|(c, x)| b.m.concentration(*c).to_si() / x - 1.0)
                .collect();
            println!("  from empty at {until} h: {rel:?}");
            for (r, c) in rel.iter().zip(pools) {
                assert!(r.abs() < 1e-12, "{} at {until} h: {r:e}", b.m.label(c));
            }
        }
        // The liver's extraction, measured off the state: arterial in, venous out.
        let e =
            1.0 - b.m.venous_concentration(b.liver).to_si() / b.m.concentration(b.blood).to_si();
        println!("  measured E = {e:.15}");
        assert!((e - 0.4).abs() < 1e-12, "E = {e}");
    }

    /// **PBPK closed form 2: `V_ss = V_b + Σ K_p V_t`, from the state and from the moments.**
    ///
    /// The liver's clearance is set to zero here, so only blood eliminates and `V_ss` is exact:
    /// `5 + 4·1.8 + 2·30 + 5·10 = 122.2 L`.
    ///
    /// **From the state**, at steady state under an infusion: `body_burden/C_b`. Each tissue holds
    /// `K_p C_b V_t`, so the ratio is `V_ss` exactly once the transient has gone.
    ///
    /// **From the moments**, after a bolus into the blood: `V_ss = CL·AUMC/AUC`, the mean
    /// residence time times the clearance. The stepped curve has these moments **exactly**, with
    /// one correction, and that is derived rather than measured. With `P = I + hM` and the left
    /// sums `AUC_d = h Σ C_n`, `AUMC_d = h Σ (nh) C_n`, the geometric series of matrices give
    ///
    /// ```text
    /// h Σ Pⁿ = h(I − P)⁻¹ = −M⁻¹                       (the continuous AUC, at any h)
    /// h² Σ n Pⁿ = h² P (I − P)⁻² = M⁻² + h M⁻¹          (the continuous AUMC, less h·AUC)
    /// ```
    ///
    /// so `AUMC_d/AUC_d = MRT − h` to the last digit, and `CL·(AUMC_d/AUC_d + h) = V_ss` at any
    /// step — not first order, exact. Summed to 1 200 h, where the truncated tail is below the
    /// AUC's own measured agreement with `D/CL`. The bound is `n·ε` for `n` steps of rounding and
    /// an `n`-term sum, the AUC test's bound; measured at 0.002 to 0.003 of it, 3.9e-14 and
    /// 9.8e-14 at `h` = 60 and 30 s, and the state's ratio at 1.2e-14.
    ///
    /// The `+ h` is not slack: without it the moment estimate is off by `h·CL/V_ss`, 8.2e-4 at
    /// a 60 s step and 4.1e-4 at 30 s — first order, and so exactly what a reader who had not
    /// derived the correction would have taken for discretisation error and loosened a tolerance
    /// to. It is asserted as the derived number.
    #[test]
    fn the_steady_state_volume_is_the_blood_plus_kp_times_each_tissue() {
        let v_ss = 5.0 + 4.0 * 1.8 + 2.0 * 30.0 + 5.0 * 10.0;
        let b = body(0.0);
        assert!(
            (b.m.volume_of_distribution().in_litres() / v_ss - 1.0).abs() < 4.0 * f64::EPSILON,
            "{} L",
            b.m.volume_of_distribution().in_litres()
        );

        // From the state.
        let mut b = body(0.0);
        b.m.infuse(b.blood, Mass::from_si(1e-3), Time::ZERO, hours(3000.0))
            .unwrap();
        let h = 60.0;
        let steps = (1500.0 * HOUR / h) as usize;
        run(&mut b.m, Time::ZERO, Time::s(h), steps);
        let ratio = b.m.body_burden().to_si() / b.m.concentration(b.blood).to_si();
        let rel = ratio / Volume::litres(v_ss).to_si() - 1.0;
        println!("  from the state at 1500 h: V_ss off by {rel:.2e}");
        assert!(rel.abs() < 1e-12, "{} L", ratio * 1e3);

        // From the moments.
        let cl = l_h(6.0).to_si();
        for h in [60.0, 30.0] {
            let mut b = body(0.0);
            let dose = 1e-4;
            b.m.bolus(b.blood, Mass::from_si(dose)).unwrap();
            let steps = (1200.0 * HOUR / h) as usize;
            let (mut auc, mut aumc) = (0.0, 0.0);
            let mut bus = Exchange::new();
            for n in 0..steps {
                let c = b.m.concentration(b.blood).to_si();
                auc += c * h;
                aumc += n as f64 * h * c * h;
                b.m.step(Time::s(n as f64 * h), Time::s(h), &mut bus)
                    .unwrap();
            }
            let auc_rel = auc * cl / dose - 1.0;
            let mrt = aumc / auc;
            let from_moments = cl * (mrt + h) / Volume::litres(v_ss).to_si() - 1.0;
            let uncorrected = cl * mrt / Volume::litres(v_ss).to_si() - 1.0;
            println!(
                "  h = {h} s: AUC·CL/D off by {auc_rel:.2e}, CL·(MRT_d + h)/V_ss by \
                 {from_moments:.2e} ({:.3} n eps), without the + h by {uncorrected:.2e}",
                from_moments.abs() / (steps as f64 * f64::EPSILON)
            );
            let bound = steps as f64 * f64::EPSILON;
            assert!(auc_rel.abs() < bound, "AUC = D/CL: {auc_rel:e}");
            assert!(
                from_moments.abs() < bound,
                "V_ss from the moments: {from_moments:e}"
            );
            // And the correction is the derived one: `−h·CL/V_ss`, to the same bound.
            let expected = -h * cl / Volume::litres(v_ss).to_si();
            assert!(
                (uncorrected - expected).abs() < bound,
                "{uncorrected:e} against {expected:e}"
            );
        }
    }

    /// Blood and one tissue: `V_b = 5 L` eliminating `CL = 10 L/h`, `V_t = 20 L` at `K_p = 3`,
    /// perfused at `Q = 30 L/h`. Micro constants `k10 = CL/V_b = 2/h`, `k12 = Q/V_b = 6/h` and
    /// `k21 = Q/(K_p V_t) = 0.5/h`. A `K_p` on the blood's side would give `k12 = 2/h` and
    /// `k21 = 1.5/h`, which is the point of picking `K_p ≠ 1`.
    fn blood_and_tissue() -> (CompartmentModel, Compartment, Compartment) {
        let mut m = CompartmentModel::new("subject");
        let blood = m.eliminating("blood", Volume::litres(5.0), l_h(10.0));
        let tissue = m.compartment("tissue", Volume::litres(20.0));
        m.perfuse(blood, tissue, l_h(30.0), 3.0).unwrap();
        (m, blood, tissue)
    }

    /// `[k10, k12, k21, α, β]` for [`blood_and_tissue`], per second, from the closed form.
    ///
    /// `β` as `k10 k21/α` rather than `(sum − disc)/2`: the two roots differ by a factor of 30,
    /// and the subtraction would lose five bits of the small one.
    fn blood_and_tissue_roots() -> [f64; 5] {
        let (k10, k12, k21) = (2.0 / HOUR, 6.0 / HOUR, 0.5 / HOUR);
        let sum = k10 + k12 + k21;
        let alpha = 0.5 * (sum + (sum * sum - 4.0 * k10 * k21).sqrt());
        [k10, k12, k21, alpha, k10 * k21 / alpha]
    }

    /// **PBPK closed form 3: blood and one tissue after a bolus — the bi-exponential, for the
    /// asymmetric link, and the exact solution of the recurrence that steps it.**
    ///
    /// # The derivation
    ///
    /// In amounts, with the bolus `D` into the blood,
    ///
    /// ```text
    /// d/dt [A_b]   [ −(Q + CL)/V_b    Q/(K_p V_t) ] [A_b]     [ −(k10 + k12)   k21 ]
    ///      [A_t] = [    Q/V_b        −Q/(K_p V_t) ] [A_t]  =  [     k12       −k21 ]
    /// ```
    ///
    /// which is the classical two-compartment matrix with `k21 = Q/(K_p V_t)`. So `α` and `β` are
    /// the roots of `λ² − (k10 + k12 + k21)λ + k10 k21`, and
    ///
    /// ```text
    /// A_b(t) = D [ (α − k21)/(α − β) e^{−αt} + (k21 − β)/(α − β) e^{−βt} ]
    /// A_t(t) = D k12/(α − β) [ e^{−βt} − e^{−αt} ]
    /// ```
    ///
    /// Here `α + β = 8.5/h` and `αβ = 1/h²`, so `α = 8.38068/h` and `β = 0.11932/h`.
    ///
    /// **The recurrence.** Euler steps `xₙ₊₁ = (I + hM) xₙ`. `I + hM` has the eigenvectors of `M`
    /// and the eigenvalues `1 − hα` and `1 − hβ`, and the bolus decomposes on those eigenvectors
    /// exactly as it does in the continuous solution, so the stepped trajectory is **the same
    /// bi-exponential with `e^{−λt}` replaced by `(1 − hλ)ⁿ`** — computed here by repeated
    /// multiplication, never by `powi`. That is the exact solution of the scheme the model steps,
    /// and the model must agree with it at every point to rounding: a wrong error *constant*
    /// still converges at first order, and this is what pins it.
    ///
    /// # The bound, and what it catches
    ///
    /// `(3n + 40)ε` times the sum of the magnitudes of the four terms, the envelope a cancelling
    /// difference is bounded by: the model rounds about twice per step on what it carries, the
    /// reference's powers once a step, and `α`, `β` and the coefficients a handful of ulps,
    /// amplified by `n` through the powers — the oral recurrence's reasoning, without its depot.
    /// Measured below at 72 to 9 s. A `K_p` on the wrong side changes `k12` and `k21` by factors
    /// of three and is off at the first step by `≈ hQ·C_b`; a perfusion with no return flow
    /// has no `k21` at all.
    ///
    /// Then convergence on the continuous bi-exponential at first order, each compartment
    /// separately, at `h` = 18, 9 and 4.5 s, with the tri-exponential test's halving band.
    #[test]
    fn blood_and_one_tissue_follow_the_bi_exponential_and_its_recurrence() {
        let [k10, k12, k21, alpha, beta] = blood_and_tissue_roots();
        println!(
            "  alpha = {:.5}/h, beta = {:.5}/h",
            alpha * HOUR,
            beta * HOUR
        );
        assert!((alpha + beta - (k10 + k12 + k21)).abs() < 4.0 * f64::EPSILON * alpha);
        // The model's own rate constants are the closed form's, in both directions.
        let (m, blood, tissue) = blood_and_tissue();
        for (got, want, which) in [
            (m.transfer_rate(blood, tissue).to_si(), k12, "k12"),
            (m.transfer_rate(tissue, blood).to_si(), k21, "k21"),
            (m.elimination_rate(blood).to_si(), k10, "k10"),
        ] {
            assert!(
                (got / want - 1.0).abs() < 2.0 * f64::EPSILON,
                "{which}: {} /h against {} /h",
                got * HOUR,
                want * HOUR
            );
        }

        let dose = 1e-4;
        let coeff_b = [
            (alpha - k21) / (alpha - beta),
            (k21 - beta) / (alpha - beta),
        ];
        let coeff_t = [-k12 / (alpha - beta), k12 / (alpha - beta)];
        // Before any code under test: A_b(0) = D, A_t(0) = 0, A_t'(0) = k12 D.
        assert!((coeff_b[0] + coeff_b[1] - 1.0).abs() < 2.0 * f64::EPSILON);
        assert!((coeff_t[0] + coeff_t[1]).abs() < 1e-30);
        assert!(((-alpha * coeff_t[0] - beta * coeff_t[1]) / k12 - 1.0).abs() < 4.0 * f64::EPSILON);

        // The recurrence's exact solution, at four steps.
        for h in [72.0, 36.0, 18.0, 9.0] {
            let (mut m, blood, tissue) = blood_and_tissue();
            m.bolus(blood, Mass::from_si(dose)).unwrap();
            let (mu_a, mu_b) = (1.0 - h * alpha, 1.0 - h * beta);
            let (mut pa, mut pb) = (1.0f64, 1.0f64); // (1 − hα)ⁿ, (1 − hβ)ⁿ
            let mut bus = Exchange::new();
            let mut worst = 0.0f64;
            let steps = (12.0 * HOUR / h) as usize;
            for n in 0..steps {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                pa *= mu_a;
                pb *= mu_b;
                let n = (n + 1) as f64;
                let want_b = dose * (coeff_b[0] * pa + coeff_b[1] * pb);
                let want_t = dose * (coeff_t[0] * pa + coeff_t[1] * pb);
                let envelope = dose
                    * (coeff_b[0].abs() * pa.abs()
                        + coeff_b[1].abs() * pb
                        + coeff_t[0].abs() * pa.abs()
                        + coeff_t[1].abs() * pb);
                let bound = (3.0 * n + 40.0) * f64::EPSILON * envelope;
                for (got, want, which) in [
                    (m.amount(blood).to_si(), want_b, "blood"),
                    (m.amount(tissue).to_si(), want_t, "tissue"),
                ] {
                    let ratio = (got - want).abs() / bound;
                    assert!(
                        ratio < 1.0,
                        "h = {h} s, step {n}, {which}: {got:e} against the recurrence's {want:e}"
                    );
                    worst = worst.max(ratio);
                }
            }
            println!("  h = {h:>4} s: worst {worst:.3} of the rounding bound");
        }

        // And the continuous closed form, at first order.
        let exact_b =
            |t: f64| dose * (coeff_b[0] * (-alpha * t).exp() + coeff_b[1] * (-beta * t).exp());
        let exact_t =
            |t: f64| dose * (coeff_t[0] * (-alpha * t).exp() + coeff_t[1] * (-beta * t).exp());
        let (mut errors_b, mut errors_t) = (Vec::new(), Vec::new());
        for h in [18.0, 9.0, 4.5] {
            let (mut m, blood, tissue) = blood_and_tissue();
            m.bolus(blood, Mass::from_si(dose)).unwrap();
            let mut bus = Exchange::new();
            let (mut wb, mut wt) = (0.0f64, 0.0f64);
            for n in 0..(12.0 * HOUR / h) as usize {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                let t = (n + 1) as f64 * h;
                wb = wb.max((m.amount(blood).to_si() - exact_b(t)).abs());
                wt = wt.max((m.amount(tissue).to_si() - exact_t(t)).abs());
            }
            errors_b.push(wb / dose);
            errors_t.push(wt / dose);
        }
        for (which, errors) in [("blood", &errors_b), ("tissue", &errors_t)] {
            let measured = orders(errors);
            println!("  bi-exponential, {which}: errors {errors:?}, orders {measured:.4?}");
            assert!(errors[0] < 0.01, "{which} near the closed form at all");
            let mut h = 18.0;
            for p in &measured {
                let band = 0.035 * h / 18.0;
                h /= 2.0;
                assert!(
                    (p - 1.0).abs() < band,
                    "{which}: first order, measured {p:.4}"
                );
            }
            assert!(
                (measured[1] - 1.0).abs() < (measured[0] - 1.0).abs(),
                "{which}: the order should approach 1 as h falls: {measured:.4?}"
            );
            assert!(*errors.last().unwrap() > 1e-6, "{which} is above the floor");
        }
    }

    /// **PBPK closed form 4: as `Q → ∞` the model is one compartment of `V_b + Σ K_p V_t`, and the
    /// deviation falls as `1/Q`.**
    ///
    /// Blood (5 L, `CL` 6 L/h), a muscle (30 L, `K_p` 2) and a fat (10 L, `K_p` 5), so
    /// `V_ss = 115 L` and the limit's rate is `k = CL/V_ss = 0.0522/h`. Both blood flows are
    /// scaled by `s` = 1, 2, 4, 8, 16 from 60 and 50 L/h, where each tissue already turns over at
    /// 1/h — twenty times `k`.
    ///
    /// **What converges, and why it is first order in `1/Q`.** The distribution modes are as fast
    /// as `Q` and are gone; what is left is the slow mode, `a(Q)·e^{−β(Q)t}`, and both its
    /// amplitude and its rate are analytic in `1/Q` about the one-compartment `D/V_ss` and `k`.
    /// So the deviation from one compartment is `O(1/Q)`, and its order is measured the way a step
    /// is: `log₂` of the ratio between successive doublings of `Q`.
    ///
    /// **Stepped at one `h` for every `Q`, against the one-compartment recurrence at that `h`**,
    /// `C_n = (D/V_ss)(1 − hk)ⁿ` by repeated multiplication. That is what keeps the measurement
    /// about `Q` and not about `h`: the eigenvalues of `I + hM` are `1 + hλ` exactly, so the
    /// stepped slow mode is `(1 − hβ(Q))ⁿ` and differs from `(1 − hk)ⁿ` only by `β(Q) − k`. Compared
    /// with the *continuous* exponential instead, the Euler error `t k² h/2` would not fall with
    /// `Q` and would put a floor under the order. `h = 5 s` keeps the turnover number under
    /// 0.5 at the largest flow.
    ///
    /// **Sampled at one time, 48 h, and not as a maximum over a window**, relative to `D/V_ss`.
    /// The first version took the worst point from 24 to 72 h, and its orders were 0.9997,
    /// 1.0005, 1.0004, 1.0003 — first order, and with departures that would not halve, so no
    /// band could be earned from them. The continuous solution, by eigendecomposition outside
    /// this crate, said why: the worst point sat at 24 h, where the `1/Q²` coefficient of the
    /// deviation nearly vanishes, and what was left was the `1/Q³` share and the `s = 1`
    /// distribution mode. At 48 h the same eigendecomposition gives orders 0.8949, 0.9488,
    /// 0.9747, 0.9874, departures halving from 0.105, and that is what is asserted: a band of
    /// `0.2·(Q₀/Q)` on each, halving with the `1/Q²` share it bounds, and the last departure
    /// under a quarter of the first.
    #[test]
    fn fast_perfusion_collapses_to_one_compartment() {
        let dose = 1e-4;
        let v_ss = Volume::litres(5.0 + 2.0 * 30.0 + 5.0 * 10.0).to_si();
        let cl = l_h(6.0).to_si();
        let h = 5.0;
        let at = (48.0 * HOUR / h) as usize;
        let mut deviations = Vec::new();
        for s in [1.0, 2.0, 4.0, 8.0, 16.0] {
            let mut m = CompartmentModel::new("subject");
            let blood = m.eliminating("blood", Volume::litres(5.0), l_h(6.0));
            let muscle = m.compartment("muscle", Volume::litres(30.0));
            let fat = m.compartment("fat", Volume::litres(10.0));
            m.perfuse(blood, muscle, l_h(60.0 * s), 2.0).unwrap();
            m.perfuse(blood, fat, l_h(50.0 * s), 5.0).unwrap();
            assert!(m.turnover_number(Time::s(h)) < 0.5);
            assert!((m.volume_of_distribution().to_si() / v_ss - 1.0).abs() < 4.0 * f64::EPSILON);
            m.bolus(blood, Mass::from_si(dose)).unwrap();
            let r = 1.0 - h * cl / v_ss;
            let mut one = dose / v_ss; // (D/V_ss) rⁿ
            let mut bus = Exchange::new();
            for n in 0..at {
                m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
                one *= r;
            }
            deviations.push((m.concentration(blood).to_si() - one).abs() / (dose / v_ss));
        }
        let measured = orders(&deviations);
        println!("  fast perfusion: deviations {deviations:?}, orders {measured:.4?}");
        for (i, p) in measured.iter().enumerate() {
            let band = 0.2 / (1u32 << i) as f64;
            assert!(
                (p - 1.0).abs() < band,
                "the deviation should fall as 1/Q: order {p:.4} at s = {}",
                1u32 << i
            );
        }
        assert!(
            (measured[measured.len() - 1] - 1.0).abs() < 0.25 * (measured[0] - 1.0).abs(),
            "the order should approach 1 as Q grows: {measured:.4?}"
        );
        assert!(*deviations.last().unwrap() > 1e-6, "above the floor");
    }

    /// **`K_p = 1` is the symmetric link, to the bit.**
    ///
    /// Not a bound: `x/1.0` is `x` in IEEE arithmetic, and the perfusion's flux, the venous
    /// concentration elimination acts on, the turnover rate and the rate constants are written
    /// so that with `K_p = 1` they are the link's operations in the link's order. So two models
    /// — one linked, one perfused at `K_p = 1`, each with a second, symmetric link beside it, a
    /// clearance from the blood *and* the tissue, a bolus, an infusion and a depot into the
    /// tissue — must agree in every amount, the books, the ledger and the readings at every step,
    /// and in the step limit and every rate reported, with `assert_eq`.
    #[test]
    fn a_partition_of_one_is_the_symmetric_link_to_the_bit() {
        let build = |perfused: bool| {
            let mut m = CompartmentModel::new("subject");
            let blood = m.eliminating("blood", Volume::litres(5.0), l_h(7.0));
            let tissue = m.eliminating("tissue", Volume::litres(13.0), l_h(3.3));
            let other = m.compartment("other", Volume::litres(21.0));
            if perfused {
                m.perfuse(blood, tissue, l_h(37.0), 1.0).unwrap();
            } else {
                m.link(blood, tissue, l_h(37.0)).unwrap();
            }
            m.link(other, blood, l_h(11.0)).unwrap();
            m.bolus(blood, Mass::from_si(1e-4)).unwrap();
            m.infuse(other, Mass::from_si(2e-4), hours(0.3), hours(2.7))
                .unwrap();
            let gut = m
                .depot("gut", tissue, per_hour(1.1), 0.7, hours(0.2))
                .unwrap();
            m.oral(gut, Mass::from_si(5e-5), hours(0.05)).unwrap();
            (m, blood, tissue, other)
        };
        let (mut a, ab, at, ao) = build(false);
        let (mut b, bb, bt, bo) = build(true);
        assert_eq!(a.max_stable_dt(Time::ZERO), b.max_stable_dt(Time::ZERO));
        assert_eq!(a.volume_of_distribution(), b.volume_of_distribution());
        assert_eq!(a.total_clearance(), b.total_clearance());
        for (x, y) in [(ab, at), (at, ab), (ab, ao), (at, ao)] {
            assert_eq!(a.transfer_rate(x, y), b.transfer_rate(x, y));
        }
        assert_eq!(a.elimination_rate(at), b.elimination_rate(bt));
        let h = 30.0;
        let mut bus = Exchange::new();
        for n in 0..(8.0 * HOUR / h) as usize {
            let t = Time::s(n as f64 * h);
            a.step(t, Time::s(h), &mut bus).unwrap();
            b.step(t, Time::s(h), &mut bus).unwrap();
            for (x, y) in [(ab, bb), (at, bt), (ao, bo)] {
                assert_eq!(a.amount(x), b.amount(y), "step {n}, {}", a.label(x));
                assert_eq!(a.concentration(x), b.venous_concentration(y));
            }
            assert_eq!(a.transfer(ab, at), b.transfer(bb, bt), "step {n}");
            assert_eq!(a.transfer(at, ab), b.transfer(bt, bb), "step {n}");
            assert_eq!(a.eliminated_mass(), b.eliminated_mass());
            assert_eq!(a.ledger(), b.ledger());
        }
        assert_eq!(a.readings(), b.readings());
        assert!(a.eliminated_mass().to_si() > 1e-4, "not vacuous");
    }

    /// **The well-stirred liver: first-pass extraction, `F_H = 1 − E`, exactly, at any step.**
    ///
    /// A dose that enters the body *through* the liver — a depot into it, as an oral dose
    /// absorbed into the portal blood is — reaches the blood less what the liver takes on the
    /// way. With a renal clearance on the blood as well, so the answer is not trivially the
    /// organ's own:
    ///
    /// ```text
    /// AUC_b(dose into the liver) / AUC_b(dose into the blood) = Q_h/(Q_h + CL) = 1 − E
    /// ```
    ///
    /// whatever `CL_r` is. The derivation is the steady state's: `h Σ xₙ = −M⁻¹U` for a total
    /// input `U` (the AUC test's identity, for any input schedule), and `−M⁻¹e_l` is the steady
    /// state under a unit infusion into the liver, whose blood balance is
    /// `Q_h(C_l/K − C_b) = CL_r C_b` and liver balance `R = CL·C_l/K + Q_h(C_l/K − C_b)`; so
    /// `C_b = R/(CL + CL_r(Q_h + CL)/Q_h)`, against `R/(CL_r + Q_h E)` into the blood, and the
    /// ratio is `Q_h/(Q_h + CL)`.
    ///
    /// So it is **exact, not first order** — the stepped AUC is the continuous one — and it is
    /// asserted at the AUC test's `n·ε`. Both doses through a depot with the same `kₐ`, so the
    /// timing of the input is the same and only where it lands differs. `E = 0.4` here. Taking
    /// elimination off `C_t` instead of `C_t/K_p` makes the extraction `K CL/(Q + K CL)`, 0.73,
    /// and `F_H` 0.27 rather than 0.6. And the blood's AUC from a dose into it is `D/(CL_r + Q_h E)`.
    #[test]
    fn the_liver_extracts_what_the_well_stirred_model_says_on_the_first_pass() {
        let (q_h, cl, cl_r) = (l_h(90.0).to_si(), l_h(60.0).to_si(), l_h(6.0).to_si());
        let e = cl / (q_h + cl);
        let auc = |into_liver: bool| {
            let mut b = body(60.0);
            let target = if into_liver { b.liver } else { b.blood };
            let gut =
                b.m.depot("gut", target, per_hour(2.0), 1.0, Time::ZERO)
                    .unwrap();
            b.m.oral(gut, Mass::from_si(ORAL), Time::ZERO).unwrap();
            let h = 60.0;
            let steps = (1500.0 * HOUR / h) as usize;
            let c = trace(&mut b.m, b.blood, h, steps);
            (c[..steps].iter().map(|x| x * h).sum::<f64>(), steps)
        };
        let ((via_liver, n), (direct, _)) = (auc(true), auc(false));
        let bound = n as f64 * f64::EPSILON;
        let f_h = via_liver / direct;
        let systemic = direct * (cl_r + q_h * e) / ORAL;
        println!(
            "  F_H = {f_h:.15} against 1 - E = {}, off by {:.2e}; AUC·CL/D off by {:.2e}",
            1.0 - e,
            f_h / (1.0 - e) - 1.0,
            systemic - 1.0
        );
        assert!((f_h / (1.0 - e) - 1.0).abs() < bound, "F_H = {f_h}");
        assert!((systemic - 1.0).abs() < bound, "the systemic clearance");
    }

    /// **The step limit sees the tissue's own rate, `Q/(K_p V_t)`, and bounds the true one.**
    ///
    /// A small tissue that partitions weakly: blood 5 L with `Q` 30 L/h, a tissue of 0.5 L at
    /// `K_p = 0.25`. The blood turns over at `(30 + 10)/5 = 8/h`; the tissue at
    /// `30/(0.25·0.5) = 240/h`, thirty times faster, and that is the limit — 15 s. A limit taken
    /// from the blood alone would be 450 s, at which the tissue's diagonal factor `1 − h r` is −29.
    ///
    /// Against the closed form, as for the classical pair: the two eigenvalues are the roots of
    /// `λ² − (k10 + k12 + k21)λ + k10 k21`, and the limit must not exceed `2/α` and is within a
    /// factor of two of it. The step past it is refused naming the tissue, and at the limit no
    /// pool goes negative.
    #[test]
    fn the_step_limit_sees_the_tissue() {
        let mut m = CompartmentModel::new("subject");
        let blood = m.eliminating("blood", Volume::litres(5.0), l_h(10.0));
        let tissue = m.compartment("tissue", Volume::litres(0.5));
        m.perfuse(blood, tissue, l_h(30.0), 0.25).unwrap();
        let limit = m.max_stable_dt(Time::ZERO).to_si();
        let want = 1.0 / (l_h(30.0).to_si() / (0.25 * Volume::litres(0.5).to_si()));
        assert!(
            (limit / want - 1.0).abs() < 4.0 * f64::EPSILON,
            "{limit} s against {want} s"
        );
        assert!((limit - 15.0).abs() < 1e-9, "{limit} s");

        let (k10, k12, k21) = (2.0 / HOUR, 6.0 / HOUR, 240.0 / HOUR);
        let sum = k10 + k12 + k21;
        let alpha = 0.5 * (sum + (sum * sum - 4.0 * k10 * k21).sqrt());
        println!("  limit {limit} s, 2/alpha = {} s", 2.0 / alpha);
        assert!(limit <= 2.0 / alpha, "inside the true limit");
        assert!(limit > 1.0 / alpha, "and within a factor of two of it");

        m.bolus(blood, Mass::from_si(1e-4)).unwrap();
        let mut bus = Exchange::new();
        let err = m
            .step(Time::ZERO, Time::s(limit * 1.05), &mut bus)
            .expect_err("past the tissue's limit");
        assert_eq!(err.quantity, "compartment turnover number");
        assert!(err.site.contains("tissue"), "{}", err.site);
        for n in 0..2000 {
            m.step(Time::s(n as f64 * limit), Time::s(limit), &mut bus)
                .unwrap();
            for (c, label) in m.handles() {
                assert!(m.amount(c).to_si() >= 0.0, "{label} negative at step {n}");
            }
        }
        assert!(m.eliminated_mass().to_si() > 0.5e-4);
    }

    /// Run a model started at a closed-form steady state and return the worst drift of any
    /// listed compartment's concentration, as a fraction of `(n + 2)ε` after `n` steps — the
    /// bound [`under_an_infusion_every_tissue_settles_at_kp_times_the_blood`] derives: the fixed
    /// point is the continuous one, so only rounding moves it, and `2ε` is writing `C·V` in and
    /// reading `A/V` back.
    fn drift_from(
        m: &mut CompartmentModel,
        pools: &[(Compartment, f64)],
        h: f64,
        steps: usize,
    ) -> f64 {
        let mut bus = Exchange::new();
        let mut ratio = 0.0f64;
        for n in 0..steps {
            m.step(Time::s(n as f64 * h), Time::s(h), &mut bus).unwrap();
            for (c, want) in pools {
                let rel = (m.concentration(*c).to_si() / want - 1.0).abs();
                ratio = ratio.max(rel / ((n as f64 + 2.0) * f64::EPSILON));
            }
        }
        ratio
    }

    /// **A perfused tissue with a clearance *and* a link: it sets the step limit, and holds its
    /// steady state** — the branches every other perfused tissue here leaves unreached.
    ///
    /// Blood 5 L clearing 10 L/h; a tissue of 0.5 L at `K_p = 0.25` perfused at 30 L/h, clearing
    /// 15 L/h on `C_t/K_p`, and linked symmetrically at 6 L/h to a third pool of 2 L that only
    /// distributes. The tissue turns over at `6/0.5 + (30 + 15)/(0.25·0.5) = 12 + 360 = 372/h`,
    /// so the limit is `3600/372 = 9.677 s` — with the clearance left out of it 252/h and
    /// 14.3 s, with the link left out 360/h and 10.0 s, and with the perfusion left out 132/h.
    ///
    /// At steady state under an infusion `R` into the blood, `E = 15/45 = 1/3`,
    /// `C_b = R/(10 + 30·E) = R/(20 L/h)`, the tissue's venous `C_t/K_p = (1 − E) C_b`, so
    /// `C_t = C_b/6`; and the third pool, joined by a **symmetric** link, sits at `C_t` itself —
    /// not at `C_t/K_p`, which is what a link that read the venous concentration would give.
    /// Started there, it must stay to `(n + 2)ε`, as the body's steady state does.
    ///
    /// The readers too, at that state: `k = CL/(K_p V) = 120/h` out of the tissue,
    /// `Q/(K_p V) = 240/h` back to the blood and `Q/V_b = 6/h` into it, and a perfusion carrying
    /// `Q(C_b − C_t/K_p) = Q E C_b` — positive into the tissue, the negation the other way.
    #[test]
    fn a_tissue_with_a_clearance_and_a_link_sets_the_limit_and_holds_its_steady_state() {
        let build = || {
            let mut m = CompartmentModel::new("subject");
            let blood = m.eliminating("blood", Volume::litres(5.0), l_h(10.0));
            let tissue = m.eliminating("tissue", Volume::litres(0.5), l_h(15.0));
            let other = m.compartment("other", Volume::litres(2.0));
            m.perfuse(blood, tissue, l_h(30.0), 0.25).unwrap();
            m.link(tissue, other, l_h(6.0)).unwrap();
            (m, blood, tissue, other)
        };
        let (mut m, blood, tissue, other) = build();
        let limit = m.max_stable_dt(Time::ZERO).to_si();
        let want = 3600.0 / 372.0;
        println!("  limit {limit} s against {want} s");
        assert!((limit / want - 1.0).abs() < 4.0 * f64::EPSILON, "{limit} s");

        let rate = 1e-3 / (1000.0 * HOUR);
        let c_b = rate / l_h(20.0).to_si();
        let c_t = c_b / 6.0;
        let pools = [(blood, c_b), (tissue, c_t), (other, c_t)];
        for (c, x) in pools {
            m.bolus(c, Mass::from_si(x * m.volume(c).to_si())).unwrap();
        }
        m.infuse(blood, Mass::from_si(1e-3), Time::ZERO, hours(1000.0))
            .unwrap();
        for (got, want, which) in [
            (
                m.elimination_rate(tissue).to_si(),
                120.0 / HOUR,
                "k out of the tissue",
            ),
            (
                m.transfer_rate(tissue, blood).to_si(),
                240.0 / HOUR,
                "tissue to blood",
            ),
            (
                m.transfer_rate(blood, tissue).to_si(),
                6.0 / HOUR,
                "blood to tissue",
            ),
            (
                m.transfer_rate(tissue, other).to_si(),
                12.0 / HOUR,
                "tissue to other",
            ),
            (
                m.transfer(blood, tissue).to_si(),
                l_h(30.0).to_si() * c_b / 3.0,
                "the perfusion's flux, Q E C_b",
            ),
            (
                m.transfer(tissue, blood).to_si(),
                -l_h(30.0).to_si() * c_b / 3.0,
                "and its negation",
            ),
        ] {
            assert!(
                (got / want - 1.0).abs() < 8.0 * f64::EPSILON,
                "{which}: {got:e} against {want:e}"
            );
        }
        assert!(
            m.transfer(tissue, other).to_si().abs() < 1e-30,
            "no net flux on the link"
        );
        let ratio = drift_from(&mut m, &pools, 0.98 * limit, 4000);
        println!("  started at the closed form: {ratio:.3} of (n + 2) eps");
        assert!(
            ratio < 1.0,
            "the steady state moved: {ratio} of (n + 2) eps"
        );

        // Refused just past the limit, naming the tissue; never negative at it.
        let (mut m, blood, _, _) = build();
        m.bolus(blood, Mass::from_si(1e-4)).unwrap();
        let mut bus = Exchange::new();
        let err = m
            .step(Time::ZERO, Time::s(1.05 * limit), &mut bus)
            .expect_err("past the tissue's limit");
        assert!(err.site.contains("tissue"), "{}", err.site);
        for n in 0..3000 {
            m.step(Time::s(n as f64 * limit), Time::s(limit), &mut bus)
                .unwrap();
            for (c, label) in m.handles() {
                assert!(m.amount(c).to_si() >= 0.0, "{label} negative at step {n}");
            }
        }
    }

    /// **A link and a perfusion on the same pair add, in the readers and in the step.**
    ///
    /// Blood 5 L and a tissue of 4 L at `K_p = 3`, perfused at 20 L/h and also linked at 2 L/h,
    /// each holding some drug. `transfer_rate` must be `(20 + 2)/5` out of the blood and
    /// `2/4 + 20/(3·4)` out of the tissue, `transfer` `2(C_b − C_t) + 20(C_b − C_t/3)`, and one
    /// step must move exactly `transfer·h` into the tissue — the reader and the step agreeing, to
    /// the rounding of one product.
    #[test]
    fn a_link_beside_a_perfusion_adds_to_it() {
        let mut m = CompartmentModel::new("subject");
        let blood = m.compartment("blood", Volume::litres(5.0));
        let tissue = m.compartment("tissue", Volume::litres(4.0));
        m.perfuse(blood, tissue, l_h(20.0), 3.0).unwrap();
        m.link(tissue, blood, l_h(2.0)).unwrap();
        m.bolus(blood, Mass::from_si(1e-4)).unwrap();
        m.bolus(tissue, Mass::from_si(3e-5)).unwrap();
        let (c_b, c_t) = (1e-4 / 5e-3, 3e-5 / 4e-3);
        let (q_p, q_l) = (l_h(20.0).to_si(), l_h(2.0).to_si());
        let flux = q_l * (c_b - c_t) + q_p * (c_b - c_t / 3.0);
        for (got, want, which) in [
            (
                m.transfer_rate(blood, tissue).to_si(),
                22.0 / 5.0 / HOUR,
                "blood to tissue",
            ),
            (
                m.transfer_rate(tissue, blood).to_si(),
                (2.0 / 4.0 + 20.0 / 12.0) / HOUR,
                "tissue to blood",
            ),
            (m.transfer(blood, tissue).to_si(), flux, "the flux"),
            (m.transfer(tissue, blood).to_si(), -flux, "its negation"),
        ] {
            assert!(
                (got / want - 1.0).abs() < 8.0 * f64::EPSILON,
                "{which}: {got:e} against {want:e}"
            );
        }
        let h = 10.0;
        let before = m.amount(tissue).to_si();
        run(&mut m, Time::ZERO, Time::s(h), 1);
        let moved = m.amount(tissue).to_si() - before;
        assert!(
            (moved / (flux * h) - 1.0).abs() < 1e-12,
            "the step moved {moved:e}, the reader says {:e}",
            flux * h
        );
    }

    /// **Two blood pools, each with its own tissue, are two models.**
    ///
    /// Nothing else here has more than one blood pool, so a step that sent every perfusion's
    /// drug to the first pool's blood would pass all of it. Pool A: 5 L clearing 6 L/h, a tissue
    /// of 10 L at `K_p = 2` perfused at 40 L/h; pool B: 3 L clearing 3 L/h, a tissue of 7 L at
    /// `K_p = 7` perfused at 15 L/h; infused at different rates. Each blood settles at `R/CL` of
    /// its own and each tissue at its own `K_p` times it — started there and held to
    /// `(n + 2)ε`. That alone passed a step that sent every perfusion to the first pool's blood,
    /// because a perfusion carries nothing at steady state; the bolus half below is what sees it.
    #[test]
    fn two_blood_pools_are_two_models() {
        let mut m = CompartmentModel::new("subject");
        let a = m.eliminating("blood A", Volume::litres(5.0), l_h(6.0));
        let ta = m.compartment("tissue A", Volume::litres(10.0));
        let b = m.eliminating("blood B", Volume::litres(3.0), l_h(3.0));
        let tb = m.compartment("tissue B", Volume::litres(7.0));
        m.perfuse(a, ta, l_h(40.0), 2.0).unwrap();
        m.perfuse(b, tb, l_h(15.0), 7.0).unwrap();
        let (ra, rb) = (1e-3 / (1000.0 * HOUR), 2e-4 / (1000.0 * HOUR));
        let (ca, cb) = (ra / l_h(6.0).to_si(), rb / l_h(3.0).to_si());
        let pools = [(a, ca), (ta, 2.0 * ca), (b, cb), (tb, 7.0 * cb)];
        for (c, x) in pools {
            m.bolus(c, Mass::from_si(x * m.volume(c).to_si())).unwrap();
        }
        m.infuse(a, Mass::from_si(1e-3), Time::ZERO, hours(1000.0))
            .unwrap();
        m.infuse(b, Mass::from_si(2e-4), Time::ZERO, hours(1000.0))
            .unwrap();
        let h = 0.98 * m.max_stable_dt(Time::ZERO).to_si();
        let ratio = drift_from(&mut m, &pools, h, 4000);
        println!("  two pools, started at the closed form: {ratio:.3} of (n + 2) eps");
        assert!(ratio < 1.0, "a pool moved: {ratio} of (n + 2) eps");

        // **And away from it.** At steady state every perfusion's flux is zero, so a step that
        // sent one pool's perfusion to the other's blood moves nothing and the check above
        // passes it — measured: it did. So a bolus into each pool in turn: the other must stay
        // empty **to the bit**, and the dosed one must be a model of that pool alone, to the bit,
        // because its arithmetic is the same operations on the same numbers.
        let pool = |which: u8, with_both: bool| {
            let mut m = CompartmentModel::new("subject");
            let mut handles = Vec::new();
            for (k, (vb, cl, vt, q, kp)) in
                [(5.0, 6.0, 10.0, 40.0, 2.0), (3.0, 3.0, 7.0, 15.0, 7.0)]
                    .into_iter()
                    .enumerate()
            {
                if !with_both && k as u8 != which {
                    continue;
                }
                let b = m.eliminating(format!("blood {k}"), Volume::litres(vb), l_h(cl));
                let t = m.compartment(format!("tissue {k}"), Volume::litres(vt));
                m.perfuse(b, t, l_h(q), kp).unwrap();
                handles.push((k as u8, b, t));
            }
            (m, handles)
        };
        for dosed in [0u8, 1] {
            let (mut both, hb) = pool(dosed, true);
            let (mut alone, ha) = pool(dosed, false);
            let (_, db, dt) = *hb.iter().find(|x| x.0 == dosed).unwrap();
            let (_, ob, ot) = *hb.iter().find(|x| x.0 != dosed).unwrap();
            let (_, ab, at) = ha[0];
            both.bolus(db, Mass::from_si(1e-4)).unwrap();
            alone.bolus(ab, Mass::from_si(1e-4)).unwrap();
            let h = 30.0;
            let mut bus = Exchange::new();
            for n in 0..(24.0 * HOUR / h) as usize {
                let t = Time::s(n as f64 * h);
                both.step(t, Time::s(h), &mut bus).unwrap();
                alone.step(t, Time::s(h), &mut bus).unwrap();
                assert_eq!(
                    (both.amount(ob).to_si(), both.amount(ot).to_si()),
                    (0.0, 0.0)
                );
                assert_eq!(both.amount(db), alone.amount(ab), "pool {dosed}, step {n}");
                assert_eq!(both.amount(dt), alone.amount(at), "pool {dosed}, step {n}");
            }
            assert!(both.amount(dt).to_si() > 0.0, "not vacuous");
        }
    }

    /// **Mass balance with perfusion, hepatic and renal elimination, a depot, a bolus and an
    /// infusion, at every step, to machine precision** — by the accessors and by
    /// [`Domain::ledger`].
    ///
    /// The ledger is blind to perfusions by construction, as it is to links: `+q` and `−q` are in
    /// one sum. What this check sees is a perfusion that moves drug out of one pool and not into
    /// the other, and an elimination that removes drug without booking it. A depot into the liver
    /// with a lag, a symmetric link to a compartment that is not a tissue, a bolus and an
    /// infusion, so every entry of the books moves. The bound is the existing `1e-12`, for the
    /// existing reason.
    #[test]
    fn mass_balance_holds_with_perfusion_and_hepatic_elimination() {
        let mut b = body(60.0);
        let gut_wall = b.m.compartment("gut wall", Volume::litres(2.0));
        b.m.link(gut_wall, b.blood, l_h(4.0)).unwrap();
        let gut =
            b.m.depot("gut", b.liver, per_hour(1.4), 0.55, hours(0.3))
                .unwrap();
        let (bolus, infused, oral) = (1e-4, 2e-4, 3e-4);
        b.m.bolus(b.blood, Mass::from_si(bolus)).unwrap();
        b.m.infuse(b.muscle, Mass::from_si(infused), hours(0.41), hours(3.3))
            .unwrap();
        b.m.oral(gut, Mass::from_si(oral), hours(0.77)).unwrap();
        let committed = bolus + infused + oral;
        let h = 60.0;
        let steps = (48.0 * HOUR / h) as usize;
        let mut bus = Exchange::new();
        let (mut worst, mut worst_ledger) = (0.0f64, 0.0f64);
        for n in 0..steps {
            b.m.step(Time::s(n as f64 * h), Time::s(h), &mut bus)
                .unwrap();
            let books = b.m.body_burden().to_si()
                + b.m.eliminated_mass().to_si()
                + b.m.depot_mass().to_si()
                + b.m.presystemic_loss().to_si()
                + b.m.pending_dose().to_si();
            worst = worst.max((books - committed).abs() / committed);
            let ledger = b.m.ledger().get(quantity::MASS).unwrap();
            worst_ledger = worst_ledger.max((ledger - committed).abs() / committed);
        }
        println!("  books {worst:.2e}, ledger {worst_ledger:.2e} over {steps} steps");
        assert!(worst < 1e-12, "the books drifted by {worst:e}");
        assert!(
            worst_ledger < 1e-12,
            "the ledger drifted by {worst_ledger:e}"
        );
        assert!(
            b.m.eliminated_mass().to_si() > 0.3 * committed,
            "not vacuous"
        );
        assert!(b.m.amount(b.fat).to_si() > 0.0 && b.m.amount(gut_wall).to_si() > 0.0);
    }

    /// **The same bits however a perfused run is cut**, and after a `restore`.
    ///
    /// A perfusion holds no state of its own — its `Q` and `K_p` are configuration — so a
    /// checkpoint needs nothing new, and this is the check that that is true: seventy-two
    /// scheduler windows of 100 s against thirty-six of 200, which the multirate schedule cuts
    /// into the same 100 s substeps under a limit of 115 s (asserted, and the scheduler audits the
    /// ledger at `1e-12` every window), and a run rewound part-way through and run again, which
    /// must also be the scheduler's run.
    #[test]
    fn a_perfused_run_is_the_same_bits_however_it_is_cut() {
        let build = || {
            let mut b = body(60.0);
            b.m.bolus(b.blood, Mass::from_si(5e-5)).unwrap();
            let gut =
                b.m.depot("gut", b.liver, per_hour(1.7), 0.75, Time::s(500.0))
                    .unwrap();
            b.m.oral(gut, Mass::from_si(ORAL), Time::s(130.0)).unwrap();
            b.m.infuse(b.fat, Mass::from_si(4e-5), Time::s(700.0), Time::s(5100.0))
                .unwrap();
            b.m
        };
        let state = |m: &CompartmentModel| -> Vec<f64> {
            let mut v: Vec<f64> = m.handles().map(|(c, _)| m.amount(c).to_si()).collect();
            v.push(m.depot_mass().to_si());
            v.push(m.eliminated_mass().to_si());
            v.push(m.pending_dose().to_si());
            v
        };
        let cut = |window: f64, windows: usize| {
            let mut sim = Simulation::new(Schedule::Multirate)
                .conservation_tolerance(1e-12)
                .with(build());
            let mut substeps = 0;
            for _ in 0..windows {
                substeps = sim.advance(Time::s(window)).unwrap().substeps[0].1;
            }
            let m: &CompartmentModel = sim.domain_as("subject").unwrap();
            (state(m), substeps)
        };
        let limit = build().max_stable_dt(Time::ZERO).to_si();
        // The blood binds: `(90 + 50 + 10 + 6)/5 = 31.2/h`.
        assert!(
            (limit / (3600.0 / 31.2) - 1.0).abs() < 4.0 * f64::EPSILON,
            "{limit} s"
        );
        let (fine, n_fine) = cut(100.0, 72);
        let (coarse, n_coarse) = cut(200.0, 36);
        println!("  limit {limit:.1} s; substeps {n_fine} and {n_coarse}");
        assert_eq!((n_fine, n_coarse), (1, 2), "both cuts take 100 s substeps");
        assert_eq!(fine, coarse, "two cuts of one run");

        let dt = 100.0;
        let steps = 72;
        let mut straight = build();
        run(&mut straight, Time::ZERO, Time::s(dt), steps);
        let mut rewound = build();
        run(&mut rewound, Time::ZERO, Time::s(dt), steps / 3);
        rewound.checkpoint();
        let mid = state(&rewound);
        let resume = Time::s((steps / 3) as f64 * dt);
        run(&mut rewound, resume, Time::s(dt), steps - steps / 3);
        rewound.restore();
        assert_eq!(state(&rewound), mid, "restore");
        run(&mut rewound, resume, Time::s(dt), steps - steps / 3);
        assert_eq!(state(&rewound), state(&straight));
        assert_eq!(
            state(&straight),
            fine,
            "and the scheduler ran the same steps"
        );
        assert!(
            mid[5] > 0.0 && mid[4] > 0.0,
            "the run is not vacuous: {mid:?}"
        );
    }

    /// **A bad perfusion is refused where it is declared, and named.** A negative, infinite or
    /// NaN blood flow; a `K_p` that is zero, negative, infinite or NaN; a tissue perfused from
    /// itself; a tissue perfused twice; a tissue that is another's blood; a blood pool that is a
    /// tissue; a handle from another model. And nothing refused was half-built.
    #[test]
    fn a_bad_perfusion_is_refused_and_named() {
        let mut m = CompartmentModel::new("subject");
        let blood = m.eliminating("blood", Volume::litres(5.0), l_h(6.0));
        let liver = m.compartment("liver", Volume::litres(1.8));
        let fat = m.compartment("fat", Volume::litres(10.0));
        for q in [-1.0, f64::INFINITY, f64::NAN] {
            let err = m
                .perfuse(blood, liver, Clearance::from_si(q), 2.0)
                .expect_err("a bad blood flow");
            assert!(err.quantity.contains("blood flow"), "{err}");
            assert!(err.site.contains("blood-liver"), "{err}");
        }
        for k in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            let err = m
                .perfuse(blood, liver, l_h(90.0), k)
                .expect_err("a bad partition");
            assert!(err.quantity.contains("partition coefficient"), "{err}");
        }
        let err = m
            .perfuse(liver, liver, l_h(90.0), 2.0)
            .expect_err("from itself");
        assert!(err.quantity.contains("itself"), "{err}");
        let mut other = CompartmentModel::new("another subject");
        let elsewhere = other.compartment("liver", Volume::litres(1.8));
        let err = m
            .perfuse(blood, elsewhere, l_h(90.0), 2.0)
            .expect_err("another model");
        assert!(err.quantity.contains("different model"), "{err}");
        // Nothing refused above was half-built.
        assert_eq!(m.partition(liver), 1.0);
        assert_eq!(m.blood_flow(liver).to_si(), 0.0);
        assert_eq!(m.extraction_ratio(liver), None);

        // A zero flow is allowed — a tissue cut off — and its extraction is a number, not 0/0.
        m.perfuse(blood, liver, Clearance::ZERO, 2.0).unwrap();
        assert_eq!(m.extraction_ratio(liver), Some(0.0));
        let err = m
            .perfuse(blood, liver, l_h(90.0), 3.0)
            .expect_err("perfused twice");
        assert!(err.quantity.contains("already perfused"), "{err}");
        let err = m
            .perfuse(liver, fat, l_h(10.0), 5.0)
            .expect_err("a tissue as a blood pool");
        assert!(err.quantity.contains("cannot be a blood pool"), "{err}");
        let err = m
            .perfuse(fat, blood, l_h(10.0), 5.0)
            .expect_err("a blood pool as a tissue");
        assert!(err.quantity.contains("perfused as a tissue"), "{err}");
        assert_eq!(
            m.partition(liver),
            2.0,
            "the refused second perfusion left K_p alone"
        );
        assert_eq!(m.partition(blood), 1.0);
        assert_eq!(m.partition(fat), 1.0);
    }
}
