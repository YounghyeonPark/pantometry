//! One molecule in, one verdict out: absorption, binding and distribution read together.
//!
//! This is the one thing the facade implements, and it is here because it has to be: the verdict
//! reads a [`Component`] from `pantometry-forcefield` and a tissue partition from
//! `pantometry-pharmacokinetic`, and no domain may depend on another.
//!
//! ```
//! use pantometry::candidate::{assess, Call, Evidence, Measured, Rule, Thresholds};
//! use pantometry::forcefield::Component;
//!
//! // The dictionary entry's text is the caller's to supply. This path is this repository's.
//! let benzene = Component::from_ccd(include_str!(
//!     "../../pantometry-forcefield/components/BNZ.cif"
//! ))
//! .unwrap();
//! let evidence = Evidence {
//!     log_p: Some(Measured::new(1.92, 0.23)),
//!     ..Evidence::default()
//! };
//! // The thresholds are the caller's. See `Thresholds` for why none is built in.
//! let thresholds = Thresholds { violations_for_alert: 2, ligand_efficiency: 0.3 };
//! let a = assess(&benzene, &evidence, &thresholds).unwrap();
//! assert_eq!(a.absorption.call, Call::Pass);
//! assert_eq!(a.absorption.limit(Rule::LogP).call, Call::Pass);
//! assert_eq!(a.verdict, Call::NotMeasured); // no binding free energy was given
//! ```
//!
//! # What it reads, and what it does not compute
//!
//! From the molecule itself, exactly: the molecular weight, the hydrogen-bond donors and
//! acceptors, and the heavy atoms ([`Descriptors`]). Everything else is **evidence the caller
//! brings**, each with its standard error, because each costs a free-energy calculation of a day
//! or more and this function is not the place to start one:
//!
//! - **log P**, octanol/water. In this repository's own measurements, benzene is 1.92 ± 0.23 and
//!   p-xylene 3.06 ± 0.22 against experiment's 2.13 and 3.15.
//! - **the standard binding free energy** `ΔG°_bind`, kcal/mol at a stated temperature. Benzene's
//!   in T4 lysozyme L99A is −8.975 ± 0.279 with UFF throughout and −4.61 ± 0.29 with GAFF,
//!   ff14SB and AM1-BCC, against experiment's −5.19 ± 0.16. Those two are a factor of 1 583 apart
//!   in `K_d`, and benzene passes the binding screen with either, because it is small and both are
//!   efficient per heavy atom. The verdict cannot know which method made its evidence.
//! - **the fraction unbound in plasma** `f_u`, which a tissue partition needs and nothing here can
//!   compute.
//!
//! Anything not given is [`Call::NotMeasured`], never a default. That is a statement about each
//! value, not about each screen: a screen whose missing value could not change its outcome is
//! still decided. See [`Absorption::call`].
//!
//! # What it refuses
//!
//! **A molecule whose atoms are not its stated formula**, when the entry states one. A component
//! read without its hydrogens would otherwise come out lighter and with no donors at all, and pass.
//! A component with no stated formula, such as one assembled from a PDB entry, is not checked.
//!
//! # How an uncertain number is called
//!
//! A value is the interval `value ± 2σ` ([`INTERVAL_SIGMAS`]). It **passes** a limit if the whole
//! interval is on the right side, **fails** if the whole interval is on the wrong side, and is
//! **undecided** if the limit falls inside it. Two standard errors is the interval this
//! repository's pre-registered log P gates have used. A count of atoms has no uncertainty, so it
//! is never undecided.
//!
//! # The criteria, and which of them is checked against its source
//!
//! - **Absorption: Lipinski's rule of five.** The limits are those of the abstract of Lipinski,
//!   Lombardo, Dominy and Feeney, *Adv Drug Deliv Rev* 2001;**46**:3–26 (the 1997 paper,
//!   reprinted), doi:10.1016/S0169-409X(00)00129-0: poor absorption or permeation is more likely
//!   "when there are more than 5 H-bond donors, 10 H-bond acceptors, the molecular weight (MWT) is
//!   greater than 500 and the calculated Log P (CLogP) is greater than 5 (or MlogP > 4.15)". That
//!   text was read, from Europe PMC's copy of the abstract. Its log P is a calculated one, CLogP,
//!   and this module applies the CLogP limit of 5 to whatever log P the caller measured.
//! - **How the donors and acceptors are counted is this module's reading**, not a quotation. The
//!   paper's text was not available to check. A donor here is each hydrogen bonded to a nitrogen
//!   or an oxygen, so an NH₂ is two; an acceptor is each nitrogen and each oxygen. That is the
//!   common reading, and it is the one [`Descriptors`] documents.
//! - **How many exceeded limits raise the alert** is [`Thresholds::violations_for_alert`]. It is
//!   usually quoted as two, but the abstract does not say it and the paper was not read, so it is
//!   not built in.
//! - **Binding: ligand efficiency**, `LE = −ΔG°_bind / N_heavy` (Hopkins, Groom and Alex, *Drug
//!   Discov Today* 2004;**9**:430–431, doi:10.1016/S1359-6446(04)03069-7). The least acceptable
//!   value is [`Thresholds::ligand_efficiency`]. It is usually quoted as about 0.3 kcal/mol per
//!   heavy atom, but the paper was not read, so it is not built in.
//! - **Distribution is reported and not judged.** `K_p` for rat liver, kidney and muscle comes from
//!   [`tissue_partition`] across the log P interval. No threshold on `K_p` is published that this
//!   module could cite. Adipose is left out, because Rodgers and Rowland give it a vegetable
//!   oil:water coefficient and not log P.
//! - **That equation is for a neutral molecule.** A component with a net formal charge gets
//!   [`Distribution::Charged`] and no `K_p`. A molecule with no formal charge that ionises at
//!   physiological pH, as an acid such as aspirin does, cannot be recognised from its atoms. Its
//!   `K_p` here would be the neutral equation's, which Rodgers and Rowland do not use for it, and
//!   not giving it a fraction unbound is the caller's way to leave it out.
//!
//! The **verdict** is the worst of the absorption call and the binding call, in [`Call`]'s order.
//! A molecule passes only if both were measured and both passed. Nothing here says it is a drug,
//! since passing two screens is only passing two screens.

use pantometry_core::math;
use pantometry_forcefield::uff::KCAL_PER_MOL;
use pantometry_forcefield::{Component, Element};
use pantometry_pharmacokinetic::partition::{
    tissue_partition, TissueComposition, RAT_KIDNEY, RAT_LIVER, RAT_MUSCLE, RAT_PLASMA_LIPIDS,
};
use pantometry_units::{Temperature, BOLTZMANN};
use std::fmt;

/// How many standard errors either side of a value its interval spans: two.
pub const INTERVAL_SIGMAS: f64 = 2.0;

/// The standard concentration a standard binding free energy is referred to: 1 mol/L.
pub const STANDARD_CONCENTRATION: f64 = 1.0;

/// A measured value and its standard error, in whatever unit the field that holds it names.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Measured {
    /// The value.
    pub value: f64,
    /// Its standard error, `σ`. Zero for a value taken as exact.
    pub sigma: f64,
}

impl Measured {
    /// A value with its standard error.
    pub fn new(value: f64, sigma: f64) -> Measured {
        Measured { value, sigma }
    }

    /// `[value − 2σ, value + 2σ]`.
    pub fn interval(self) -> [f64; 2] {
        let half = INTERVAL_SIGMAS * self.sigma;
        [self.value - half, self.value + half]
    }
}

impl fmt::Display for Measured {
    /// Two decimal places unless the formatter is given a precision: `{:.3}` gives three.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let p = f.precision().unwrap_or(2);
        write!(f, "{:.p$} ± {:.p$}", self.value, self.sigma)
    }
}

/// What a criterion, a screen or the whole verdict came to, ordered from best to worst, so the
/// worst of several is their `max`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Call {
    /// The whole interval is on the acceptable side.
    Pass,
    /// No value was given, and the value could change the outcome.
    NotMeasured,
    /// The limit falls inside the interval, so the measurement cannot decide it.
    Undecided,
    /// The whole interval is on the wrong side.
    Fail,
}

impl fmt::Display for Call {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Call::Pass => "pass",
            Call::Fail => "fail",
            Call::Undecided => "undecided",
            Call::NotMeasured => "not measured",
        })
    }
}

/// What is read off the molecule itself, with no uncertainty.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Descriptors {
    /// The sum of the atoms' standard atomic weights ([`Element::atomic_weight`]), g/mol.
    pub molecular_weight: f64,
    /// Each hydrogen bonded to a nitrogen or an oxygen. An NH₂ is two, a hydroxyl one.
    pub donors: usize,
    /// Each nitrogen and each oxygen, whatever its bonds.
    pub acceptors: usize,
    /// Every atom that is not hydrogen.
    pub heavy_atoms: usize,
}

impl Descriptors {
    /// Read them off a component, which must carry its hydrogens as the dictionary does.
    ///
    /// A component without its hydrogens gives a molecular weight too low and no donors at all.
    /// That is not detected here. [`assess`] detects it when the entry states a formula.
    pub fn of(component: &Component) -> Descriptors {
        let atoms = component.atoms();
        let polar = |e: Element| matches!(e, Element::N | Element::O);
        let donors = component
            .bonds()
            .iter()
            .filter(|b| {
                let [a, z] = b.atoms.map(|i| atoms[i].element);
                (a == Element::H && polar(z)) || (z == Element::H && polar(a))
            })
            .count();
        Descriptors {
            molecular_weight: atoms.iter().map(|a| a.element.atomic_weight()).sum(),
            donors,
            acceptors: atoms.iter().filter(|a| polar(a.element)).count(),
            heavy_atoms: atoms.iter().filter(|a| a.element != Element::H).count(),
        }
    }
}

/// The thresholds this module does not build in, because their sources were not read.
///
/// There is deliberately no `Default`. Both are usually quoted, as two violations and about 0.3
/// kcal/mol per heavy atom, but a default would state a source's number without the source. A
/// caller who uses them is making that choice, and the call site says so.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Thresholds {
    /// How many of the rule of five's four limits must be exceeded before absorption fails.
    /// One to four.
    pub violations_for_alert: usize,
    /// The least ligand efficiency a binder must reach, kcal/mol per heavy atom. Positive.
    pub ligand_efficiency: f64,
}

/// A standard binding free energy and the temperature it was measured at.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BindingFreeEnergy {
    /// `ΔG°_bind`, kcal/mol, referred to [`STANDARD_CONCENTRATION`]. Negative binds. A value in
    /// kJ/mol is 4.184 times too large and is not detected.
    pub kcal_per_mol: Measured,
    /// The temperature, as a type, so 25 °C is `Temperature::celsius(25.0)` and not 25 K.
    pub temperature: Temperature,
}

/// What the caller measured. `Evidence::default()` is nothing measured.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Evidence {
    /// log P, octanol/water.
    pub log_p: Option<Measured>,
    /// The standard binding free energy to the target.
    pub binding: Option<BindingFreeEnergy>,
    /// The fraction unbound in plasma, in `(0, 1]`, taken as exact. Give it only for a molecule
    /// that is neutral at physiological pH: see the [module documentation](self).
    pub fraction_unbound: Option<f64>,
}

/// Why [`assess`] refused. Each variant names what it refused.
#[derive(Clone, Debug, PartialEq)]
pub enum AssessError {
    /// A value or a standard error that is not finite, or a standard error below zero.
    NotAMeasurement {
        /// Which input.
        what: &'static str,
        /// What it was.
        value: Measured,
    },
    /// A fraction unbound outside `(0, 1]`.
    FractionUnbound(f64),
    /// A temperature that is not finite and positive, in kelvin.
    Temperature(f64),
    /// A ligand-efficiency threshold that is not finite and positive.
    LigandEfficiency(f64),
    /// An alert count outside `1..=4`. Zero would fail every molecule, and more than four can
    /// never be reached.
    Violations(usize),
    /// A molecule with no heavy atom, whose ligand efficiency would divide by zero.
    NoHeavyAtoms,
    /// The atoms read are not the formula the entry states, most often because the hydrogens are
    /// missing.
    Formula {
        /// `_chem_comp.formula`, as the file has it.
        stated: String,
        /// The atoms read, in the same notation.
        read: String,
    },
}

impl fmt::Display for AssessError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            AssessError::NotAMeasurement { what, value } => write!(
                f,
                "{what} is {} ± {}, which is not a measurement: both must be finite and the \
                 standard error at least zero",
                value.value, value.sigma
            ),
            AssessError::FractionUnbound(fu) => {
                write!(f, "a fraction unbound of {fu} is outside (0, 1]")
            }
            AssessError::Temperature(t) => {
                write!(f, "a temperature of {t} K is not finite and positive")
            }
            AssessError::LigandEfficiency(x) => {
                write!(
                    f,
                    "a ligand-efficiency threshold of {x} is not finite and positive"
                )
            }
            AssessError::Violations(n) => write!(
                f,
                "an alert at {n} violations is outside 1..=4: the rule of five has four limits"
            ),
            AssessError::NoHeavyAtoms => {
                f.write_str("the molecule has no heavy atom, so it has no ligand efficiency")
            }
            AssessError::Formula { stated, read } => write!(
                f,
                "the entry states the formula {stated} and its atoms are {read}; a component \
                 without its hydrogens would be judged on the wrong counts"
            ),
        }
    }
}

impl std::error::Error for AssessError {}

/// The rule of five's four limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Rule {
    /// At most five H-bond donors.
    Donors,
    /// At most ten H-bond acceptors.
    Acceptors,
    /// A molecular weight of at most 500.
    MolecularWeight,
    /// A log P of at most 5.
    LogP,
}

impl Rule {
    /// The four, in the order [`Absorption::limits`] holds them.
    pub const ALL: [Rule; 4] = [
        Rule::Donors,
        Rule::Acceptors,
        Rule::MolecularWeight,
        Rule::LogP,
    ];

    /// The most that passes. Exceeding it is a violation.
    pub fn most(self) -> f64 {
        match self {
            Rule::Donors => 5.0,
            Rule::Acceptors => 10.0,
            Rule::MolecularWeight => 500.0,
            Rule::LogP => 5.0,
        }
    }
}

impl fmt::Display for Rule {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Rule::Donors => "H-bond donors",
            Rule::Acceptors => "H-bond acceptors",
            Rule::MolecularWeight => "molecular weight",
            Rule::LogP => "log P",
        })
    }
}

/// One of the rule of five's limits, and where the molecule stands against it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Limit {
    /// Which limit.
    pub rule: Rule,
    /// The molecule's value, if it was measured. A count has a standard error of zero.
    pub value: Option<Measured>,
    /// Pass if the whole interval is at or below [`Rule::most`], fail if the whole interval is
    /// above it.
    pub call: Call,
}

/// The absorption screen: the rule of five.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Absorption {
    /// In [`Rule::ALL`]'s order. [`Absorption::limit`] finds one by name.
    pub limits: [Limit; 4],
    /// How many limits are certainly exceeded.
    pub exceeded: usize,
    /// How many are undecided or not measured, and so might be exceeded.
    pub open: usize,
    /// Fail if `exceeded` alone reaches the alert. Pass if `exceeded + open` cannot. Otherwise
    /// undecided, or not measured if every open limit is unmeasured.
    ///
    /// **A pass does not mean every limit was looked at.** With the alert at two, a molecule
    /// under the three counted limits passes with no log P at all, because one open limit cannot
    /// reach two. `open` says how many were not decided, and the log P limit's own call says
    /// whether it was measured.
    pub call: Call,
}

impl Absorption {
    /// The limit for one rule.
    pub fn limit(&self, rule: Rule) -> &Limit {
        &self.limits[rule as usize]
    }
}

/// The binding screen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Affinity {
    /// `K_d = c° exp(ΔG°/RT)`, mol/L, at the measured value.
    pub dissociation_constant: f64,
    /// `K_d` at `ΔG° − 2σ` and at `ΔG° + 2σ`, in that order, so tighter first.
    pub dissociation_interval: [f64; 2],
    /// `−ΔG° / N_heavy`, kcal/mol per heavy atom, with its standard error.
    pub ligand_efficiency: Measured,
    /// Pass if the whole interval reaches the threshold, fail if the whole interval falls short.
    pub call: Call,
}

/// The three tissues a partition is reported for, all rat.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Tissue {
    /// [`RAT_LIVER`].
    Liver,
    /// [`RAT_KIDNEY`].
    Kidney,
    /// [`RAT_MUSCLE`].
    Muscle,
}

impl Tissue {
    /// The three, in the order [`Distribution::Reported`] holds them.
    pub const ALL: [Tissue; 3] = [Tissue::Liver, Tissue::Kidney, Tissue::Muscle];

    /// Its composition, from `pantometry-pharmacokinetic`.
    pub fn composition(self) -> TissueComposition {
        match self {
            Tissue::Liver => RAT_LIVER,
            Tissue::Kidney => RAT_KIDNEY,
            Tissue::Muscle => RAT_MUSCLE,
        }
    }
}

impl fmt::Display for Tissue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Tissue::Liver => "liver",
            Tissue::Kidney => "kidney",
            Tissue::Muscle => "muscle",
        })
    }
}

/// One tissue's partition coefficient across the log P interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Partition {
    /// Which tissue.
    pub tissue: Tissue,
    /// `K_p`, tissue to plasma, at log P.
    pub kp: f64,
    /// The least and the greatest `K_p` over `log P ± 2σ`.
    pub kp_interval: [f64; 2],
}

/// What the distribution came to. Reported, never judged.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Distribution {
    /// Liver, kidney and muscle, in [`Tissue::ALL`]'s order.
    Reported([Partition; 3]),
    /// log P or the fraction unbound, or both, was not given.
    NotMeasured,
    /// The component carries this net formal charge, and the equation is for a neutral molecule.
    Charged(i32),
}

impl Distribution {
    /// The partition for one tissue, if any was reported.
    pub fn partition(&self, tissue: Tissue) -> Option<&Partition> {
        match self {
            Distribution::Reported(p) => Some(&p[tissue as usize]),
            _ => None,
        }
    }
}

/// Everything [`assess`] found, and the verdict.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Assessment {
    /// What was read off the molecule.
    pub descriptors: Descriptors,
    /// The rule of five.
    pub absorption: Absorption,
    /// The binding screen, if a binding free energy was given.
    pub affinity: Option<Affinity>,
    /// The tissue partitions. Not judged.
    pub distribution: Distribution,
    /// The worst of the absorption call and the binding call.
    pub verdict: Call,
}

fn check(what: &'static str, m: Measured) -> Result<Measured, AssessError> {
    if m.value.is_finite() && m.sigma.is_finite() && m.sigma >= 0.0 {
        Ok(m)
    } else {
        Err(AssessError::NotAMeasurement { what, value: m })
    }
}

fn limit(rule: Rule, value: Option<Measured>) -> Limit {
    let most = rule.most();
    let call = match value {
        None => Call::NotMeasured,
        Some(m) => {
            let [lo, hi] = m.interval();
            if hi <= most {
                Call::Pass
            } else if lo > most {
                Call::Fail
            } else {
                Call::Undecided
            }
        }
    };
    Limit { rule, value, call }
}

/// The antilog, `10^x`, through the kernel's `exp` so it is the same bits everywhere.
fn antilog(x: f64) -> f64 {
    math::exp(x * std::f64::consts::LN_10)
}

/// A formula as `(element, count)` in element order, or `None` if a word is not one.
fn formula_counts(formula: &str) -> Option<Vec<(Element, usize)>> {
    let mut out = Vec::new();
    for word in formula.split_whitespace() {
        let split = word
            .find(|c: char| c.is_ascii_digit())
            .unwrap_or(word.len());
        let (symbol, n) = word.split_at(split);
        let element = Element::from_symbol(symbol)?;
        let n = if n.is_empty() { 1 } else { n.parse().ok()? };
        out.push((element, n));
    }
    out.sort();
    Some(out)
}

/// The atoms read, in the dictionary's notation: `C6 H6`.
fn written(counts: &[(Element, usize)]) -> String {
    counts
        .iter()
        .map(|&(e, n)| {
            if n == 1 {
                e.symbol().to_uppercase()
            } else {
                format!("{}{n}", e.symbol().to_uppercase())
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

/// Judge one molecule on what the caller measured.
///
/// The component must carry its hydrogens. When the entry states a formula, the atoms are
/// checked against it, and a component without its hydrogens is refused.
///
/// # Errors
///
/// Any [`AssessError`]: a measurement that is not finite or has a negative standard error, a
/// fraction unbound outside `(0, 1]`, a temperature that is not finite and positive, a
/// ligand-efficiency threshold that is not finite and positive, an alert count outside `1..=4`,
/// atoms that are not the stated formula, or a molecule with no heavy atom given a binding free
/// energy.
pub fn assess(
    component: &Component,
    evidence: &Evidence,
    thresholds: &Thresholds,
) -> Result<Assessment, AssessError> {
    let le_threshold = thresholds.ligand_efficiency;
    if !(le_threshold.is_finite() && le_threshold > 0.0) {
        return Err(AssessError::LigandEfficiency(le_threshold));
    }
    if !(1..=4).contains(&thresholds.violations_for_alert) {
        return Err(AssessError::Violations(thresholds.violations_for_alert));
    }
    let log_p = evidence.log_p.map(|m| check("log P", m)).transpose()?;
    if let Some(b) = evidence.binding {
        check("the binding free energy", b.kcal_per_mol)?;
        let t = b.temperature.to_si();
        if !(t.is_finite() && t > 0.0) {
            return Err(AssessError::Temperature(t));
        }
    }
    if let Some(fu) = evidence.fraction_unbound {
        if !(fu > 0.0 && fu <= 1.0) {
            return Err(AssessError::FractionUnbound(fu));
        }
    }
    if let Some(stated) = component.formula() {
        let read: Vec<(Element, usize)> = component.composition().into_iter().collect();
        if formula_counts(stated).as_deref() != Some(read.as_slice()) {
            return Err(AssessError::Formula {
                stated: stated.to_string(),
                read: written(&read),
            });
        }
    }

    let d = Descriptors::of(component);
    let exact = |x: f64| Some(Measured::new(x, 0.0));
    let limits = [
        limit(Rule::Donors, exact(d.donors as f64)),
        limit(Rule::Acceptors, exact(d.acceptors as f64)),
        limit(Rule::MolecularWeight, exact(d.molecular_weight)),
        limit(Rule::LogP, log_p),
    ];
    let exceeded = limits.iter().filter(|l| l.call == Call::Fail).count();
    let open: Vec<Call> = limits
        .iter()
        .map(|l| l.call)
        .filter(|c| matches!(c, Call::Undecided | Call::NotMeasured))
        .collect();
    let alert = thresholds.violations_for_alert;
    let absorption_call = if exceeded >= alert {
        Call::Fail
    } else if exceeded + open.len() < alert {
        Call::Pass
    } else if open.contains(&Call::Undecided) {
        Call::Undecided
    } else {
        Call::NotMeasured
    };
    let absorption = Absorption {
        limits,
        exceeded,
        open: open.len(),
        call: absorption_call,
    };

    let affinity = match evidence.binding {
        None => None,
        Some(b) => {
            if d.heavy_atoms == 0 {
                return Err(AssessError::NoHeavyAtoms);
            }
            // RT in kcal/mol: k_B T over one kcal/mol per molecule.
            let rt = BOLTZMANN.to_si() * b.temperature.to_si() / KCAL_PER_MOL;
            let kd = |dg: f64| STANDARD_CONCENTRATION * math::exp(dg / rt);
            let [lo, hi] = b.kcal_per_mol.interval();
            let n = d.heavy_atoms as f64;
            let le = Measured::new(-b.kcal_per_mol.value / n, b.kcal_per_mol.sigma / n);
            let [le_lo, le_hi] = le.interval();
            let call = if le_lo >= le_threshold {
                Call::Pass
            } else if le_hi < le_threshold {
                Call::Fail
            } else {
                Call::Undecided
            };
            Some(Affinity {
                dissociation_constant: kd(b.kcal_per_mol.value),
                dissociation_interval: [kd(lo), kd(hi)],
                ligand_efficiency: le,
                call,
            })
        }
    };

    let charge: i32 = component.atoms().iter().map(|a| a.charge).sum();
    let distribution = match (log_p, evidence.fraction_unbound) {
        (Some(lp), Some(fu)) if charge == 0 => {
            let [lo, hi] = lp.interval();
            Distribution::Reported(Tissue::ALL.map(|tissue| {
                let t = tissue.composition();
                let kp = |x: f64| tissue_partition(&t, &RAT_PLASMA_LIPIDS, antilog(x), fu);
                let (a, z) = (kp(lo), kp(hi));
                Partition {
                    tissue,
                    kp: kp(lp.value),
                    kp_interval: [a.min(z), a.max(z)],
                }
            }))
        }
        (Some(_), Some(_)) => Distribution::Charged(charge),
        _ => Distribution::NotMeasured,
    };

    let verdict = absorption
        .call
        .max(affinity.map_or(Call::NotMeasured, |a| a.call));
    Ok(Assessment {
        descriptors: d,
        absorption,
        affinity,
        distribution,
        verdict,
    })
}
