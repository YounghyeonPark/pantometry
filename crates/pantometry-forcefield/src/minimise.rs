//! Energy minimisation: L-BFGS with a backtracking line search, one iteration at a time.
//!
//! # The method, and why
//!
//! **Limited-memory BFGS** (Nocedal, *Math. Comp.* **35**, 773 (1980)), keeping the last
//! [`MEMORY`] position and gradient differences, with an **Armijo backtracking line search**:
//! a step `α d` along the L-BFGS direction `d` is accepted only if
//! `E(x + α d) ≤ E(x) + 10⁻⁴ α (g · d)`, halving α from 1 until it is.
//!
//! - Against **steepest descent with the same line search**, which is the simplest monotone
//!   method: a drug molecule's Hessian spans from a bond's ~700 kcal mol⁻¹ Å⁻² to a methyl
//!   rotor's ~1, and steepest descent's progress per step is set by that ratio. L-BFGS builds
//!   the curvature from the steps it has taken and does not pay it.
//! - Against **FIRE**, a damped dynamics: its energy is not monotone by construction (it
//!   overshoots and then quenches), and its convergence depends on a time step and five tuning
//!   constants. The line search gives a guarantee instead of a tuning.
//!
//! # What is guaranteed
//!
//! **The energy never increases.** A step is taken only when the Armijo condition holds, and
//! `g · d < 0` is checked before the search (the history is discarded and steepest descent used
//! when it is not), so every accepted step lowers the energy, by at least `10⁻⁴ α |g · d|`. When no
//! step length of 2⁻³⁰ or more lowers it — which happens near a tight minimum, where the decrease
//! a step can make falls below the rounding of the energy — the minimiser takes no step and
//! reports [`Status::Stalled`], with the force it stopped at. It never moves to a higher energy
//! to make progress.
//!
//! **Converged** means the largest force on any atom, as a vector, is at most the tolerance —
//! in newtons in this API ([`KCAL_PER_MOL_ANGSTROM`] converts). The criterion is checked before
//! each step, so a converged minimiser does not move.
//!
//! **Deterministic.** The unrestrained minimiser itself uses only `+ − × ÷` and `sqrt`, in a
//! fixed order, over `Vec`s: no clock, no randomness, no hash order. A [`DihedralRestraint`]
//! adds the platform's `atan2` on every step. The energy it minimises calls the platform's `cos`
//! three times per torsion evaluation (`cos nφ₀`, and the switch's constant edge once for each
//! central angle) and `ln`, `sin`, `cos` and `asin` when the force field is built. Those are
//! not correctly rounded and need not be the same function on two platforms, so the same start
//! gives the same bits on one machine at any optimisation level, and the crate does not claim
//! the same bits across platforms; see the crate documentation.
//!
//! # Restraints
//!
//! A [`DihedralRestraint`] adds `½ k (φ − φ₀)²`, the difference wrapped to (−π, π], to the energy
//! the minimiser sees and nothing else: [`Progress::energy`] is the force field's, and
//! [`Progress::restraint`] the restraint's, so a constrained scan reports the physical energy at a
//! dihedral held to within `τ / k` of its target, τ the torque the molecule puts on it. To second
//! order, the energy of the exactly-constrained minimum is `energy + 2 · restraint`, so the
//! restraint energy is also the size of the error a scan makes by not subtracting it.

use crate::ccd::ANGSTROM;
use crate::energy::ForceField;
use crate::uff::KCAL_PER_MOL;

/// One kcal mol⁻¹ Å⁻¹, the unit forces are reported in, in newtons.
pub const KCAL_PER_MOL_ANGSTROM: f64 = KCAL_PER_MOL / ANGSTROM;

/// How many past steps L-BFGS keeps.
pub const MEMORY: usize = 8;

/// The largest distance any atom moves in one step, ångström. A cap, not a trust region: it
/// stops a first step along a large clash force from throwing an atom across the molecule.
pub const MAX_STEP: f64 = 0.2;

/// The Armijo constant.
const ARMIJO: f64 = 1e-4;

/// How many times a step is halved before the search gives up: 2⁻³⁰ of the first trial.
const HALVINGS: usize = 30;

/// The first step's inverse curvature, Å² per kcal/mol: a guess of 100 kcal mol⁻¹ Å⁻², a stiff
/// bend's, used only until there is a measured step to scale by.
const FIRST_INVERSE_CURVATURE: f64 = 0.01;

/// Whether a minimisation is going, done, or cannot go on.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    /// Not converged yet, and the last step lowered the energy.
    Running,
    /// The largest force on any atom is at most the tolerance.
    Converged,
    /// No step lowers the energy, in floating point, and the force is still above the
    /// tolerance. The positions are the lowest the minimiser found.
    Stalled,
}

/// Where a minimisation is, after a step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Progress {
    /// The force field's energy at the current positions, joules per molecule — without any
    /// restraint.
    pub energy: f64,
    /// The restraints' energy, joules per molecule; zero without restraints.
    pub restraint: f64,
    /// The largest force on any atom, restraints included, newtons.
    pub max_force: f64,
    /// The root-mean-square of the per-atom forces, restraints included, newtons.
    pub rms_force: f64,
    /// How many steps have been taken.
    pub steps: usize,
    /// Running, converged or stalled.
    pub status: Status,
}

/// A harmonic restraint on the dihedral I–J–K–L: `½ k (φ − φ₀)²`, `φ − φ₀` wrapped to (−π, π].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DihedralRestraint {
    /// I, J, K, L as indices into the component.
    pub atoms: [usize; 4],
    /// φ₀, radians, in the sign convention of [`dihedral`].
    pub target: f64,
    /// k, joules per molecule per radian².
    pub stiffness: f64,
}

/// The dihedral I–J–K–L at positions `at`, radians in (−π, π]: `atan2(|b₂| b₁·(b₂×b₃),
/// (b₁×b₂)·(b₂×b₃))` with `b₁ = r_J − r_I`, `b₂ = r_K − r_J`, `b₃ = r_L − r_K` — 0 eclipsed, ±π
/// anti, positive clockwise looking down J→K (the IUPAC convention).
pub fn dihedral(at: &[[f64; 3]], [i, j, k, l]: [usize; 4]) -> f64 {
    let b1 = sub(at[j], at[i]);
    let b2 = sub(at[k], at[j]);
    let b3 = sub(at[l], at[k]);
    let m = cross(b1, b2);
    let n = cross(b2, b3);
    (norm(b2) * dot(b1, n)).atan2(dot(m, n))
}

impl DihedralRestraint {
    /// The restraint's energy at positions `at` (metres), joules per molecule.
    pub fn energy_at(&self, at: &[[f64; 3]]) -> f64 {
        let d = wrap(dihedral(at, self.atoms) - self.target);
        0.5 * self.stiffness * d * d
    }

    /// Adds the restraint's forces to `forces` and returns its energy. The dihedral's gradient,
    /// with `m = b₁ × b₂`, `n = b₂ × b₃`, `p = b₁·b₂/|b₂|²` and `q = b₃·b₂/|b₂|²`:
    /// `∂φ/∂r_I = −|b₂| m/|m|²`, `∂φ/∂r_L = |b₂| n/|n|²`,
    /// `∂φ/∂r_J = −(p + 1) ∂φ/∂r_I + q ∂φ/∂r_L` and `∂φ/∂r_K = −(q + 1) ∂φ/∂r_L + p ∂φ/∂r_I`,
    /// which sum to zero. Checked against a central difference in the tests.
    pub fn add_forces(&self, at: &[[f64; 3]], forces: &mut [[f64; 3]]) -> f64 {
        let [i, j, k, l] = self.atoms;
        let b1 = sub(at[j], at[i]);
        let b2 = sub(at[k], at[j]);
        let b3 = sub(at[l], at[k]);
        let m = cross(b1, b2);
        let n = cross(b2, b3);
        let (mm, nn, lb2) = (dot(m, m), dot(n, n), norm(b2));
        if mm == 0.0 || nn == 0.0 {
            return self.energy_at(at);
        }
        let d = wrap((lb2 * dot(b1, n)).atan2(dot(m, n)) - self.target);
        let de = self.stiffness * d;
        let gi = scale(m, -lb2 / mm);
        let gl = scale(n, lb2 / nn);
        let (p, q) = (dot(b1, b2) / (lb2 * lb2), dot(b3, b2) / (lb2 * lb2));
        let gj = add(scale(gi, -p - 1.0), scale(gl, q));
        let gk = add(scale(gl, -q - 1.0), scale(gi, p));
        for (a, g) in [(i, gi), (j, gj), (k, gk), (l, gl)] {
            for c in 0..3 {
                forces[a][c] -= de * g[c];
            }
        }
        0.5 * self.stiffness * d * d
    }
}

/// Wraps an angle difference to (−π, π].
fn wrap(d: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let mut d = d % TAU;
    if d > PI {
        d -= TAU;
    } else if d <= -PI {
        d += TAU;
    }
    d
}

/// One evaluation, in the minimiser's own units: positions in Å, energy in kcal/mol, gradient in
/// kcal mol⁻¹ Å⁻¹.
#[derive(Clone, Debug, PartialEq)]
struct Point {
    x: Vec<f64>,
    /// The positions evaluated, in metres, exactly as written back: what a later step's input is
    /// compared with, since metres → Å → metres does not round-trip bit for bit.
    at: Vec<[f64; 3]>,
    energy: f64,
    physical: f64,
    restraint: f64,
    gradient: Vec<f64>,
}

/// An L-BFGS minimiser that takes one step at a time, so that each step can be a frame. Holds
/// its curvature history between steps; [`Minimiser::step`] notices when the positions it is
/// given are not the ones it left, and starts its history again.
#[derive(Clone, Debug, PartialEq)]
pub struct Minimiser {
    tolerance: f64,
    memory: usize,
    history: Vec<(Vec<f64>, Vec<f64>)>,
    last: Option<Point>,
    steps: usize,
    status: Status,
}

impl Minimiser {
    /// A minimiser that converges when the largest force on any atom is at most `tolerance`
    /// newtons.
    pub fn new(tolerance: f64) -> Minimiser {
        Minimiser::with_memory(tolerance, MEMORY)
    }

    /// A minimiser keeping `memory` past steps instead of [`MEMORY`]. With 0 it is steepest descent
    /// with the same line search — the baseline that shows what the history buys.
    pub fn with_memory(tolerance: f64, memory: usize) -> Minimiser {
        Minimiser {
            tolerance,
            memory,
            history: Vec::new(),
            last: None,
            steps: 0,
            status: Status::Running,
        }
    }

    /// The tolerance, newtons.
    pub fn tolerance(&self) -> f64 {
        self.tolerance
    }

    /// The status after the last step.
    pub fn status(&self) -> Status {
        self.status
    }

    /// How many steps have been taken.
    pub fn steps(&self) -> usize {
        self.steps
    }

    /// How many past steps the curvature estimate is currently built from, at most [`MEMORY`].
    /// It resets when the minimiser is handed positions it did not leave, or when its direction
    /// fails to go downhill.
    pub fn history_len(&self) -> usize {
        self.history.len()
    }

    /// One iteration: evaluates `ff` (plus `restraints`) at `at`, and unless converged or
    /// stalled, moves `at` along the L-BFGS direction by the longest step of 1, ½, ¼, … that
    /// satisfies the Armijo condition. Returns where it then is.
    pub fn step(
        &mut self,
        ff: &ForceField,
        restraints: &[DihedralRestraint],
        at: &mut [[f64; 3]],
    ) -> Progress {
        let here = match self.last.take() {
            Some(p) if p.at == at => p,
            _ => {
                let x = flatten(at);
                self.history.clear();
                if self.status == Status::Stalled {
                    self.status = Status::Running;
                }
                evaluate(ff, restraints, &x)
            }
        };
        if max_atom(&here.gradient) <= self.tolerance / KCAL_PER_MOL_ANGSTROM {
            self.status = Status::Converged;
        }
        if self.status != Status::Running {
            let progress = self.progress(&here);
            self.last = Some(here);
            return progress;
        }
        let next = match self.search(ff, restraints, &here, true) {
            Some(p) => Some(p),
            None if !self.history.is_empty() => {
                self.history.clear();
                self.search(ff, restraints, &here, false)
            }
            None => None,
        };
        let now = match next {
            Some(p) => {
                let s: Vec<f64> = p.x.iter().zip(&here.x).map(|(a, b)| a - b).collect();
                let y: Vec<f64> = p
                    .gradient
                    .iter()
                    .zip(&here.gradient)
                    .map(|(a, b)| a - b)
                    .collect();
                // Keep a pair only if it has positive curvature, which BFGS needs to stay
                // positive definite; the Armijo search alone does not ensure it.
                if dot_n(&s, &y) > 0.0 && self.memory > 0 {
                    if self.history.len() == self.memory {
                        self.history.remove(0);
                    }
                    self.history.push((s, y));
                }
                self.steps += 1;
                at.copy_from_slice(&p.at);
                if max_atom(&p.gradient) <= self.tolerance / KCAL_PER_MOL_ANGSTROM {
                    self.status = Status::Converged;
                }
                p
            }
            None => {
                self.status = Status::Stalled;
                here
            }
        };
        let progress = self.progress(&now);
        self.last = Some(now);
        progress
    }

    /// The L-BFGS direction (or steepest descent, when `use_history` is false or the direction
    /// is not downhill), capped, and the Armijo search along it. `None` if no step length
    /// lowers the energy.
    fn search(
        &self,
        ff: &ForceField,
        restraints: &[DihedralRestraint],
        here: &Point,
        use_history: bool,
    ) -> Option<Point> {
        let g = &here.gradient;
        let mut d = if use_history && !self.history.is_empty() {
            self.direction(g)
        } else {
            g.iter().map(|v| -FIRST_INVERSE_CURVATURE * v).collect()
        };
        let mut slope = dot_n(g, &d);
        if slope >= 0.0 {
            d = g.iter().map(|v| -FIRST_INVERSE_CURVATURE * v).collect();
            slope = dot_n(g, &d);
        }
        if slope >= 0.0 {
            return None;
        }
        let longest = max_atom(&d);
        if longest > MAX_STEP {
            let f = MAX_STEP / longest;
            for v in &mut d {
                *v *= f;
            }
            slope *= f;
        }
        let mut alpha = 1.0;
        for _ in 0..=HALVINGS {
            let x: Vec<f64> = here.x.iter().zip(&d).map(|(a, b)| a + alpha * b).collect();
            let trial = evaluate(ff, restraints, &x);
            if trial.energy <= here.energy + ARMIJO * alpha * slope && trial.energy < here.energy {
                return Some(trial);
            }
            alpha *= 0.5;
        }
        None
    }

    /// The two-loop recursion: `−H g`, with `H₀ = (s·y / y·y) I` from the newest pair.
    fn direction(&self, g: &[f64]) -> Vec<f64> {
        two_loop(&self.history, g)
    }

    fn progress(&self, p: &Point) -> Progress {
        let n = p.gradient.len() / 3;
        let sum: f64 = (0..n)
            .map(|a| {
                let v = &p.gradient[3 * a..3 * a + 3];
                v[0] * v[0] + v[1] * v[1] + v[2] * v[2]
            })
            .sum();
        Progress {
            energy: p.physical * KCAL_PER_MOL,
            restraint: p.restraint * KCAL_PER_MOL,
            max_force: max_atom(&p.gradient) * KCAL_PER_MOL_ANGSTROM,
            rms_force: (sum / n.max(1) as f64).sqrt() * KCAL_PER_MOL_ANGSTROM,
            steps: self.steps,
            status: self.status,
        }
    }
}

impl ForceField {
    /// Minimises the energy from positions `at` (metres) in place, for at most `max_steps`
    /// steps or until the largest force on any atom is at most `tolerance` newtons, and says
    /// where it stopped. See [`crate::minimise`] for the method and what it guarantees.
    pub fn minimise(&self, at: &mut [[f64; 3]], max_steps: usize, tolerance: f64) -> Progress {
        self.minimise_restrained(at, &[], max_steps, tolerance)
    }

    /// [`ForceField::minimise`] with `restraints` added to the energy minimised.
    /// [`Progress::energy`] is still the force field's alone.
    pub fn minimise_restrained(
        &self,
        at: &mut [[f64; 3]],
        restraints: &[DihedralRestraint],
        max_steps: usize,
        tolerance: f64,
    ) -> Progress {
        let mut m = Minimiser::new(tolerance);
        let mut p = m.step(self, restraints, at);
        while p.status == Status::Running && p.steps < max_steps {
            p = m.step(self, restraints, at);
        }
        p
    }
}

/// The L-BFGS product `−H g` by the two-loop recursion (Nocedal 1980), with `H₀ = (s·y / y·y) I`
/// from the newest pair: the inverse-Hessian estimate that applying the BFGS update
/// `H ← (I − ρ s yᵀ) H (I − ρ y sᵀ) + ρ s sᵀ`, `ρ = 1/(y·s)`, for each pair in `pairs` from the
/// oldest to the newest would build, times `−g`. Public so that it can be checked against that
/// dense update directly; every pair must have `s · y > 0`, which [`Minimiser`] ensures.
///
/// # Panics
///
/// If `pairs` is empty.
pub fn two_loop(pairs: &[(Vec<f64>, Vec<f64>)], g: &[f64]) -> Vec<f64> {
    let mut q = g.to_vec();
    let mut alphas = Vec::with_capacity(pairs.len());
    for (s, y) in pairs.iter().rev() {
        let rho = 1.0 / dot_n(y, s);
        let a = rho * dot_n(s, &q);
        for (qi, yi) in q.iter_mut().zip(y) {
            *qi -= a * yi;
        }
        alphas.push((rho, a));
    }
    let (s, y) = pairs.last().expect("called with a history");
    let gamma = dot_n(s, y) / dot_n(y, y);
    for qi in &mut q {
        *qi *= gamma;
    }
    for ((s, y), (rho, a)) in pairs.iter().zip(alphas.iter().rev()) {
        let b = rho * dot_n(y, &q);
        for (qi, si) in q.iter_mut().zip(s) {
            *qi += (a - b) * si;
        }
    }
    q.iter().map(|v| -v).collect()
}

fn evaluate(ff: &ForceField, restraints: &[DihedralRestraint], x: &[f64]) -> Point {
    let mut at = vec![[0.0; 3]; x.len() / 3];
    unflatten(x, &mut at);
    let mut ev = ff.evaluate(&at);
    let mut restraint = 0.0;
    for r in restraints {
        restraint += r.add_forces(&at, &mut ev.forces);
    }
    let at = at;
    let physical = ev.energy.total / KCAL_PER_MOL;
    let restraint = restraint / KCAL_PER_MOL;
    let gradient = ev
        .forces
        .iter()
        .flat_map(|f| f.iter().map(|v| -v / KCAL_PER_MOL_ANGSTROM))
        .collect();
    Point {
        x: x.to_vec(),
        at,
        energy: physical + restraint,
        physical,
        restraint,
        gradient,
    }
}

fn flatten(at: &[[f64; 3]]) -> Vec<f64> {
    at.iter()
        .flat_map(|p| p.iter().map(|v| v / ANGSTROM))
        .collect()
}

fn unflatten(x: &[f64], at: &mut [[f64; 3]]) {
    for (a, p) in at.iter_mut().enumerate() {
        *p = [
            x[3 * a] * ANGSTROM,
            x[3 * a + 1] * ANGSTROM,
            x[3 * a + 2] * ANGSTROM,
        ];
    }
}

/// The largest per-atom vector norm in a flat 3N vector.
fn max_atom(v: &[f64]) -> f64 {
    (0..v.len() / 3)
        .map(|a| {
            let c = &v[3 * a..3 * a + 3];
            (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()
        })
        .fold(0.0, f64::max)
}

fn dot_n(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}
