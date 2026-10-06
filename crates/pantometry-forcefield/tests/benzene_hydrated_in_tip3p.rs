//! **Benzene decoupled from TIP3P water**: the alchemical Hamiltonian in a periodic box with
//! constrained dynamics, checked against exact identities, and — ignored, in release — benzene's
//! hydration free energy set beside experiment. Step W3 of the explicit-water track.
//!
//! Default, unoptimised:
//!
//! - **the end states are the systems they say**: fully coupled, benzene in the box less its own
//!   images, its own pairs in vacuum; fully decoupled, the water alone in its box plus benzene
//!   alone in vacuum, energy and every force;
//! - **`∂U/∂λ_e` and `∂U/∂λ_v` are the energy's derivatives** at a configuration constrained
//!   dynamics produced, against central differences of the whole energy;
//! - **every state of `couplings` is that state's `coupling` to the bit**, and the differences
//!   between states are the whole energy's;
//! - **a charge held at a distance from a frozen one decouples to its closed form**: the sphere
//!   average of the cross Ewald energy is `k_e q Q (1/d + ξ/L + 2πd²/3L³)`, Wigner's constant and
//!   the background's curvature, and the free energy the windows sample is the quadrature's;
//! - **constrained windows are the same bits however they are cut**, `Windows::new` is
//!   `Windows::with_dynamics` unconstrained to the bit, and the measurement's own sampler is the
//!   windows' to the bit.

// The sums over atoms and axes are written as the formulas' index loops.
#![allow(clippy::needless_range_loop)]

mod hydration;
mod protein;

use hydration::*;
use pantometry_forcefield::ewald::COULOMB;
use pantometry_forcefield::free_energy::{bennett_chain, window_seed};
use pantometry_forcefield::uff::{self, KCAL_PER_MOL};
use pantometry_forcefield::water::Settle;
use pantometry_forcefield::{
    Alchemical, AtLambda, Bath, Binding, Component, Coupling, Estimate, EwaldParameters,
    ForceField, Lambda, MolecularDynamics, PeriodicBox, PeriodicDecoupling, PeriodicEnergy,
    PeriodicForceField, Protocol, Quadrature, Shake, Windows,
};
use pantometry_units::BOLTZMANN;

const EPS: f64 = f64::EPSILON;
/// FreeSolv's temperature, and the experiment's.
const KELVIN: f64 = 298.15;
/// The Wigner constant of the simple cubic lattice, as `tests/the_ewald_sum_against_closed_forms.rs`
/// checks it.
const XI: f64 = -2.837_297_479_480_62;

fn kcal(j: f64) -> f64 {
    j / KCAL_PER_MOL
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn len(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}

/// The magnitudes the energy is a sum of, for its rounding.
fn magnitude(e: &PeriodicEnergy) -> f64 {
    let w = &e.ewald;
    e.bond.abs()
        + e.angle.abs()
        + e.torsion.abs()
        + e.inversion.abs()
        + e.van_der_waals.abs()
        + e.dispersion_correction.abs()
        + w.real.abs()
        + w.reciprocal.abs()
        + w.self_energy.abs()
        + w.excluded.abs()
        + w.background.abs()
}

/// A central difference and the tolerance it earns: three times the truncation term is
/// `|D(2h) − D(h)|`, and the energy's rounding, at most `64 ε` of the sum of its parts'
/// magnitudes, divided by `h`. W1's.
fn finite_difference(energy: impl Fn(f64) -> (f64, f64), h: f64) -> (f64, f64) {
    let d = |s: f64| {
        let (plus, scale) = energy(s);
        let (minus, _) = energy(-s);
        ((plus - minus) / (2.0 * s), scale)
    };
    let (d1, scale) = d(h);
    let (d2, _) = d(2.0 * h);
    (d1, 2.0 * (d2 - d1).abs() / 3.0 + 64.0 * EPS * scale / h)
}

/// Benzene in 27 lattice waters less those it overlaps, 9.3 Å across, `r_c` = 4.6 Å, δ = 10⁻⁵.
fn small_box(seed: u64) -> (Benzene, Solvated) {
    let b = Benzene::new(None);
    let s = Solvated::new(&b, 3, seed, 2.3 * ANGSTROM, 4.6 * ANGSTROM, 1e-5);
    (b, s)
}

/// The small box after `steps` of constrained BAOAB at 2 fs at `lambda`.
fn moved(b: &Benzene, s: &Solvated, lambda: Lambda, steps: usize) -> Vec<[f64; 3]> {
    let d = s.decoupling();
    let mut at = s.at.clone();
    let mut md = s
        .dynamics(b)
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 5e12,
            seed: 17,
        })
        .thermalised(&at, KELVIN, 17);
    md.run(
        &AtLambda {
            hamiltonian: &d,
            lambda,
        },
        &mut at,
        2.0 * FS,
        steps,
    );
    at
}

/// Two windows' samples are the same bits: every field, the `NaN` of a missing neighbour
/// included, which `==` would call unequal to itself.
fn same_samples(a: &[pantometry_forcefield::Sample], b: &[pantometry_forcefield::Sample]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert_eq!(x.step, y.step);
        assert_eq!(x.gradient.map(f64::to_bits), y.gradient.map(f64::to_bits));
        assert_eq!(x.to_previous.to_bits(), y.to_previous.to_bits());
        assert_eq!(x.to_next.to_bits(), y.to_next.to_bits());
    }
}

/// The largest force component, for a scale.
fn largest_force(f: &[[f64; 3]]) -> f64 {
    f.iter().flatten().fold(0.0f64, |m, x| m.max(x.abs()))
}

/// **The end states are the systems they say**, at a configuration constrained dynamics left.
/// Fully coupled: the solvated force field less benzene's interaction with its own images, and
/// with its own pairs in vacuum — the solvated field, less benzene alone in the same box, plus
/// benzene alone in vacuum. Fully decoupled: the water alone in its box plus benzene alone in
/// vacuum, the energy and every force, each atom's from the system it belongs to. Each within
/// `10⁻¹²` of the parts' magnitudes: the same terms, summed in other orders.
#[test]
fn the_end_states_are_the_systems_they_say() {
    let (b, s) = small_box(0xE5D);
    let at = moved(&b, &s, Lambda::new(0.0, 0.5, 0.7), 20);
    let d = s.decoupling();
    let n = at.len();
    let mut f = vec![[0.0; 3]; n];
    let coupled = d.energy_and_forces(&at, Lambda::COUPLED, &mut f);
    let whole = s.field.energy(&at);
    let alone_in_the_box = PeriodicForceField::new(
        &b.component,
        &b.types,
        s.cell(),
        s.field.ewald().parameters(),
    )
    .unwrap()
    .with_charges(b.charges.clone())
    .energy(&at[..12]);
    let vacuum = b.vacuum.evaluate(&at[..12]);
    let expected = whole.total - alone_in_the_box.total + vacuum.energy.total;
    let scale = magnitude(&whole);
    println!(
        "coupled {:.10e} J against {expected:.10e}; benzene's own images were {:+.4e} kcal/mol",
        coupled,
        kcal(alone_in_the_box.total - vacuum.energy.total)
    );
    assert!((coupled - expected).abs() <= 1e-12 * scale);

    // The long-range correction's cross part, written out: every ordered pair of a benzene atom
    // and a water oxygen, twice, by UFF's geometric rule; a water hydrogen has no well.
    let rc = s.field.ewald().parameters().cutoff;
    let (rc3, v) = (rc * rc * rc, s.cell().volume());
    let [xo, dox] = pantometry_forcefield::water::oxygen_vdw();
    let mut tail = 0.0;
    for t in &b.types {
        let p = t.parameters();
        let (x, dd) = ((p.vdw_distance * xo).sqrt(), (p.vdw_energy * dox).sqrt());
        let x3 = x * x * x;
        tail += s.waters as f64
            * dd
            * (x3 * x3 * x3 * x3 / (9.0 * rc3 * rc3 * rc3) - 2.0 * x3 * x3 / (3.0 * rc3));
    }
    let tail = 2.0 * 2.0 * std::f64::consts::PI / v * tail;
    let theirs = d.cross_dispersion_correction();
    println!(
        "the cross correction {:+.6e} kcal/mol, written out {:+.6e}",
        kcal(theirs),
        kcal(tail)
    );
    assert!((theirs - tail).abs() <= 1e-13 * tail.abs());

    let decoupled = d.energy_and_forces(&at, Lambda::new(0.0, 0.0, 0.0), &mut f);
    let water = s.water_alone().evaluate(&at[12..]);
    let expected = water.energy.total + vacuum.energy.total;
    println!("decoupled {decoupled:.10e} J against {expected:.10e}");
    assert!((decoupled - expected).abs() <= 1e-12 * scale);
    let fscale = largest_force(&f);
    for i in 0..n {
        let want = if i < 12 {
            vacuum.forces[i]
        } else {
            water.forces[i - 12]
        };
        for k in 0..3 {
            assert!(
                (f[i][k] - want[k]).abs() <= 1e-12 * fscale,
                "atom {i} axis {k}: {:e} against {:e}",
                f[i][k],
                want[k]
            );
        }
    }
}

/// **`∂U/∂λ_e` and `∂U/∂λ_v` are the energy's derivatives**, at a configuration constrained
/// dynamics produced in the box and at four states — the charges part-way, both part-way, the van
/// der Waals part-way with the charges off, and near the decoupled end, where the soft core is
/// softest; each λ where its central difference stays inside [0, 1] —
/// against central differences of the **whole** energy, `energy_and_forces`, not of the coupling
/// alone, each to the tolerance the difference earns.
#[test]
fn the_lambda_derivatives_are_the_whole_energys() {
    let (b, s) = small_box(0xDE1);
    let d = s.decoupling();
    for (lambda, steps) in [
        (Lambda::new(0.0, 0.5, 1.0), 20),
        (Lambda::new(0.0, 0.4, 0.8), 20),
        (Lambda::new(0.0, 0.0, 0.5), 20),
        (Lambda::new(0.0, 0.0, 0.1), 40),
    ] {
        let at = moved(&b, &s, lambda, steps);
        let c = d.coupling(&at, lambda);
        let scale = magnitude(&s.field.energy(&at));
        for (k, h) in [(1usize, 1e-3), (2, 1e-4)] {
            let energy = |e: f64| {
                let mut l = lambda;
                if k == 1 {
                    l.electrostatics += e;
                } else {
                    l.van_der_waals += e;
                }
                let mut f = vec![[0.0; 3]; at.len()];
                (d.energy_and_forces(&at, l, &mut f), scale)
            };
            // A step that would leave [0, 1] is not taken: each derivative is checked where its
            // central difference is inside.
            let at_k = lambda.components()[k];
            if at_k < 2.0 * h || at_k > 1.0 - 2.0 * h {
                continue;
            }
            let (dl, tolerance) = finite_difference(energy, h);
            println!(
                "{lambda:?}: ∂U/∂λ[{k}] = {:+.10e} J, by differences {dl:+.10e} ± {tolerance:.1e}",
                c.gradient[k]
            );
            assert!((c.gradient[k] - dl).abs() <= tolerance, "{lambda:?} {k}");
            assert!(
                tolerance < 1e-3 * c.gradient[k].abs(),
                "the check can see a part in 10³"
            );
        }
    }
}

/// **Every state `couplings` returns is that state's `coupling`, to the bit**, whichever states
/// share the call; and the difference between two states' coupling energies is the difference of
/// the whole energies, to the whole energy's rounding.
#[test]
fn every_state_of_couplings_is_its_own_coupling() {
    let (b, s) = small_box(0xC0B);
    let d = s.decoupling();
    let at = moved(&b, &s, Lambda::new(0.0, 0.0, 0.6), 20);
    let states = grid();
    let all = d.couplings(&at, &states);
    let some = d.couplings(&at, &[states[7], states[2], states[7]]);
    for (k, l) in states.iter().enumerate() {
        let one = d.coupling(&at, *l);
        assert_eq!(one.energy.to_bits(), all[k].energy.to_bits(), "{l:?}");
        assert_eq!(
            one.gradient.map(f64::to_bits),
            all[k].gradient.map(f64::to_bits)
        );
    }
    assert_eq!(some[0].energy.to_bits(), all[7].energy.to_bits());
    assert_eq!(some[1].energy.to_bits(), all[2].energy.to_bits());
    let scale = magnitude(&s.field.energy(&at));
    let mut f = vec![[0.0; 3]; at.len()];
    let u = |l: Lambda, f: &mut [[f64; 3]]| d.energy_and_forces(&at, l, f);
    let base = u(states[0], &mut f);
    for (k, l) in states.iter().enumerate().skip(1) {
        let whole = u(*l, &mut f) - base;
        let part = all[k].energy - all[0].energy;
        assert!(
            (whole - part).abs() <= 64.0 * EPS * scale,
            "{l:?}: {whole:e} against {part:e}"
        );
    }
}

/// **A solute takes out exactly the waters it overlaps**: every water with an atom within the
/// clearance of a benzene atom by minimum image, and no other, found here by brute force over the
/// 27 nearest images of each pair; the survivors keep their order and their bits.
#[test]
fn a_solute_takes_out_exactly_the_waters_it_overlaps() {
    let b = Benzene::new(None);
    let lattice = pantometry_forcefield::water::WaterBox::lattice(3, 0x0E1);
    let side = lattice.cell().lengths()[0];
    // Off centre, so that some of benzene's neighbours are across a face.
    let solute = b.centred_at([0.1 * side, 0.5 * side, 0.95 * side]);
    let clearance = 2.0 * ANGSTROM;
    let kept = lattice.without_overlaps(&solute, clearance);
    let mut want = Vec::new();
    for k in 0..lattice.count() {
        let w = &lattice.positions()[3 * k..3 * k + 3];
        let mut near = false;
        for p in w {
            for s in &solute {
                for i in -1..=1 {
                    for j in -1..=1 {
                        for k in -1..=1 {
                            let image = [
                                s[0] + f64::from(i) * side,
                                s[1] + f64::from(j) * side,
                                s[2] + f64::from(k) * side,
                            ];
                            near |= len(sub(*p, image)) < clearance;
                        }
                    }
                }
            }
        }
        if !near {
            want.extend_from_slice(w);
        }
    }
    assert_eq!(kept.positions(), &want[..]);
    assert!(kept.count() < lattice.count() && kept.count() > 0);
    println!("{} of {} waters kept", kept.count(), lattice.count());
}

/// **A restraint in the box is scaled by its own λ**: benzene restrained to two waters by a
/// Boresch restraint taken at the start and evaluated after constrained dynamics, where it is not
/// zero. `∂U/∂λ_r` is the restraint's energy, against a central difference of the whole energy in
/// `λ_r` and against `Boresch::energy` to the bit; every state of `couplings` with the restraint is
/// its own `coupling` to the bit; and the coupling energies' differences across `λ_r`, `λ_e` and
/// `λ_v` are the whole energy's. The branch W4's complex leg needs, run.
#[test]
fn a_restraint_in_the_box_is_scaled_by_its_lambda() {
    let (b, s) = small_box(0xB05);
    let o = s.field.rigid_waters();
    let restraint = pantometry_forcefield::Boresch::at(
        &s.at,
        [o[0][0], o[0][1], o[1][0]],
        [0, 1, 2],
        [
            20.0 * KCAL_PER_MOL / (ANGSTROM * ANGSTROM),
            20.0 * KCAL_PER_MOL,
            20.0 * KCAL_PER_MOL,
            20.0 * KCAL_PER_MOL,
            20.0 * KCAL_PER_MOL,
            20.0 * KCAL_PER_MOL,
        ],
    );
    let d = s.decoupling().with_restraint(restraint);
    let at = moved(&b, &s, Lambda::new(0.0, 0.5, 1.0), 20);
    let u_b = restraint.energy(&at);
    assert!(u_b > 1e-3 * KCAL_PER_MOL, "the restraint is stretched");
    let states = [
        Lambda::new(0.0, 1.0, 1.0),
        Lambda::new(0.3, 1.0, 1.0),
        Lambda::new(1.0, 0.5, 1.0),
        Lambda::new(1.0, 0.0, 0.4),
        Lambda::new(0.7, 0.2, 0.9),
    ];
    let all = d.couplings(&at, &states);
    let scale = magnitude(&s.field.energy(&at));
    let mut f = vec![[0.0; 3]; at.len()];
    let base = d.energy_and_forces(&at, states[0], &mut f);
    for (k, l) in states.iter().enumerate() {
        let one = d.coupling(&at, *l);
        assert_eq!(one.energy.to_bits(), all[k].energy.to_bits(), "{l:?}");
        assert_eq!(
            one.gradient.map(f64::to_bits),
            all[k].gradient.map(f64::to_bits)
        );
        assert_eq!(
            all[k].gradient[0].to_bits(),
            u_b.to_bits(),
            "∂U/∂λ_r is U_B"
        );
        let whole = d.energy_and_forces(&at, *l, &mut f) - base;
        let part = all[k].energy - all[0].energy;
        assert!(
            (whole - part).abs() <= 64.0 * EPS * scale,
            "{l:?}: {whole:e} against {part:e}"
        );
    }
    let lambda = Lambda::new(0.5, 0.5, 1.0);
    let energy = |e: f64| {
        let mut l = lambda;
        l.restraint += e;
        (
            d.energy_and_forces(&at, l, &mut vec![[0.0; 3]; at.len()]),
            scale,
        )
    };
    let (dl, tolerance) = finite_difference(energy, 1e-3);
    let c = d.coupling(&at, lambda);
    println!(
        "∂U/∂λ_r = {:+.10e} J, by differences {dl:+.10e} ± {tolerance:.1e}; U_B = {:.4} kcal/mol",
        c.gradient[0],
        kcal(u_b)
    );
    assert!((c.gradient[0] - dl).abs() <= tolerance);
    assert!(tolerance < 1e-3 * u_b);
}

// ---------------------------------------------------------------------------------------------
// A charge held at a distance: the closed form

/// One chloride-typed ion and one water whose oxygen carries the other charge and whose
/// hydrogens carry none, in a cube `side` across with `r_c` = `side/2`: the ion is the group, the
/// water is the rest, frozen, and the ion is held at `d` from the oxygen by SHAKE. A hydrogen has
/// no van der Waals term and no charge, so the ion's only interaction is with the oxygen, at the
/// fixed distance: its van der Waals is a constant, and its Coulomb energy depends on the
/// direction only through the periodic images.
struct HeldCharge {
    decoupling: PeriodicDecoupling,
    start: Vec<[f64; 3]>,
    dynamics: MolecularDynamics,
    centre: [f64; 3],
    d: f64,
    q: f64,
    big_q: f64,
    side: f64,
}

fn held_charge(accuracy: f64) -> HeldCharge {
    let side = 9.0 * ANGSTROM;
    let d = 3.3 * ANGSTROM;
    let (q, big_q) = (1.0, -1.0);
    let ion = Component::from_ccd(
        "data_ION\n_chem_comp.id ION\nloop_\n_chem_comp_atom.comp_id\n_chem_comp_atom.atom_id\n\
         _chem_comp_atom.type_symbol\n_chem_comp_atom.charge\n_chem_comp_atom.pdbx_aromatic_flag\n\
         _chem_comp_atom.pdbx_model_Cartn_x_ideal\n_chem_comp_atom.pdbx_model_Cartn_y_ideal\n\
         _chem_comp_atom.pdbx_model_Cartn_z_ideal\nION CL Cl 0 N 0.0 0.0 0.0\n",
    )
    .expect("one ion");
    let cell = PeriodicBox::cubic(side);
    let p = EwaldParameters::for_accuracy(&cell, 0.5 * side, accuracy);
    let field = PeriodicForceField::solvated(&ion, &uff::assign(&ion), 1, cell, p)
        .unwrap()
        .with_charges(vec![q, big_q, 0.0, 0.0]);
    let centre = [0.5 * side; 3];
    let body = pantometry_forcefield::water::body_frame();
    let mut start = vec![[centre[0] + d, centre[1], centre[2]]];
    for p in body {
        start.push([centre[0] + p[0], centre[1] + p[1], centre[2] + p[2]]);
    }
    // The oxygen is not at the water's centre of mass: put it at the centre instead.
    let shift = sub(centre, start[1]);
    for p in &mut start[1..] {
        *p = [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]];
    }
    let decoupling = PeriodicDecoupling::new(&field, &[true, false, false, false]).unwrap();
    let h = pantometry_forcefield::Element::H.mass();
    let dynamics = MolecularDynamics::new(vec![h, 16.0 * h, h, h])
        .with_frozen(vec![false, true, true, true])
        .with_bond_constraints(Shake::new(vec![[0, 1]], vec![d]));
    HeldCharge {
        decoupling,
        start,
        dynamics,
        centre,
        d,
        q,
        big_q,
        side,
    }
}

/// Gauss–Legendre nodes and weights on [−1, 1], by Newton's method on `P_n` from the usual
/// first guesses.
fn gauss_legendre(n: usize) -> Vec<(f64, f64)> {
    (0..n)
        .map(|i| {
            let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
            let mut dp = 0.0;
            for _ in 0..100 {
                let (mut p0, mut p1) = (1.0, x);
                for k in 2..=n {
                    let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                    p0 = p1;
                    p1 = p2;
                }
                dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
                let step = p1 / dp;
                x -= step;
                if step.abs() < 1e-16 {
                    break;
                }
            }
            (x, 2.0 / ((1.0 - x * x) * dp * dp))
        })
        .collect()
}

/// The cross Coulomb energy at full strength with the ion at each direction of a Gauss–Legendre ×
/// trapezoid product rule on the sphere, and each direction's weight (summing to 1).
fn on_the_sphere(h: &HeldCharge, n: usize) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    let m = 2 * n;
    for (mu, w) in gauss_legendre(n) {
        let s = (1.0 - mu * mu).sqrt();
        for j in 0..m {
            let phi = 2.0 * std::f64::consts::PI * j as f64 / m as f64;
            let mut at = h.start.clone();
            at[0] = [
                h.centre[0] + h.d * s * phi.cos(),
                h.centre[1] + h.d * s * phi.sin(),
                h.centre[2] + h.d * mu,
            ];
            let e = h.decoupling.coupling(&at, Lambda::COUPLED).gradient[1];
            out.push((e, 0.5 * w / m as f64));
        }
    }
    out
}

/// **The sphere average of the cross Ewald energy is a closed form.** The periodic potential of a
/// unit charge with its neutralising background, `ψ`, satisfies `∇²ψ = −4π (δ − 1/V)` inside the
/// ball of radius `L` about the charge, which holds no image. So `ψ − 1/r − 2πr²/3V` is harmonic
/// there, and its average over a sphere of radius `d < L` is its value at the centre — the Wigner
/// term `ξ/L` that makes a lone charge's energy `q²ξ/2L`. Hence
/// `⟨E⟩ = k_e q Q (1/d + ξ/L + 2πd²/3L³)` over directions, whatever the charges' images do to
/// each direction alone. By a 32-point Gauss–Legendre rule in `cos θ` and 64 in `φ` at δ = 10⁻¹⁰,
/// to four times Kolafa and Perram's estimate of the sum's RMS error; the spread over directions
/// is printed beside it, and required, to show the average is not of a constant.
#[test]
fn a_charge_on_a_sphere_averages_to_its_closed_form() {
    let h = held_charge(1e-10);
    let points = on_the_sphere(&h, 32);
    let mean: f64 = points.iter().map(|(e, w)| e * w).sum();
    let total: f64 = points.iter().map(|(_, w)| w).sum();
    let (lo, hi) = points
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), (e, _)| {
            (a.min(*e), b.max(*e))
        });
    let l = h.side;
    let exact = COULOMB
        * h.q
        * h.big_q
        * (1.0 / h.d + XI / l + 2.0 * std::f64::consts::PI * h.d * h.d / (3.0 * l * l * l));
    println!(
        "⟨E⟩ over the sphere {:.12} kcal/mol against the closed form {:.12}; directions span \
         {:.4} to {:.4}",
        kcal(mean),
        kcal(exact),
        kcal(lo),
        kcal(hi)
    );
    assert!((total - 1.0).abs() < 1e-14);
    // Four times Kolafa and Perram's RMS estimate of the sum's error for these two charges: the
    // quadrature's own error is far below it (a smooth integrand on a 32 × 64 product rule).
    let cell = h.decoupling.field().cell();
    let bound = 4.0
        * h.decoupling
            .field()
            .ewald()
            .parameters()
            .energy_error(&cell, 2.0);
    println!(
        "off by {:.3e} kcal/mol against a bound of {:.3e}",
        kcal(mean - exact),
        kcal(bound)
    );
    assert!((mean - exact).abs() <= bound);
    assert!(hi - lo > 1e-3 * exact.abs(), "the directions differ");
}

/// **The free energy of decoupling a held charge is the quadrature's.** The ion moves on the sphere
/// of radius `d` about the frozen oxygen; every other term is a constant there, so at coupling
/// `λ_e` its directions are distributed as `exp(−λ_e E(n)/k_BT)` and
/// `F(0) − F(1) = k_BT ln ⟨exp(−E/k_BT)⟩` over directions, uniform, by the product rule above. The
/// windows — five states, the ion held by SHAKE to a frozen atom — must give it by BAR and by TI,
/// each within four of its own standard errors; TI against the trapezoid of the exact `⟨E⟩_λ` at
/// the same five states, which is its own bias. And the decoupled window's mean `∂U/∂λ_e` is
/// the uniform average, within four standard errors. δ = 10⁻⁶ for the windows' speed, the
/// quadrature at the same δ, so that it is the same Hamiltonian.
#[test]
fn a_held_charge_decouples_to_its_quadrature() {
    let h = held_charge(1e-6);
    let kt = BOLTZMANN.to_si() * KELVIN;
    let points = on_the_sphere(&h, 24);
    // ⟨·⟩_λ over directions with weight exp(−λE/kT), shifted for range.
    let e0 = points.iter().map(|(e, w)| e * w).sum::<f64>();
    let average_at = |l: f64| {
        let (mut z, mut ez) = (0.0, 0.0);
        for (e, w) in &points {
            let b = w * (-l * (e - e0) / kt).exp();
            z += b;
            ez += b * e;
        }
        (z, ez / z)
    };
    let exact = kt * average_at(1.0).0.ln() - e0;
    let schedule: Vec<Lambda> = [1.0, 0.75, 0.5, 0.25, 0.0]
        .iter()
        .map(|&e| Lambda::new(0.0, e, 1.0))
        .collect();
    let mut trapezoid = 0.0;
    for k in 0..schedule.len() - 1 {
        let (a, b) = (schedule[k].electrostatics, schedule[k + 1].electrostatics);
        trapezoid += 0.5 * (b - a) * (average_at(a).1 + average_at(b).1);
    }
    let protocol = Protocol {
        time_step: 2.0 * FS,
        temperature: KELVIN,
        friction: 5e12,
        equilibration: 500,
        stride: 10,
        samples: 800,
        seed: 0x4E1D,
    };
    let mut w = Windows::with_dynamics(&h.decoupling, schedule, &h.start, &h.dynamics, protocol);
    assert!(w.run(&h.decoupling, u64::MAX));
    let bar = w.bennett_total();
    let ti = w.thermodynamic_integration(Quadrature::Trapezoid).unwrap();
    let last = w.windows().last().unwrap();
    let grad: Vec<f64> = last.samples().iter().map(|s| s.gradient[1]).collect();
    let g = Estimate::of(&grad);
    println!(
        "decoupling: BAR {:+.4} ± {:.4}, TI {:+.4} ± {:.4} kcal/mol; exact {:+.4}, the trapezoid \
         of the exact ⟨E⟩ {:+.4}; ⟨∂U/∂λ_e⟩ decoupled {:+.4} ± {:.4} against the uniform {:+.4}; \
         −⟨E⟩₀ alone would be {:+.4}",
        kcal(bar.value),
        kcal(bar.error),
        kcal(ti.value),
        kcal(ti.error),
        kcal(exact),
        kcal(trapezoid),
        kcal(g.mean),
        kcal(g.error),
        kcal(e0),
        kcal(-e0)
    );
    // The ion stayed at its distance; the oxygen where it was.
    let at = last.positions();
    assert_eq!(at[1], h.start[1]);
    assert!((len(sub(at[0], at[1])) - h.d).abs() <= 16.0 * EPS * 6.0 * ANGSTROM);
    assert!(((bar.value - exact) / bar.error).abs() < 4.0, "BAR");
    assert!(((ti.value - trapezoid) / ti.error).abs() < 4.0, "TI");
    assert!(
        ((g.mean - e0) / g.error).abs() < 4.0,
        "the decoupled average"
    );
    // The Boltzmann weighting is visible: BAR is held within 4σ of the exact answer, so for the
    // check to exclude the unweighted answer −⟨E⟩₀ the two must be more than 8σ apart, the two
    // four-σ windows not touching. Measured 0.764 kcal/mol against 8σ = 0.564: a margin of 1.35.
    assert!((exact + e0).abs() > 8.0 * bar.error, "the test can tell");
}

// ---------------------------------------------------------------------------------------------
// The windows with constraints

/// Constrained windows in the box — SETTLE on the water, SHAKE on benzene, `PeriodicDecoupling`
/// — are the same bits however they are advanced: in one call, in single steps, or unevenly.
#[test]
fn constrained_windows_are_the_same_bits_however_they_are_cut() {
    let (b, s) = small_box(0xC17);
    let d = s.decoupling();
    let schedule = vec![
        Lambda::COUPLED,
        Lambda::new(0.0, 0.0, 1.0),
        Lambda::new(0.0, 0.0, 0.3),
    ];
    let protocol = Protocol {
        time_step: 2.0 * FS,
        temperature: KELVIN,
        friction: 5e12,
        equilibration: 4,
        stride: 3,
        samples: 4,
        seed: 0xC17,
    };
    let fresh = || Windows::with_dynamics(&d, schedule.clone(), &s.at, &s.dynamics(&b), protocol);
    let mut whole = fresh();
    assert!(whole.run(&d, u64::MAX));
    for cuts in [vec![1u64; 16], vec![7, 9], vec![2, 11, 3]] {
        let mut part = fresh();
        for c in cuts {
            part.run(&d, c);
        }
        assert!(part.is_complete());
        for (x, y) in part.windows().iter().zip(whole.windows()) {
            same_samples(x.samples(), y.samples());
            assert_eq!(x.positions(), y.positions());
        }
    }
    // The windows ran with the constraints: benzene's C–H at their length.
    let shake = b.shake();
    for w in whole.windows() {
        let at = w.positions();
        for (k, &[i, j]) in shake.bonds().iter().enumerate() {
            assert!(
                (len(sub(at[i], at[j])) - shake.lengths()[k]).abs() <= 16.0 * EPS * 10.0 * ANGSTROM
            );
        }
        assert!(w.dynamics().bond_constraints().is_some());
        assert!(w.dynamics().constraints().is_some());
    }
}

/// `Windows::new` is `Windows::with_dynamics` on an unconstrained template with the same masses and
/// frozen atoms, to the bit.
#[test]
fn windows_new_is_with_dynamics_unconstrained() {
    let h = held_charge(1e-5);
    let schedule = vec![Lambda::COUPLED, Lambda::new(0.0, 0.0, 1.0)];
    let protocol = Protocol {
        time_step: 1.0 * FS,
        temperature: KELVIN,
        friction: 5e12,
        equilibration: 3,
        stride: 2,
        samples: 5,
        seed: 9,
    };
    let masses = h.dynamics.masses().to_vec();
    let frozen = h.dynamics.frozen().to_vec();
    let mut a = Windows::new(
        &h.decoupling,
        schedule.clone(),
        &h.start,
        masses.clone(),
        frozen.clone(),
        protocol,
    );
    let mut b = Windows::with_dynamics(
        &h.decoupling,
        schedule,
        &h.start,
        &MolecularDynamics::new(masses).with_frozen(frozen),
        protocol,
    );
    a.run(&h.decoupling, u64::MAX);
    b.run(&h.decoupling, u64::MAX);
    for (x, y) in a.windows().iter().zip(b.windows()) {
        same_samples(x.samples(), y.samples());
        assert_eq!(x.positions(), y.positions());
    }
}

/// **The measurement's sampler is the windows' own**: [`sample_window`], which the hydration
/// measurement runs so that each window records every candidate state, gives at a window's own
/// state and its neighbours the bits `Windows` records — the gradient and both energy differences
/// of every sample — and ends where it does.
#[test]
fn the_measurements_sampler_is_the_windows_own() {
    let (b, s) = small_box(0x5A3);
    let d = s.decoupling();
    let g = grid();
    let picks = [0usize, 4, 12];
    let schedule: Vec<Lambda> = picks.iter().map(|&k| g[k]).collect();
    let protocol = Protocol {
        time_step: 2.0 * FS,
        temperature: KELVIN,
        friction: 5e12,
        equilibration: 3,
        stride: 2,
        samples: 3,
        seed: 0x5A3,
    };
    let template = s.dynamics(&b);
    let mut w = Windows::with_dynamics(&d, schedule.clone(), &s.at, &template, protocol);
    w.run(&d, u64::MAX);
    for (pos, window) in w.windows().iter().enumerate() {
        let mut rows = Vec::new();
        let end = sample_window(
            &d,
            &s.at,
            &template,
            protocol,
            pos,
            schedule[pos],
            &g,
            |_, c| rows.push(c.to_vec()),
        );
        assert_eq!(end, window.positions());
        for (row, sample) in rows.iter().zip(window.samples()) {
            let own = &row[0];
            assert_eq!(own.gradient, sample.gradient);
            let at = |k: usize| row[1 + picks[k]].energy - own.energy;
            if pos > 0 {
                assert_eq!(at(pos - 1).to_bits(), sample.to_previous.to_bits());
            }
            if pos + 1 < picks.len() {
                assert_eq!(at(pos + 1).to_bits(), sample.to_next.to_bits());
            }
        }
        assert_eq!(rows.len(), window.samples().len());
    }
}

// ---------------------------------------------------------------------------------------------
// The measurement

/// The candidate states, in path order: the charges off at full van der Waals, then the van der
/// Waals off with no charge. Every window of the measurement is one of them and records its
/// energy at all of them, so a window inserted between two that ran needs neither run again.
fn grid() -> Vec<Lambda> {
    let mut g = Vec::new();
    for e in [1.0, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0.0] {
        g.push(Lambda::new(0.0, e, 1.0));
    }
    for k in 1..=20 {
        // 0.95, 0.9, …, 0.05, 0: twenty steps of 0.05, written as integers over 20 so that each
        // state is the same double every time.
        g.push(Lambda::new(0.0, 0.0, f64::from(20 - k) / 20.0));
    }
    g
}

/// One window of the measurement, as [`Windows`] would run it at schedule position `index`: the
/// template's dynamics in the protocol's bath, thermalised from `window_seed(protocol.seed,
/// index)`, at `lambda`, from `start`. After the equilibration, every `stride` steps, `on` is
/// given the step and the couplings at `lambda` followed by every one of `states`. Returns where
/// it ended.
#[allow(clippy::too_many_arguments)]
fn sample_window(
    h: &PeriodicDecoupling,
    start: &[[f64; 3]],
    template: &MolecularDynamics,
    p: Protocol,
    index: usize,
    lambda: Lambda,
    states: &[Lambda],
    mut on: impl FnMut(u64, &[Coupling]),
) -> Vec<[f64; 3]> {
    let seed = window_seed(p.seed, index);
    let mut md = template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: p.temperature,
            friction: p.friction,
            seed,
        })
        .thermalised(start, p.temperature, seed);
    let potential = AtLambda {
        hamiltonian: h,
        lambda,
    };
    let mut at = start.to_vec();
    md.prepare(&potential, &at);
    let mut lambdas = vec![lambda];
    lambdas.extend_from_slice(states);
    let mut taken = 0;
    // Clippy on current stable suggests `u64::is_multiple_of`, stabilised in 1.87; this crate
    // builds on 1.78.
    #[allow(clippy::manual_is_multiple_of)]
    while taken < p.samples {
        md.step(&potential, &mut at, p.time_step);
        let s = md.steps();
        if s > p.equilibration && (s - p.equilibration) % p.stride == 0 {
            on(s, &h.couplings(&at, &lambdas));
            taken += 1;
        }
    }
    at
}

/// Where the measurement writes: `PANTOMETRY_W3_DIR`, or `pantometry-w3` in the system's temporary
/// directory. Never the repository.
fn results_dir() -> std::path::PathBuf {
    let dir = std::env::var_os("PANTOMETRY_W3_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("pantometry-w3"));
    std::fs::create_dir_all(&dir).expect("the results directory");
    dir
}

/// A lock that makes a second run in the same directory refuse, rather than write into the files
/// the first is writing. Removed when dropped; a run that was killed leaves it, and says so.
struct Lock(std::path::PathBuf);

impl Lock {
    fn take(dir: &std::path::Path, name: &str) -> Lock {
        let path = dir.join(format!("{name}.lock"));
        std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap_or_else(|e| {
                panic!(
                    "{} exists ({e}): another run is writing here, or one was killed — remove it \
                     by hand once no run is",
                    path.display()
                )
            });
        Lock(path)
    }
}

impl Drop for Lock {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

/// Appends a line to the log and prints it.
fn log(dir: &std::path::Path, line: &str) {
    use std::io::Write;
    println!("{line}");
    let mut f = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(dir.join("hydration.log"))
        .expect("the log");
    writeln!(f, "{line}").expect("the log");
}

/// The measurement's protocol: 2 fs, a 1 ps⁻¹ bath at 298.15 K, 20 ps discarded and then 2000
/// samples, one every 100 fs, 200 ps a window. `the_cost_measured` put a window at 0.79 h before
/// any ran.
const PROTOCOL: Protocol = Protocol {
    time_step: 2e-15,
    temperature: KELVIN,
    friction: 1e12,
    equilibration: 10_000,
    stride: 50,
    samples: 2000,
    seed: 0x3_3A7E,
};

/// The first schedule, as indices into [`grid`]: the charges off in two steps and the van der
/// Waals in ten, every other candidate. Thirteen windows.
const SCHEDULE: [usize; 13] = [0, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28];

/// An interval whose overlap is below this gets a window at the candidate between its ends: 3c's
/// rule and threshold, three times the 0.03 Klimovich, Shirts and Mobley find tolerable.
const OVERLAP_THRESHOLD: f64 = 0.1;

/// What one window recorded, read back from its file.
#[derive(Clone, Debug)]
struct Recorded {
    grid: usize,
    lambda: Lambda,
    gradients: Vec<[f64; 3]>,
    /// Per sample, the coupling energy at every candidate, joules.
    energies: Vec<Vec<f64>>,
    seconds: f64,
}

impl Recorded {
    /// Reduced `u(to) − u(own)` per sample.
    fn delta(&self, to: usize, kt: f64) -> Vec<f64> {
        self.energies
            .iter()
            .map(|e| (e[to] - e[self.grid]) / kt)
            .collect()
    }

    fn slice(&self, range: std::ops::Range<usize>) -> Recorded {
        Recorded {
            gradients: self.gradients[range.clone()].to_vec(),
            energies: self.energies[range].to_vec(),
            ..self.clone()
        }
    }
}

/// BAR between neighbouring windows, reduced, the delta method's variance: the forward set
/// `u_b − u_a` on `a`'s samples and the reverse set the same difference on `b`'s.
fn interval_bar(a: &Recorded, b: &Recorded, kt: f64) -> pantometry_forcefield::Bar {
    let reverse: Vec<f64> = b.delta(a.grid, kt).iter().map(|x| -x).collect();
    pantometry_forcefield::free_energy::bar_correlated(&a.delta(b.grid, kt), &reverse)
}

/// BAR along `records`, reduced, by the delta method window by window.
fn chain(records: &[Recorded], kt: f64) -> pantometry_forcefield::FreeEnergy {
    let n = records.len();
    let windows: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|j| {
            let prev = if j > 0 {
                records[j].delta(records[j - 1].grid, kt)
            } else {
                Vec::new()
            };
            let next = if j + 1 < n {
                records[j].delta(records[j + 1].grid, kt)
            } else {
                Vec::new()
            };
            (prev, next)
        })
        .collect();
    bennett_chain(&windows)
}

/// TI on the trapezoid along `records`' path, kcal/mol, with its error.
fn trapezoid(records: &[Recorded]) -> (f64, f64, Vec<Estimate>) {
    let n = records.len();
    let terms: Vec<Estimate> = (0..n)
        .map(|k| {
            let lo = records[k.saturating_sub(1)].lambda.components();
            let hi = records[(k + 1).min(n - 1)].lambda.components();
            let w = [0, 1, 2].map(|c| 0.5 * (hi[c] - lo[c]));
            let ys: Vec<f64> = records[k]
                .gradients
                .iter()
                .map(|g| (0..3).map(|c| g[c] * w[c]).sum())
                .collect();
            Estimate::of(&ys)
        })
        .collect();
    (
        kcal(terms.iter().map(|t| t.mean).sum()),
        kcal(terms.iter().map(|t| t.error * t.error).sum::<f64>().sqrt()),
        terms,
    )
}

/// `ln mean(exp(x))` without overflow.
fn log_mean_exp(xs: &[f64]) -> f64 {
    let m = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    m + (xs.iter().map(|x| (x - m).exp()).sum::<f64>() / xs.len() as f64).ln()
}

/// The measurement's system and its windows' common inputs.
struct Measurement {
    benzene: Benzene,
    solvated: Solvated,
    decoupling: PeriodicDecoupling,
    template: MolecularDynamics,
}

/// Benzene's QEq charges as 3c-3's solvent leg had them — QEq on the ligand of 181L at the crystal
/// pose — after checking that 181L's ligand is the dictionary's benzene atom for atom: the same
/// elements in the same order, and the same bonds.
fn crystal_charges() -> Vec<f64> {
    let system = protein::system();
    let b = Binding::new(&system, 3.0 * ANGSTROM).expect("181L at 3 Å");
    let bnz = Component::from_ccd(BNZ).unwrap();
    let elements: Vec<_> = bnz.atoms().iter().map(|a| a.element).collect();
    assert_eq!(&b.elements()[b.ligand_range()], &elements[..]);
    let ours = ForceField::new(&bnz, &uff::assign(&bnz)).unwrap();
    let theirs: Vec<[usize; 2]> = b
        .ligand_force_field()
        .stretches()
        .iter()
        .map(|s| s.atoms)
        .collect();
    let mine: Vec<[usize; 2]> = ours.stretches().iter().map(|s| s.atoms).collect();
    assert_eq!(
        theirs, mine,
        "181L's benzene is the dictionary's, bond for bond"
    );
    b.ligand_force_field().charges().to_vec()
}

/// The density the box is built at: 33.00 waters per nm³, 0.9872 g/cm³, where Yeh and Hummer
/// (*J. Phys. Chem. B* **108**, 15873 (2004), read in W2) measure TIP3P's pressure under Ewald at
/// −2.9 ± 2.5 bar. W2's liquid at 0.997 g/cm³ read +139 to +180 bar; a hydration free energy at
/// that pressure would carry the solute's partial molar volume times it, a few tenths of a kcal/mol,
/// so the box is at the model's own density at one atmosphere instead of experiment's.
fn density() -> f64 {
    33.0e27 * pantometry_forcefield::water::molecular_mass()
}

/// 512 lattice waters at [`density`] (24.94 Å, the smallest lattice box past `2 r_c` plus
/// benzene's 5.0 Å span at `r_c` = 9 Å), less every one with an atom within 1.4 Å of a benzene
/// atom — five, about benzene's volume in waters (liquid benzene's 89 cm³/mol is 4.9 of water's
/// 18.07) — and Ewald at δ = 10⁻⁵.
fn measurement() -> Measurement {
    let benzene = Benzene::new(Some(crystal_charges()));
    let solvated = Solvated::at_density(
        &benzene,
        8,
        density(),
        0x3_3A7E,
        1.4 * ANGSTROM,
        9.0 * ANGSTROM,
        1e-5,
    );
    let decoupling = solvated.decoupling();
    let template = solvated.dynamics(&benzene);
    Measurement {
        benzene,
        solvated,
        decoupling,
        template,
    }
}

/// Writes positions one atom a line, exactly (`{:e}` round-trips).
fn write_positions(path: &std::path::Path, at: &[[f64; 3]]) {
    let text: String = at
        .iter()
        .map(|p| format!("{:e} {:e} {:e}\n", p[0], p[1], p[2]))
        .collect();
    std::fs::write(path, text).expect("positions");
}

fn read_positions(path: &std::path::Path) -> Option<Vec<[f64; 3]>> {
    let text = std::fs::read_to_string(path).ok()?;
    Some(
        text.lines()
            .map(|l| {
                let v: Vec<f64> = l.split_whitespace().map(|x| x.parse().unwrap()).collect();
                [v[0], v[1], v[2]]
            })
            .collect(),
    )
}

/// The coupled box after melting the lattice with benzene grown into it: 0.5 ps at 0.5 fs in a
/// 50 ps⁻¹ bath at each of `λ_v` = 0.5 with no charges — the soft core finite where a lattice
/// water overlaps benzene — `λ_v` = 1 with no charges, and fully coupled; then 20 ps at 2 fs in a
/// 5 ps⁻¹ bath, fully coupled. W2's melt with the soft core in front of it; every constraint on
/// throughout. Read from `start.txt` if a run has written it.
fn equilibrated(m: &Measurement, dir: &std::path::Path) -> Vec<[f64; 3]> {
    let path = dir.join("start.txt");
    if let Some(at) = read_positions(&path) {
        assert_eq!(
            at.len(),
            m.solvated.at.len(),
            "start.txt is another system's"
        );
        return at;
    }
    let t = std::time::Instant::now();
    let potential = AtLambda {
        hamiltonian: &m.decoupling,
        lambda: Lambda::COUPLED,
    };
    let mut at = m.solvated.at.clone();
    let bath = |friction: f64, seed: u64| Bath::Langevin {
        temperature: KELVIN,
        friction,
        seed,
    };
    let mut md = m
        .template
        .clone()
        .with_bath(bath(50e12, 0xE01))
        .thermalised(&at, KELVIN, 0xE01);
    for (k, lambda) in [Lambda::new(0.0, 0.0, 0.5), Lambda::new(0.0, 0.0, 1.0)]
        .into_iter()
        .enumerate()
    {
        let growing = AtLambda {
            hamiltonian: &m.decoupling,
            lambda,
        };
        md = md.with_bath(bath(50e12, 0xE10 + k as u64));
        md.forget_forces();
        md.run(&growing, &mut at, 0.5 * FS, 1000);
    }
    md.forget_forces();
    md.run(&potential, &mut at, 0.5 * FS, 1000);
    let mut md = md.with_bath(bath(5e12, 0xE02));
    md.run(&potential, &mut at, 2.0 * FS, 10_000);
    write_positions(&path, &at);
    log(
        dir,
        &format!(
            "equilibrated in {:.0} s: 1.5 ps at 0.5 fs growing benzene in, then 20 ps at 2 fs",
            t.elapsed().as_secs_f64()
        ),
    );
    at
}

fn header(g: usize, l: Lambda, p: &Protocol, atoms: usize) -> String {
    format!(
        "# candidate {g} lambda {:e} {:e} {:e} protocol {:e} {:e} {:e} {} {} {} {} atoms {atoms} \
         states {}",
        l.restraint,
        l.electrostatics,
        l.van_der_waals,
        p.time_step,
        p.temperature,
        p.friction,
        p.equilibration,
        p.stride,
        p.samples,
        p.seed,
        grid().len()
    )
}

fn parse(text: &str) -> Recorded {
    let mut lines = text.lines();
    let head: Vec<&str> = lines.next().expect("a header").split_whitespace().collect();
    let field = |name: &str| head.iter().position(|&w| w == name).expect(name) + 1;
    let grid_index: usize = head[field("candidate")].parse().unwrap();
    let l = field("lambda");
    let lambda = Lambda::new(
        head[l].parse().unwrap(),
        head[l + 1].parse().unwrap(),
        head[l + 2].parse().unwrap(),
    );
    let states: usize = head[field("states")].parse().unwrap();
    let (mut gradients, mut energies, mut seconds) = (Vec::new(), Vec::new(), f64::NAN);
    for line in lines {
        if let Some(rest) = line.strip_prefix("# done seconds ") {
            seconds = rest.trim().parse().unwrap();
            continue;
        }
        let v: Vec<f64> = line
            .split_whitespace()
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        assert_eq!(v.len(), 3 + states, "a sample line");
        gradients.push([v[0], v[1], v[2]]);
        energies.push(v[3..].to_vec());
    }
    assert!(seconds.is_finite(), "a complete window says so");
    Recorded {
        grid: grid_index,
        lambda,
        gradients,
        energies,
        seconds,
    }
}

/// Window `g`'s record: from its file if a run completed it with this protocol, or run now —
/// written sample by sample to a `.partial` file, renamed when complete. A `.partial` left by a
/// stopped run is started again.
// Clippy on current stable suggests `usize::is_multiple_of`, stabilised in 1.87; this crate
// builds on 1.78.
#[allow(clippy::manual_is_multiple_of)]
fn window(m: &Measurement, dir: &std::path::Path, start: &[[f64; 3]], g: usize) -> Recorded {
    use std::io::Write;
    let path = dir.join(format!("window_{g:02}.txt"));
    let grid = grid();
    let want = header(g, grid[g], &PROTOCOL, start.len());
    if let Ok(text) = std::fs::read_to_string(&path) {
        assert_eq!(
            text.lines().next().unwrap_or(""),
            want,
            "{} was run differently",
            path.display()
        );
        let r = parse(&text);
        assert_eq!(r.energies.len(), PROTOCOL.samples);
        return r;
    }
    let t = std::time::Instant::now();
    let partial = path.with_extension("partial");
    let mut file = std::io::BufWriter::new(std::fs::File::create(&partial).expect("partial"));
    writeln!(file, "{want}").unwrap();
    let total = PROTOCOL.steps_per_window();
    log(
        dir,
        &format!(
            "window at candidate {g} λ = ({}, {}, {}) begins: {total} steps",
            grid[g].restraint, grid[g].electrostatics, grid[g].van_der_waals
        ),
    );
    let mut taken = 0usize;
    sample_window(
        &m.decoupling,
        start,
        &m.template,
        PROTOCOL,
        g,
        grid[g],
        &grid,
        |s, c| {
            let mut line = format!(
                "{s} {:e} {:e} {:e}",
                c[0].gradient[0], c[0].gradient[1], c[0].gradient[2]
            );
            for x in &c[1..] {
                line.push_str(&format!(" {:e}", x.energy));
            }
            writeln!(file, "{line}").unwrap();
            taken += 1;
            if taken % 150 == 0 {
                file.flush().unwrap();
                let secs = t.elapsed().as_secs_f64();
                log(
                    dir,
                    &format!(
                        "candidate {g}: {taken} of {} samples, step {s} of {total}, {secs:.0} s, \
                         {:.2} ms per step",
                        PROTOCOL.samples,
                        secs * 1e3 / s as f64
                    ),
                );
            }
        },
    );
    let seconds = t.elapsed().as_secs_f64();
    writeln!(file, "# done seconds {seconds:e}").unwrap();
    drop(file);
    std::fs::rename(&partial, &path).expect("rename");
    log(
        dir,
        &format!(
            "window at candidate {g} done in {seconds:.0} s ({:.2} ms per step)",
            seconds * 1e3 / total as f64
        ),
    );
    parse(&std::fs::read_to_string(&path).unwrap())
}

/// **Benzene's hydration free energy in TIP3P** (release, one core, ten to fifteen hours):
/// decoupling in the 507-water box, charges then van der Waals, every window recording every candidate; a window
/// inserted at the candidate between each pair whose overlap is below [`OVERLAP_THRESHOLD`], up
/// to twice; BAR with the delta method's variance, TI on the trapezoid, halves, and EXP both
/// ways. `ΔG_hyd = −ΔG_decouple`. Written to [`results_dir`]; a stopped run resumes window by
/// window. `PANTOMETRY_W3_DIR=… cargo test --release -p pantometry-forcefield --test
/// benzene_hydrated_in_tip3p the_hydration_free_energy_measured -- --ignored --nocapture`.
#[test]
#[ignore = "benzene in 507 TIP3P waters, ten to fifteen hours with --release: run with --release -- --ignored --nocapture"]
fn the_hydration_free_energy_measured() {
    let dir = results_dir();
    let _lock = Lock::take(&dir, "hydration");
    let t = std::time::Instant::now();
    let m = measurement();
    let kt = BOLTZMANN.to_si() * KELVIN;
    log(
        &dir,
        &format!(
            "benzene in {} waters, box {:.3} Å, {} wave vectors, α {:.4} Å⁻¹; charges {:?}",
            m.solvated.waters,
            m.solvated.cell().lengths()[0] / ANGSTROM,
            m.decoupling.field().ewald().wave_vectors(),
            m.decoupling.field().ewald().parameters().alpha * ANGSTROM,
            m.benzene.charges
        ),
    );
    let start = equilibrated(&m, &dir);
    let mut schedule: Vec<usize> = SCHEDULE.to_vec();
    let mut records: Vec<Recorded> = schedule
        .iter()
        .map(|&g| window(&m, &dir, &start, g))
        .collect();
    for round in 0..2 {
        let mut inserted = Vec::new();
        for k in 0..schedule.len() - 1 {
            let bar = interval_bar(&records[k], &records[k + 1], kt);
            if bar.overlap < OVERLAP_THRESHOLD {
                let between = (schedule[k] + schedule[k + 1]) / 2;
                log(
                    &dir,
                    &format!(
                        "round {round}: overlap {:.3} between candidates {} and {}: {}",
                        bar.overlap,
                        schedule[k],
                        schedule[k + 1],
                        if between > schedule[k] {
                            format!("inserting {between}")
                        } else {
                            "nothing between them".to_string()
                        }
                    ),
                );
                if between > schedule[k] {
                    inserted.push((k + 1, between));
                }
            }
        }
        if inserted.is_empty() {
            break;
        }
        for (offset, (pos, g)) in inserted.into_iter().enumerate() {
            schedule.insert(pos + offset, g);
            records.insert(pos + offset, window(&m, &dir, &start, g));
        }
    }
    analyse(&dir, &m, &records, &start, kt);
    log(
        &dir,
        &format!("the measurement took {:.0} s", t.elapsed().as_secs_f64()),
    );
}

/// Prints the windows' table, the segments, the totals, the halves and EXP, and the hydration
/// free energy beside experiment and GB.
fn analyse(
    dir: &std::path::Path,
    m: &Measurement,
    records: &[Recorded],
    start: &[[f64; 3]],
    kt: f64,
) {
    let n = records.len();
    let to_kcal = |x: f64| kcal(kt * x);
    let (ti, ti_error, terms) = trapezoid(records);
    log(
        dir,
        "| window | candidate | λ_e | λ_v | TI term (kcal/mol) | τ (samples) | wall (s) | → next: BAR (kcal/mol) | overlap | EXP fwd | EXP rev |",
    );
    log(
        dir,
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    );
    let (mut fwd, mut rev) = (0.0, 0.0);
    let mut lowest = f64::INFINITY;
    for k in 0..n {
        let r = &records[k];
        let next = if k + 1 < n {
            let s = &records[k + 1];
            let f = r.delta(s.grid, kt);
            let reverse: Vec<f64> = s.delta(r.grid, kt).iter().map(|x| -x).collect();
            let bar = interval_bar(r, s, kt);
            lowest = lowest.min(bar.overlap);
            let neg: Vec<f64> = f.iter().map(|x| -x).collect();
            let ef = -log_mean_exp(&neg);
            let er = log_mean_exp(&reverse);
            fwd += ef;
            rev += er;
            format!(
                "{:+.3} ± {:.3} | {:.3} | {:+.3} | {:+.3}",
                to_kcal(bar.delta),
                to_kcal(bar.error()),
                bar.overlap,
                to_kcal(ef),
                to_kcal(er)
            )
        } else {
            "| | |".to_string()
        };
        log(
            dir,
            &format!(
                "| {k} | {} | {} | {} | {:+.3} ± {:.3} | {:.1} | {:.0} | {next} |",
                r.grid,
                r.lambda.electrostatics,
                r.lambda.van_der_waals,
                kcal(terms[k].mean),
                kcal(terms[k].error),
                terms[k].tau,
                r.seconds
            ),
        );
    }
    let total = chain(records, kt);
    let (bar, bar_error) = (to_kcal(total.value), to_kcal(total.error));
    // The two segments: where the charges are on their way off, and after.
    let split = records
        .iter()
        .position(|r| r.lambda.electrostatics == 0.0)
        .expect("a window with the charges off");
    for (name, range) in [
        ("charges off", 0..split + 1),
        ("van der Waals off", split..n),
    ] {
        let c = chain(&records[range.clone()], kt);
        let (t, te, _) = trapezoid(&records[range]);
        log(
            dir,
            &format!(
                "{name}: BAR {:+.3} ± {:.3}, TI {t:+.3} ± {te:.3} kcal/mol",
                to_kcal(c.value),
                to_kcal(c.error)
            ),
        );
    }
    let half = records[0].energies.len() / 2;
    let first: Vec<Recorded> = records.iter().map(|r| r.slice(0..half)).collect();
    let second: Vec<Recorded> = records
        .iter()
        .map(|r| r.slice(half..r.energies.len()))
        .collect();
    let (h1, h2) = (chain(&first, kt), chain(&second, kt));
    log(
        dir,
        &format!(
            "decoupling: BAR {bar:+.3} ± {bar_error:.3}, TI {ti:+.3} ± {ti_error:.3} kcal/mol; \
             first half {:+.3} ± {:.3}, second half {:+.3} ± {:.3}; EXP forward {:+.3}, reverse \
             {:+.3}; lowest overlap {lowest:.3}",
            to_kcal(h1.value),
            to_kcal(h1.error),
            to_kcal(h2.value),
            to_kcal(h2.error),
            to_kcal(fwd),
            to_kcal(rev)
        ),
    );
    // What is in it, and what is not.
    let tail = kcal(m.decoupling.cross_dispersion_correction());
    let q = &m.benzene.charges;
    let at = &start[..12];
    let mu = [0, 1, 2].map(|k| (0..12).map(|i| q[i] * at[i][k]).sum::<f64>());
    let net: f64 = q.iter().sum();
    let v = m.solvated.cell().volume();
    let dipole_term = kcal(
        2.0 * std::f64::consts::PI * COULOMB * (mu[0] * mu[0] + mu[1] * mu[1] + mu[2] * mu[2])
            / (3.0 * v),
    );
    let alone = PeriodicForceField::new(
        &m.benzene.component,
        &m.benzene.types,
        m.solvated.cell(),
        m.solvated.field.ewald().parameters(),
    )
    .unwrap()
    .with_charges(q.clone())
    .energy(at);
    let images = kcal(alone.total - m.benzene.vacuum.energy(at).total);
    log(
        dir,
        &format!(
            "in it: the long-range van der Waals correction between benzene and the water, \
             {tail:+.4} kcal/mol, linear in λ_v and so exactly its own share of ΔG. Not in it: \
             benzene's net charge {net:+.3e} e, so no Wigner term; its dipole at the start \
             {:.3e} e Å, whose tinfoil term 2πμ²/3V is {dipole_term:.2e} kcal/mol; and its own \
             images, which the decoupling leaves out at every λ, {images:+.4} kcal/mol at the start \
             in this box (the reciprocal sum's bias and the correction's self part included)",
            len(mu) / ANGSTROM
        ),
    );
    log(
        dir,
        &format!(
            "ΔG_hydration = −ΔG_decouple = {:+.3} ± {bar_error:.3} kcal/mol (TI {:+.3} ± \
             {ti_error:.3}); FreeSolv mobley_3053621 −0.90 ± 0.20; OBC II −2.34 polar, −1.13 with \
             the nonpolar term (3c-3)",
            -bar, -ti
        ),
    );
    std::fs::write(
        dir.join("summary.txt"),
        format!("{:e} {:e} {:e} {:e}\n", -bar, bar_error, -ti, ti_error),
    )
    .expect("summary");
}

/// **The cost**, before anything long is run: one step of constrained BAOAB on the decoupled
/// Hamiltonian, the force field's evaluation alone, SETTLE and SHAKE's share, and the couplings at
/// every candidate; the projection of the measurement's wall time from them. Release, one core:
/// `cargo test --release -p pantometry-forcefield --test benzene_hydrated_in_tip3p
/// the_cost_measured -- --ignored --nocapture`.
#[test]
#[ignore = "about a minute with --release: run with --release -- --ignored --nocapture"]
fn the_cost_measured() {
    let m = measurement();
    println!(
        "benzene in {} waters ({} removed from 512), box {:.3} Å, {} wave vectors, α {:.4} Å⁻¹",
        m.solvated.waters,
        512 - m.solvated.waters,
        m.solvated.cell().lengths()[0] / ANGSTROM,
        m.decoupling.field().ewald().wave_vectors(),
        m.decoupling.field().ewald().parameters().alpha * ANGSTROM
    );
    let lambda = Lambda::new(0.0, 0.0, 0.5);
    let potential = AtLambda {
        hamiltonian: &m.decoupling,
        lambda,
    };
    let mut at = m.solvated.at.clone();
    let mut md = m
        .template
        .clone()
        .with_bath(Bath::Langevin {
            temperature: KELVIN,
            friction: 50e12,
            seed: 1,
        })
        .thermalised(&at, KELVIN, 1);
    md.run(&potential, &mut at, 0.5 * FS, 200);
    let mut md = md.with_bath(Bath::Langevin {
        temperature: KELVIN,
        friction: 1e12,
        seed: 2,
    });
    let steps = 100;
    let t = std::time::Instant::now();
    md.run(&potential, &mut at, 2.0 * FS, steps);
    let step = t.elapsed().as_secs_f64() / steps as f64;
    let t = std::time::Instant::now();
    let mut f = vec![[0.0; 3]; at.len()];
    for _ in 0..20 {
        m.decoupling.energy_and_forces(&at, lambda, &mut f);
    }
    let evaluation = t.elapsed().as_secs_f64() / 20.0;
    let settle = Settle::tip3p(m.solvated.field.rigid_waters().to_vec());
    let shake = m.benzene.shake();
    let mut v = md.velocities().to_vec();
    let t = std::time::Instant::now();
    for _ in 0..200 {
        settle.constrain_velocities(&at, &mut v);
        shake.constrain_velocities(md.masses(), md.frozen(), &at, &mut v);
    }
    let projection = t.elapsed().as_secs_f64() / 200.0;
    let states = grid();
    let t = std::time::Instant::now();
    for _ in 0..20 {
        m.decoupling.couplings(&at, &states);
    }
    let couplings = t.elapsed().as_secs_f64() / 20.0;
    let p = PROTOCOL;
    let per_window = p.steps_per_window() as f64 * step + p.samples as f64 * couplings;
    println!(
        "a step {:.2} ms, an evaluation {:.2} ms, the projections {:.3} ms; couplings at all {} \
         candidates {:.2} ms",
        step * 1e3,
        evaluation * 1e3,
        projection * 1e3,
        states.len(),
        couplings * 1e3
    );
    println!(
        "a window of {} steps and {} samples: {:.2} h; {} windows {:.1} h, {} more {:.1} h",
        p.steps_per_window(),
        p.samples,
        per_window / 3600.0,
        SCHEDULE.len(),
        SCHEDULE.len() as f64 * per_window / 3600.0,
        4,
        (SCHEDULE.len() + 4) as f64 * per_window / 3600.0
    );
}

/// **Constrained NVE in the measurement's box**: from the equilibrated start (written by the
/// measurement, or made here), 1 ps of NVE at 2 fs, 1 fs and 0.5 fs from the same velocities, fully
/// coupled. Reported, not asserted: the box's plain cutoff makes the energy step whenever a pair
/// crosses it, which no step size removes, so the `h²` the cluster shows
/// (`tests/constrained_bonds_against_closed_forms.rs`) is here only until those steps dominate.
#[test]
#[ignore = "a few minutes with --release: run with --release -- --ignored --nocapture"]
fn the_box_in_nve_measured() {
    let m = measurement();
    let dir = results_dir();
    let start = equilibrated(&m, &dir);
    let potential = AtLambda {
        hamiltonian: &m.decoupling,
        lambda: Lambda::COUPLED,
    };
    let md0 = m.template.clone().thermalised(&start, KELVIN, 0x4E5);
    for (dt, steps) in [(2.0 * FS, 500usize), (1.0 * FS, 1000), (0.5 * FS, 2000)] {
        let mut md = md0.clone();
        let mut at = start.clone();
        md.prepare(&potential, &at);
        let e0 = md.kinetic_energy() + md.potential_energy().unwrap();
        let (mut squares, mut worst, mut t_sum) = (0.0, 0.0f64, 0.0);
        let mut series = Vec::new();
        for _ in 0..steps {
            md.step(&potential, &mut at, dt);
            let e = md.kinetic_energy() + md.potential_energy().unwrap() - e0;
            squares += e * e;
            worst = worst.max(e.abs());
            t_sum += md.temperature();
            series.push(e);
        }
        let k = series.len();
        let tail: f64 = series[k - k / 10..].iter().sum::<f64>() / (k / 10) as f64;
        let head: f64 = series[..k / 10].iter().sum::<f64>() / (k / 10) as f64;
        println!(
            "dt {:.1} fs: RMS ΔE {:.4e} kcal/mol, largest {:.4e}, last tenth less first {:+.4e} \
             over 0.9 ps, ⟨T⟩ {:.1} K, {} degrees of freedom",
            dt / FS,
            kcal((squares / k as f64).sqrt()),
            kcal(worst),
            kcal(tail - head),
            t_sum / k as f64,
            md.degrees_of_freedom()
        );
    }
}
