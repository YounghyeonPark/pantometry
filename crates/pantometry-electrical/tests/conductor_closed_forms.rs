//! `Conductor` against the resistances a shape actually has.
//!
//! The whole claim of a field formulation is that nobody states a resistance. So every check here
//! computes one from geometry and material, and compares it against what came out of the solve —
//! never against a second solve, and never against a number the domain also produced.

use pantometry_core::units::{Length, Resistivity, Time, Voltage};
use pantometry_core::{Domain, Exchange};
use pantometry_electrical::Conductor;

const COPPER: f64 = 1.724e-8;

fn copper() -> Resistivity {
    Resistivity::ohm_m(COPPER)
}

/// A block `nx` by `ny` by `nz` cells of 1 mm, driven at 1 V.
fn block(counts: (usize, usize, usize)) -> Conductor {
    Conductor::new("bar", counts, Length::mm(1.0), copper(), Voltage::v(1.0))
}

/// **A uniform block's resistance is `ρL/A`, exactly.**
///
/// Not approximately. A cell-centred finite volume with the electrode half a spacing from the
/// first cell centre gives a total series length of `dx/2 + (n−1)dx + dx/2 = n·dx = L`, and the
/// potential in a uniform block is linear, which the discrete operator reproduces with no
/// truncation error at all. So this is machine precision, and a tolerance would be hiding
/// something.
///
/// Three aspect ratios, because `ρL/A` has `L` and `A` in different places and a formulation that
/// swapped them would agree for a cube.
#[test]
fn a_uniform_block_gives_rho_l_over_a_exactly() {
    for counts in [(8, 3, 2), (3, 8, 2), (2, 2, 2), (12, 1, 1)] {
        let mut c = block(counts);
        assert!(c.solve(1e-14), "residual {:.3e}", c.residual());

        let dx = 1e-3;
        let length = counts.0 as f64 * dx;
        let area = (counts.1 * counts.2) as f64 * dx * dx;
        let want = COPPER * length / area;
        let got = c.resistance().to_si();
        assert!(
            (got / want - 1.0).abs() < 1e-12,
            "{counts:?}: {got:.9e} ohm against rho*L/A = {want:.9e}"
        );

        // And the current is V/R, measured at the electrode rather than divided out.
        assert!(
            (c.current().to_si() * want - 1.0).abs() < 1e-12,
            "{counts:?}: I*R should be the 1 V drive, got {:.12}",
            c.current().to_si() * want
        );
        // Charge does not accumulate: what goes in comes out.
        assert!(
            c.current_balance() < 1e-10,
            "{counts:?}: the two electrodes disagree by {:.3e}",
            c.current_balance()
        );
    }
}

/// **Two materials in series add their resistances, and the mean at the interface is harmonic.**
///
/// The classic mistake in a finite-volume conductivity is an arithmetic face mean, which is
/// invisible for a uniform block and wrong the moment two materials meet. So the two here are
/// **four orders of magnitude apart** — copper against something like graphite — where an
/// arithmetic mean would be wrong by nearly a factor of two at the interface.
///
/// The closed form is `ρ₁L₁/A + ρ₂L₂/A`, and it is exact for the discrete problem too, because
/// the potential is piecewise linear and the operator reproduces that.
#[test]
fn two_materials_in_series_add_their_resistances() {
    let (nx, ny, nz) = (10, 3, 3);
    let mut c = block((nx, ny, nz));
    let other = Resistivity::ohm_m(COPPER * 1e4);
    c.set_region(|i, _, _| i >= 4, other);
    assert!(c.solve(1e-14), "residual {:.3e}", c.residual());

    let dx = 1e-3;
    let area = (ny * nz) as f64 * dx * dx;
    let want = COPPER * 4.0 * dx / area + COPPER * 1e4 * 6.0 * dx / area;
    let got = c.resistance().to_si();
    assert!(
        (got / want - 1.0).abs() < 1e-10,
        "series: {got:.9e} against {want:.9e}"
    );

    // An arithmetic face mean would give a different number. Stated so the test is known to
    // discriminate: (sa+sb)/2 at the interface is half the copper conductivity, where the
    // harmonic mean is twice the resistive material's — a factor of 5000 on that one face.
    let arithmetic_face = 0.5 * (1.0 / COPPER + 1.0 / (COPPER * 1e4));
    let harmonic_face = 2.0 / (COPPER + COPPER * 1e4);
    assert!(
        arithmetic_face / harmonic_face > 1000.0,
        "the two means must differ enough for this test to mean anything"
    );
}

/// **Two materials side by side add their conductances.**
///
/// The other half of the same statement, and the one that catches a formulation that treated the
/// transverse direction differently from the axial one. Split across `y`, so the current runs
/// along `x` through both halves independently and the total is `G₁ + G₂`.
#[test]
fn two_materials_in_parallel_add_their_conductances() {
    let (nx, ny, nz) = (8, 4, 3);
    let mut c = block((nx, ny, nz));
    let other = Resistivity::ohm_m(COPPER * 25.0);
    c.set_region(|_, j, _| j >= 2, other);
    assert!(c.solve(1e-14), "residual {:.3e}", c.residual());

    let dx = 1e-3;
    let length = nx as f64 * dx;
    let half_area = (2 * nz) as f64 * dx * dx;
    let g1 = half_area / (COPPER * length);
    let g2 = half_area / (COPPER * 25.0 * length);
    let want = 1.0 / (g1 + g2);
    let got = c.resistance().to_si();
    assert!(
        (got / want - 1.0).abs() < 1e-10,
        "parallel: {got:.9e} against {want:.9e}"
    );

    // The two halves must actually carry different currents, or "parallel" is a word for one
    // material. The resistive half carries 1/25 of the conductive one's density.
    let fast = c.current_density_magnitude(nx / 2, 0, 0).to_si();
    let slow = c.current_density_magnitude(nx / 2, ny - 1, 0).to_si();
    assert!(
        (fast / slow / 25.0 - 1.0).abs() < 0.02,
        "the density ratio should be the conductivity ratio: {:.3}",
        fast / slow
    );
}

/// **`∫σ|∇φ|²dV` equals `V·I`, to machine precision.**
///
/// Tellegen's theorem, and the sharpest single statement about whether the discretisation is
/// self-consistent: the power computed from the *field* must equal the power computed at the
/// *terminals*. They are different sums over different things, and a face conductance that was
/// wrong anywhere would break the identity even where it did not visibly change the resistance.
///
/// Checked on an inhomogeneous block, because on a uniform one a great many wrong formulations
/// still balance.
#[test]
fn the_field_power_equals_the_terminal_power() {
    let mut c = block((7, 4, 3));
    c.set_region(
        |i, j, _| i >= 3 && j >= 2,
        Resistivity::ohm_m(COPPER * 500.0),
    );
    c.set_region(
        |i, j, k| i == 1 && j == 1 && k == 1,
        Resistivity::ohm_m(COPPER * 1e6),
    );
    assert!(c.solve(1e-14), "residual {:.3e}", c.residual());

    let terminal = c.drive().to_si() * c.current().to_si();
    let field = c.dissipation().to_si();
    assert!(
        (field / terminal - 1.0).abs() < 1e-10,
        "Tellegen: field {field:.9e} W against terminals {terminal:.9e} W"
    );
    assert!(terminal > 0.0, "a driven resistor dissipates");

    // And `V²/R` agrees, which ties the reported resistance to the reported power.
    let from_r = c.drive().to_si().powi(2) / c.resistance().to_si();
    assert!((from_r / terminal - 1.0).abs() < 1e-10);
}

/// **Current crowds where the geometry pinches, and the resistance rises above `ρL/A`.**
///
/// The thing a lumped resistor cannot say. A block with an insulating obstruction forces the
/// current through what is left, so its resistance must **exceed** the `ρL/A` of the full section
/// and exceed even the `ρL/A` of the narrowed section — because the current also has to spread
/// back out, and spreading costs.
///
/// A bound rather than a value, because a constriction of this shape has no closed form; the
/// closed forms that exist (`ρ/4a` for a circular contact into a half-space) are limits of
/// geometries this is not. Stating a bound that is *provable* beats quoting a formula that does
/// not apply.
#[test]
fn a_constriction_costs_more_than_its_own_cross_section() {
    let (nx, ny, nz) = (9, 5, 5);
    let mut open = block((nx, ny, nz));
    assert!(open.solve(1e-14));

    let mut pinched = block((nx, ny, nz));
    // A wall across the middle with a one-cell hole in it.
    pinched.set_region(
        |i, j, k| i == nx / 2 && !(j == ny / 2 && k == nz / 2),
        Resistivity::ohm_m(COPPER * 1e12),
    );
    assert!(pinched.solve(1e-12), "residual {:.3e}", pinched.residual());

    let dx = 1e-3;
    let full = open.resistance().to_si();
    let narrow = pinched.resistance().to_si();
    assert!(
        narrow > full * 2.0,
        "a one-cell hole in a 5x5 section should cost a lot: {narrow:.4e} against {full:.4e}"
    );

    // Above the series estimate that ignores spreading: the same block with one cell's worth of
    // 1/25 section and the rest full. Spreading resistance is what the excess *is*.
    let area = (ny * nz) as f64 * dx * dx;
    let naive = COPPER * ((nx - 1) as f64 * dx) / area + COPPER * dx / (dx * dx);
    assert!(
        narrow > naive,
        "spreading should cost more than a plain series estimate: {narrow:.4e} against {naive:.4e}"
    );

    // And the current really does crowd: the density in the hole is far above the mean.
    let in_hole = pinched
        .current_density_magnitude(nx / 2, ny / 2, nz / 2)
        .to_si();
    let far = pinched.current_density_magnitude(0, 0, 0).to_si();
    assert!(
        in_hole > 5.0 * far,
        "the hole should carry the crowding: {in_hole:.3e} against {far:.3e} A/m2"
    );
}

/// **A solve that does not converge is refused, not returned.**
///
/// The failure this domain is built to avoid. An iterative solver stopped at its cap produces a
/// field shaped exactly like an answer — smooth, bounded, roughly right in the middle — and
/// nothing downstream can tell. So `step` returns a `Violation` naming the residual.
///
/// Provoked by starving the solver of iterations rather than by asking for an unreachable
/// tolerance. The first draft asked for `1e-300` and **passed anyway**: on a uniform block the
/// residual reaches exactly zero, so the tolerance was met and there was nothing to refuse. A
/// failure mode that cannot be provoked is a failure mode that is not tested, which is why
/// `with_solver` exists.
#[test]
fn a_solve_that_did_not_converge_is_refused() {
    // Starved **and** given a different problem from the one the constructor solved. Setting the
    // budget alone is not enough now that `new` solves: the stored field is already the answer,
    // so a one-iteration solve meets any tolerance immediately. Changing the material is what
    // makes the stored field wrong, and the budget is what stops it being fixed.
    let mut starved = block((6, 4, 4)).with_solver(1e-14, 1);
    starved.set_region(|i, _, _| i >= 3, Resistivity::ohm_m(COPPER * 1e8));
    let mut c = starved.clone();
    assert!(
        !c.solve(1e-14),
        "one iteration cannot re-solve 96 changed cells"
    );
    assert!(!c.converged());
    assert!(
        c.residual() > 1e-14 && c.residual().is_finite(),
        "and it reports what it did reach: {:.3e}",
        c.residual()
    );

    // What it left behind looks like an answer, which is the whole reason to refuse it: bounded,
    // between the electrodes, and wrong.
    let (nx, ny, nz) = c.counts();
    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let v = c.potential_at(i, j, k).to_si();
                assert!((0.0..=1.0).contains(&v), "({i},{j},{k}) is at {v} V");
            }
        }
    }

    // `step` refuses it, naming the residual it reached.
    let mut c = starved;
    let err = c
        .step(Time::from_si(0.0), Time::from_si(1.0), &mut Exchange::new())
        .expect_err("a domain must not publish heat it computed from a half-solve");
    assert_eq!(err.quantity, "solver residual");
    assert!(err.after > 1e-14, "the residual reached: {}", err.after);

    // And the same block with a real budget converges, so the refusal is about the budget.
    let mut ok = block((6, 4, 4));
    assert!(ok.solve(1e-12));
    assert!(ok.converged());
}

/// **Inside a simulation it publishes what it dissipates, and the books close.**
///
/// A quasi-static domain paying joules onto the bus, with its own ledger going negative to match.
/// If the ledger were empty the audit would see energy appear, which is the mistake a source
/// makes exactly once.
#[test]
fn it_pays_its_dissipation_onto_the_bus() {
    use pantometry_core::conserved::quantity;
    use pantometry_core::{Kind, Ledger, Schedule, Simulation, Violation};

    /// A sink that takes whatever is offered and holds it.
    struct Sink {
        held: f64,
    }
    impl Domain for Sink {
        fn name(&self) -> &str {
            "sink"
        }
        fn kind(&self) -> Kind {
            Kind::Evolving
        }
        fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
            self.held += bus.take_share(quantity::ENERGY, dt);
            Ok(())
        }
        fn ledger(&self) -> Ledger {
            Ledger::new().with(quantity::ENERGY, self.held)
        }
    }

    let seconds = 0.25;
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(block((6, 3, 3)))
        .with(Sink { held: 0.0 });
    sim.advance(Time::from_si(seconds))
        .expect("what it pays, the sink takes");

    let c = sim.domain_as::<Conductor>("bar").expect("still there");
    let want = c.dissipation().to_si() * seconds;
    assert!(
        (c.dissipated_energy().to_si() / want - 1.0).abs() < 1e-9,
        "{:.6e} J spent against {want:.6e}",
        c.dissipated_energy().to_si()
    );

    // The resistance is a real copper number: 6 mm of 3x3 mm section is 11.5 microhm, so 1 V
    // across it is a wildly unphysical 87 kA. That is fine — it is a test of the arithmetic, and
    // saying so beats quietly choosing a voltage that looks sensible.
    let r = c.resistance().to_si();
    assert!((r - COPPER * 6e-3 / 9e-6).abs() / r < 1e-12, "{r:.6e} ohm");
}

/// **A path that must detour through y, and one that must detour through z.**
///
/// Written because a mutation walked straight past everything above. Deleting **every
/// z-direction face** from the operator left six of the seven tests green: in all of them the
/// current runs along x and the potential is uniform across the transverse axes, so the faces
/// that were removed carried nothing.
///
/// A labyrinth fixes that. Two layers, each blocked by an insulating wall at a different `x`, so
/// the only conducting path crosses from one layer to the other and back. Without the transverse
/// faces there is no path at all and the resistance is the insulator's — twelve orders of
/// magnitude up.
///
/// Both axes, because y and z are separate arrays and separate loops, and a formulation can
/// perfectly well have one and not the other.
#[test]
fn a_detour_needs_the_transverse_faces() {
    let dx = 1e-3;
    let straight = {
        let mut c = block((9, 2, 2));
        assert!(c.solve(1e-14));
        c.resistance().to_si()
    };

    // `across` picks which transverse axis the layers are stacked along.
    let labyrinth = |across_z: bool| {
        let mut c = block((9, 2, 2));
        let wall = Resistivity::ohm_m(COPPER * 1e12);
        c.set_region(
            |i, j, k| {
                let layer = if across_z { k } else { j };
                (layer == 0 && i == 3) || (layer == 1 && i == 6)
            },
            wall,
        );
        assert!(c.solve(1e-12), "residual {:.3e}", c.residual());
        c.resistance().to_si()
    };

    for (axis, r) in [("z", labyrinth(true)), ("y", labyrinth(false))] {
        // A path exists, so the resistance is a metal's rather than an insulator's.
        assert!(
            r < 1e4 * straight,
            "detouring through {axis}: {r:.4e} ohm is an insulator, not a path — \
             the transverse faces are missing"
        );
        // And it costs something: the current is squeezed into one layer for part of the way.
        assert!(
            r > 1.5 * straight,
            "detouring through {axis} should cost more than going straight: \
             {r:.4e} against {straight:.4e}"
        );
    }

    // The two are the same problem rotated, so they must agree exactly. This is the check a
    // formulation that had y faces and not z ones fails while passing everything else.
    assert!(
        (labyrinth(true) / labyrinth(false) - 1.0).abs() < 1e-9,
        "y and z are the same physics rotated: {:.9e} against {:.9e}",
        labyrinth(true),
        labyrinth(false)
    );

    // And current really crosses between layers, or the labyrinth is a word. The transverse
    // component of J at the crossing is a real fraction of the axial one.
    let mut c = block((9, 2, 2));
    c.set_region(
        |i, _, k| (k == 0 && i == 3) || (k == 1 && i == 6),
        Resistivity::ohm_m(COPPER * 1e12),
    );
    assert!(c.solve(1e-12));
    let j = c.current_density_at(4, 0, 0);
    assert!(
        j.z.abs() > 0.05 * j.x.abs().max(1e-30),
        "the current should be crossing layers at the wall: J = {j:?}"
    );
    let _ = dx;
}

/// **A busbar's resistance is `ρL/A` at every drive, switched off included.**
///
/// The solve was at the drive, so its tolerance was relative to the drive. Turned down a
/// millionfold from a warm start it stopped with the 1 V field's error still in it, `2.5e-9` from
/// `ρL/A` — which every other test here holds to `1e-12`. Switched off, the right-hand side was
/// zero and the residual was taken against `|b|` floored at `f64::MIN_POSITIVE`: the solve
/// "converged" only because its recursive residual underflowed first, and left `1e-16 V` and
/// `−5.3e-14 A` behind it, and the resistance was `0/0`. Now the field at one volt is solved and
/// scaled, and a resistance, which is a property of the shape, has a value at every drive.
///
/// So the drive is walked down through 1 V, 1 µV, 0 and −1 mV, and then **up**, through 1 mV to
/// 1 V. At every one the field is **bit for bit** a fresh busbar's at that drive, read straight
/// after `set_drive` and again after a solve: a switch that left anything of the last drive
/// behind would show there. The way up is the half that needs it — a start a thousandth of the
/// answer is a good warm start, and a solve at the drive takes it: measured, its first cell lands
/// 55 units in the last place from a fresh busbar's. The resistance is `ρL/A` to `1e-12`
/// throughout, and at zero the potential, the current and the dissipation are exactly zero.
#[test]
fn a_busbar_has_its_resistance_at_every_drive_including_none() {
    let counts = (20, 4, 4);
    let dx = 1e-3;
    let want = COPPER * (counts.0 as f64 * dx) / ((counts.1 * counts.2) as f64 * dx * dx);
    let mut c = block(counts);
    let same_as_fresh = |c: &Conductor, volts: f64, when: &str| {
        let fresh = Conductor::new("bar", counts, Length::mm(1.0), copper(), Voltage::v(volts));
        for k in 0..counts.2 {
            for j in 0..counts.1 {
                for i in 0..counts.0 {
                    assert_eq!(
                        c.potential_at(i, j, k).to_si().to_bits(),
                        fresh.potential_at(i, j, k).to_si().to_bits(),
                        "at {volts} V {when}, cell ({i},{j},{k}) differs from a fresh busbar"
                    );
                }
            }
        }
    };
    for volts in [1.0, 1e-6, 0.0, -1e-3, 1e-3, 1.0] {
        c.set_drive(Voltage::v(volts));
        same_as_fresh(&c, volts, "before a solve");
        assert!(
            c.solve(1e-12),
            "at {volts} V: residual {:.3e}",
            c.residual()
        );
        same_as_fresh(&c, volts, "after a solve");
        let got = c.resistance().to_si();
        assert!(
            (got / want - 1.0).abs() < 1e-12,
            "at {volts} V: {got:.9e} ohm against rho*L/A = {want:.9e}"
        );
        assert!(
            c.current_balance() < 1e-10,
            "at {volts} V: the electrodes disagree by {:.3e}",
            c.current_balance()
        );
        // The resistance and the balance are read from the field at one volt and cannot see the
        // drive, so the field *at* the drive is held to its two closed forms as well — Ohm's law
        // at the electrode, to the `1e-12` the uniform block is held to, and Tellegen's
        // `∫σ|∇φ|² = V·I`, to the `1e-10` its own test uses. A field scaled by `|V|`, or a
        // current that lost its sign, passes every check above and neither of these.
        if volts != 0.0 {
            let (v, i) = (volts, c.current().to_si());
            assert!(
                (i * want / v - 1.0).abs() < 1e-12,
                "at {volts} V: I*R/V is {:.15}",
                i * want / v
            );
            let field = c.dissipation().to_si();
            assert!(
                (field / (v * i) - 1.0).abs() < 1e-10,
                "at {volts} V: field {field:.9e} W against terminals {:.9e} W",
                v * i
            );
        }
        if volts == 0.0 {
            for k in 0..counts.2 {
                for j in 0..counts.1 {
                    for i in 0..counts.0 {
                        assert_eq!(c.potential_at(i, j, k).to_si(), 0.0, "cell ({i},{j},{k})");
                    }
                }
            }
            assert_eq!(c.current().to_si(), 0.0);
            assert_eq!(c.dissipation().to_si(), 0.0);
        }
    }
}

/// **A busbar cut through at its driven end is an open circuit, and an open circuit is an answer.**
///
/// A perfect insulator across the driven layer, after a solve that conducted. The right-hand side
/// at one volt is then exactly zero — nothing reaches the electrode — and the warm start is the
/// conducting field, which answers nothing about this one. Every cell with a path to ground has
/// the zero solution, and a warm start is the one start that cannot reach it: the residual is
/// relative to `|b|`, floored at `f64::MIN_POSITIVE`. So the solve must start from zero, and then:
/// converged, no current, a resistance that is infinite rather than `NaN`, and every cell still
/// joined to the grounded electrode at exactly its potential.
#[test]
fn a_busbar_cut_at_its_driven_end_is_an_open_circuit() {
    let counts = (12, 3, 3);
    let mut c = block(counts);
    assert!(c.converged(), "the busbar conducts first");
    c.set_region(
        |i, _, _| i == counts.0 - 1,
        Resistivity::ohm_m(f64::INFINITY),
    );
    assert!(
        c.solve(1e-12),
        "an open circuit solves; residual {:.3e}",
        c.residual()
    );
    assert_eq!(c.current().to_si(), 0.0);
    assert_eq!(c.resistance().to_si(), f64::INFINITY);
    for k in 0..counts.2 {
        for j in 0..counts.1 {
            for i in 0..counts.0 - 1 {
                assert_eq!(c.potential_at(i, j, k).to_si(), 0.0, "cell ({i},{j},{k})");
            }
        }
    }
}

/// **A material that is not a number gives a field that says so, and the copper put back gives
/// the busbar back.**
///
/// One cell given a resistivity that is not a number makes `A·x` `NaN` from any start, so the
/// residual is `NaN`, and `NaN > tolerance` is false: taken as it stood, the solve stopped before
/// its first iteration and left **the last field standing** — refused, and shaped exactly like an
/// answer, which is what the committed solver did. A residual that is not a number now sends the
/// solve back to zero, where it runs, and the field it leaves is `NaN`, which is the truth about
/// it. Then the copper is put back, the warm start is that `NaN` field, and it is no start either:
/// from zero the busbar solves, and is bit for bit a fresh one.
#[test]
fn a_field_that_stopped_being_a_number_does_not_poison_the_next_solve() {
    let counts = (8, 3, 3);
    let mut c = block(counts);
    c.set_resistivity(3, 1, 1, Resistivity::ohm_m(f64::NAN));
    assert!(
        !c.solve(1e-12),
        "a material that is not a number is refused"
    );
    let poisoned = (0..counts.2).any(|k| {
        (0..counts.1).any(|j| (0..counts.0).any(|i| c.potential_at(i, j, k).to_si().is_nan()))
    });
    assert!(
        poisoned,
        "a refused material must leave a field that says so, not the last one standing"
    );

    c.set_resistivity(3, 1, 1, copper());
    assert!(
        c.solve(1e-12),
        "the copper is back and the busbar must solve; residual {:e}",
        c.residual()
    );
    let fresh = block(counts);
    for k in 0..counts.2 {
        for j in 0..counts.1 {
            for i in 0..counts.0 {
                assert_eq!(
                    c.potential_at(i, j, k).to_si().to_bits(),
                    fresh.potential_at(i, j, k).to_si().to_bits(),
                    "cell ({i},{j},{k})"
                );
            }
        }
    }
}
