//! **Every domain that says it can be rewound can be rewound more than once.**
//!
//! `Schedule::Iterative` checkpoints every domain once an advance and restores it before every
//! sweep after the first, so a coupling that takes three sweeps restores one checkpoint twice. A
//! restore that moved the checkpoint out rather than copying it gives back the right state the
//! first time, finds nothing the second, does nothing, and the third sweep begins where the second
//! ended: the domain marches two windows in one and every number it reports looks like a
//! converged coupling. `Puck` and `Waves` both did it. The puck did not claim
//! `Domain::supports_restore`, so nothing could reach it; `Waves` did.
//!
//! A domain's own tests restore once, which cannot see it. So this holds each of the eighteen
//! domains in the library that claim a restore, one of them twice, to the schedule's use of one:
//! checkpointed,
//! stepped, and restored three times, each restore giving back what the checkpoint showed and each
//! run after it the steps the first run took, to the bit; and then checkpointed again, which has
//! to replace the first. Checkpointed before the first step, as an advance from a fresh domain is,
//! and again after three.
//!
//! It found six more that did not give back what they had checkpointed.
//!
//! - `Tube`, `Room` and `Hall` put back neither the startup offset their first step sets nor the
//!   pressure a step before, which they set to the pressure now, and their energy reads both: a
//!   room restored three steps in reported 2.2% less than the same field had at its checkpoint.
//!   The room and the hall audit their own books, so `Schedule::Iterative` refused the second
//!   sweep of the first window it gave either.
//! - `Bar1D`, `Tube` and `Solid3D` kept what they had absorbed or radiated in a sweep that was
//!   thrown away, and `Fluid` the work its thermostat had done, which its ledger subtracts.
//! - A melting `Solid3D` kept the capacities, conductances and limit its last step had derived
//!   from fractions the restore had just put back.
//!
//! Three of those passed the first version of this file, because its entries never moved what
//! they had lost: no heat reached the block, and the fluid was a lattice at 1.2 K where its units
//! put a liquid at 8.7e22 K, so the bath's work was below an ulp of the ledger.
//!
//! # What it can see
//!
//! What a domain shows through the trait: its readings, its ledger, its bodies, its limit and
//! residual, and what each step left on the bus — and, where those do not show the state itself,
//! the domain's own accessors. Every step has to move something there, or the comparison would
//! pass on a domain that never moved; that is asserted, not assumed. What it does not assert is
//! that every *field* moves: an entry is built so that each total its domain keeps is fed — heat
//! on the bus, a slipping ball, a contact that damps, a bath — and a counter none of them reaches
//! is invisible here. A field restored wrongly that changes nothing shown within three steps is
//! for the domain's own tests.
//!
//! # Which domains
//!
//! The list is checked against the source: `every_claim_in_the_library_is_in_the_list` reads
//! every `impl Domain for` in `crates/` and fails on a claim the list has not got, on a listing that
//! claims nothing, and on a header or a body it cannot read, so a new one cannot skip this. The three in `app/pantometry-world` are another workspace's, and the
//! application offers no iterative schedule.

use std::any::type_name;

use glam::DVec3;
use pantometry::elastic::{Axis, Waves};
use pantometry::forcefield::Bath;
use pantometry::mechanics::{Rolling, Sphere, Surface};
use pantometry::prelude::*;
use pantometry::units::VelocityVec;

/// How to build a domain moving, and what of it to look at.
struct Entry<T> {
    build: fn() -> T,
    /// The step, when the domain's own limit is not finite.
    dt: Option<Time>,
    /// What the bus carries into every step, for a domain that moves only when fed.
    feed: fn(&mut Exchange),
    /// What the trait does not show.
    probe: fn(&T) -> Vec<String>,
}

fn claimant<T>(build: fn() -> T) -> Entry<T> {
    Entry {
        build,
        dt: None,
        feed: |_| {},
        probe: |_| Vec::new(),
    }
}

/// An entry with its type rubbed out, named by the type it builds, so an entry cannot name one
/// type and build another.
struct Claimant {
    ty: String,
    check: Box<dyn Fn() -> Result<(), String>>,
}

impl<T: Domain + 'static> Entry<T> {
    fn stepped_by(mut self, dt: Time) -> Entry<T> {
        self.dt = Some(dt);
        self
    }

    fn fed(mut self, feed: fn(&mut Exchange)) -> Entry<T> {
        self.feed = feed;
        self
    }

    fn showing(mut self, probe: fn(&T) -> Vec<String>) -> Entry<T> {
        self.probe = probe;
        self
    }

    fn listed(self) -> Claimant {
        let full = type_name::<T>();
        let krate = full.split("::").next().expect("a path").replace('_', "-");
        let last = full.rsplit("::").next().expect("a path");
        let ty = format!("{krate}::{last}");
        let name = ty.clone();
        Claimant {
            ty,
            check: Box::new(move || holds(&name, &self)),
        }
    }
}

/// A temperature in the units a Lennard-Jones fluid of unit `ε` and `σ` is stated in.
fn reduced(t: f64) -> Temperature {
    pantometry::molecular::temperature_from_reduced(t, &LennardJones::reduced())
}

fn bits(x: f64) -> String {
    format!("{x:e} ({:016x})", x.to_bits())
}

/// Every domain in the library that claims `Domain::supports_restore`, each built moving.
fn claimants() -> Vec<Claimant> {
    vec![
        // A matched end, so it radiates from the first step and the count of it moves.
        claimant(|| {
            let tube = Tube::of_air("tube", Length::m(1.0), 40, Area::from_si(1e-4));
            let matched = End::Matched(tube.impedance());
            tube.with_ends(End::Closed, matched).released_from(|x| {
                let s = (x.to_si() - 0.95) / 0.05;
                Pressure::from_si((-(s * s)).exp())
            })
        })
        .showing(|t: &Tube| {
            let mut out: Vec<String> = (0..t.cells())
                .map(|i| format!("p {i} {}", bits(t.pressure_at(i).to_si())))
                .collect();
            out.extend(
                (0..t.cells() - 1).map(|i| format!("u {i} {}", bits(t.velocity_at(i).to_si()))),
            );
            out.push(format!("radiated {}", bits(t.radiated_energy().to_si())));
            out
        })
        .listed(),
        claimant(|| {
            Room::of_air("room", Length::m(4.4), Length::m(3.1), 21).released_in_mode(
                1,
                1,
                Pressure::from_si(1.0),
            )
        })
        .listed(),
        claimant(|| {
            Hall::of_air("hall", Length::m(4.4), Length::m(3.1), Length::m(2.4), 9)
                .released_in_mode((1, 1, 1), Pressure::from_si(1.0))
        })
        .listed(),
        claimant(|| {
            let mut w = Waves::new(
                "waves",
                (1, 1, 8),
                Length::mm(1.0),
                Elastic::aluminium_6061(),
            );
            w.hold(Axis::X);
            w.hold(Axis::Y);
            w.clamp_ends(Axis::Z);
            w.release_mode(1, Axis::Z, Axis::Z, Length::from_si(1e-9));
            w
        })
        .showing(|w: &Waves| {
            let (nx, ny, nz) = w.node_counts();
            let mut out = Vec::new();
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        let u = w.displacement_at(i, j, k);
                        out.push(format!(
                            "u {i} {j} {k} {} {} {}",
                            bits(u[0]),
                            bits(u[1]),
                            bits(u[2])
                        ));
                    }
                }
            }
            out
        })
        .listed(),
        claimant(|| {
            Molecule::new(
                "aspirin",
                Component::from_ccd(include_str!(
                    "../../pantometry-forcefield/components/AIN.cif"
                ))
                .expect("AIN parses"),
            )
            .thermalised(300.0, Bath::Isolated, 1)
        })
        .listed(),
        claimant(|| {
            Rolling::new(
                "slope",
                Sphere::solid(Mass::kg(0.2), Length::m(0.05), LengthVec::ZERO),
                0.35,
                // Under the 0.104 a solid ball needs to roll on this slope, so it slips, and the
                // friction dissipates from the first step.
                Surface::new(0.0, 0.05),
            )
        })
        .stepped_by(Time::ms(1.0))
        .showing(|r: &Rolling| {
            vec![
                format!("speed {}", bits(r.speed())),
                format!("spin {}", bits(r.spin())),
                format!("distance {}", bits(r.distance().to_si())),
                format!("dissipated {}", bits(r.dissipated_energy().to_si())),
            ]
        })
        .listed(),
        claimant(|| {
            NBody::new(
                "pair",
                &[
                    Body::new(Mass::kg(6e24), LengthVec::ZERO, VelocityVec::ZERO),
                    Body::new(
                        Mass::kg(1e3),
                        LengthVec::m(7e6, 0.0, 0.0),
                        VelocityVec::m_per_s(0.0, 7.5e3, 0.0),
                    ),
                ],
            )
        })
        .stepped_by(Time::s(10.0))
        .listed(),
        claimant(|| {
            TreeNBody::new(
                "three",
                &[
                    Body::new(Mass::kg(6e24), LengthVec::ZERO, VelocityVec::ZERO),
                    Body::new(
                        Mass::kg(1e3),
                        LengthVec::m(7e6, 0.0, 0.0),
                        VelocityVec::m_per_s(0.0, 7.5e3, 0.0),
                    ),
                    Body::new(
                        Mass::kg(1e3),
                        LengthVec::m(0.0, -9e6, 0.0),
                        VelocityVec::m_per_s(6.6e3, 0.0, 0.0),
                    ),
                ],
            )
        })
        .stepped_by(Time::s(10.0))
        .listed(),
        claimant(|| {
            ContactSystem::new(
                "ball",
                // A tenth of a millimetre into the floor and still going in, so the contact
                // damps from the first step.
                &[Body::new(
                    Mass::kg(0.1),
                    LengthVec::m(0.0, 0.0, -1e-4),
                    VelocityVec::m_per_s(0.0, 0.0, -0.5),
                )],
                AccelerationVec::from_si(DVec3::new(0.0, 0.0, -9.81)),
                Ground::floor(),
                Stiffness::from_si(1e5),
                Damping::from_si(1.0),
            )
        })
        .showing(|c: &ContactSystem| {
            vec![format!(
                "dissipated {}",
                bits(c.dissipated_energy().to_si())
            )]
        })
        .listed(),
        claimant(|| {
            RigidBody::new(
                "brick",
                Inertia::solid_box(Mass::kg(1.0), DVec3::new(0.1, 0.2, 0.3)),
            )
            .spinning_body_frame(DVec3::new(1.0, 0.2, 0.1))
        })
        .stepped_by(Time::ms(1.0))
        .showing(|b: &RigidBody| {
            let q = b.orientation();
            let w = b.angular_velocity_body();
            vec![
                format!(
                    "orientation {} {} {} {}",
                    bits(q.x),
                    bits(q.y),
                    bits(q.z),
                    bits(q.w)
                ),
                format!("omega {} {} {}", bits(w.x), bits(w.y), bits(w.z)),
            ]
        })
        .listed(),
        claimant(|| {
            Fluid::lattice(
                "fluid",
                LennardJones::reduced(),
                Mass::from_si(1.0),
                3,
                0.8442,
            )
            .thermalised(reduced(1.2), 7)
            // A bath, so the work it does moves.
            .with_thermostat(Thermostat::Langevin {
                target: reduced(1.2),
                damping: 1.0,
            })
        })
        .stepped_by(Time::from_si(0.002))
        // Shown alone as well as in the ledger, where it is subtracted from an energy near -700
        // and a lattice too cold to move would have it below the ulp.
        .showing(|f: &Fluid| vec![format!("bath {}", bits(f.thermostat_work().to_si()))])
        .listed(),
        claimant(|| {
            let mut m = CompartmentModel::new("pk");
            let plasma = m.eliminating("plasma", Volume::litres(5.0), Clearance::l_per_h(5.0));
            m.bolus(plasma, Mass::kg(1e-4)).expect("a dose");
            m
        })
        .stepped_by(Time::s(60.0))
        .listed(),
        claimant(|| {
            let dx = 1e-3;
            let mut p = Puck::new(
                "puck",
                Basket {
                    counts: (8, 6, 8),
                    cell: Length::from_si(dx),
                    radius: Length::from_si(4.0 * dx),
                    porosity: 0.45,
                    grind: Grind::sieved(Length::from_si(800e-6)),
                    ..Basket::espresso()
                },
            );
            p.set_temperature(Temperature::celsius(20.0));
            p.set_inlet_temperature(Temperature::celsius(93.0));
            p
        })
        .listed(),
        claimant(|| {
            Well::new("well", Mass::kg(9.109e-31), Length::nm(2.0), 101).with_gaussian(
                Length::nm(0.8),
                Length::nm(0.15),
                Wavenumber::from_si(4.0e9),
            )
        })
        .listed(),
        claimant(|| {
            LumpedMass::new(
                "plate",
                Substance::aluminium_6061(),
                Volume::cm3(10.0),
                Length::mm(2.0),
                Temperature::celsius(80.0),
                Environment::still_air(Temperature::celsius(20.0), Area::cm2(50.0)),
            )
        })
        .stepped_by(Time::s(1.0))
        .fed(|bus| bus.publish(HEAT, 0.01))
        .listed(),
        claimant(|| {
            Bar1D::new(
                "bar",
                Substance::aluminium_6061(),
                11,
                Length::mm(1.0),
                Area::mm2(10.0),
                Temperature::celsius(20.0),
            )
        })
        .fed(|bus| bus.publish(HEAT, 0.01))
        .listed(),
        claimant(|| {
            let ambient = Temperature::celsius(25.0);
            let mut net = ThermalNetwork::new("motor");
            let w = net.node(
                "winding",
                Substance::copper(),
                Volume::cm3(18.0),
                Length::mm(2.0),
                Temperature::celsius(80.0),
            );
            let h = net.node_losing_to(
                "housing",
                Substance::aluminium_6061(),
                Volume::cm3(220.0),
                Length::mm(4.0),
                ambient,
                Environment::still_air(ambient, Area::cm2(420.0)),
            );
            net.link(w, h, Conductance::w_per_k(0.9)).expect("a link");
            net.absorbing(w).expect("an absorbing node");
            net
        })
        .stepped_by(Time::s(1.0))
        .fed(|bus| bus.publish(HEAT, 0.01))
        .showing(|n: &ThermalNetwork| {
            vec![
                format!("absorbed {}", bits(n.absorbed_energy().to_si())),
                format!("lost {}", bits(n.lost_energy().to_si())),
            ]
        })
        .listed(),
        claimant(|| {
            Solid3D::new(
                "block",
                Substance::copper(),
                (4, 3, 2),
                Length::mm(1.0),
                Temperature::celsius(80.0),
            )
            .losing_from(
                pantometry::thermal::Face::XMin,
                Environment::still_air(Temperature::celsius(20.0), Area::cm2(1.0)),
            )
        })
        .fed(|bus| bus.publish(HEAT, 0.01))
        .showing(|s: &Solid3D| {
            let (nx, ny, nz) = s.counts();
            let mut out = Vec::new();
            for k in 0..nz {
                for j in 0..ny {
                    for i in 0..nx {
                        out.push(format!(
                            "T {i} {j} {k} {}",
                            bits(s.temperature_at(i, j, k).to_si())
                        ));
                    }
                }
            }
            out
        })
        .listed(),
        // Twice, because a block with a liquid phase is another path: its step re-derives the
        // capacities, the conductances and the limit from the fractions it leaves behind.
        claimant(|| {
            let ice = Substance::ice();
            let fusion = ice.fusion.expect("ice melts");
            let two_phase = ice.with_fusion(fusion.with_liquid(ThermalProps {
                conductivity: ThermalConductivity::w_per_m_k(0.598),
                specific_heat: SpecificHeat::j_per_kg_k(4182.0),
                expansion: ThermalExpansion::ppm_per_k(207.0),
                emissivity: 0.96,
            }));
            // Ice at its melting point and fed, so every step melts some of every cell. Not set
            // with `set_melted_fraction`, which leaves the operator the constructor derived until
            // the first step: a checkpoint then would hold a limit that no restore should give.
            Solid3D::new(
                "slush",
                two_phase,
                (4, 1, 1),
                Length::mm(1.0),
                Temperature::celsius(0.0),
            )
            .losing_from(
                pantometry::thermal::Face::XMin,
                Environment::still_air(Temperature::celsius(-20.0), Area::cm2(1.0)),
            )
        })
        .fed(|bus| bus.publish(HEAT, 0.01))
        .showing(|s: &Solid3D| {
            (0..4)
                .map(|i| {
                    format!(
                        "T {i} {}, melted {}",
                        bits(s.temperature_at(i, 0, 0).to_si()),
                        bits(s.melted_fraction_at(i, 0, 0))
                    )
                })
                .collect()
        })
        .listed(),
    ]
}

/// What the test can see of a domain: through the trait, through its probe, and on the bus.
fn seen<T: Domain>(e: &Entry<T>, d: &T, now: Time, bus: &Exchange) -> Vec<String> {
    let mut out = Vec::new();
    for r in d.readings() {
        out.push(format!(
            "reading {} {} {} {}",
            r.domain,
            r.label,
            bits(r.value),
            r.unit
        ));
    }
    for (q, v) in d.ledger().quantities() {
        out.push(format!("ledger {q} {}", bits(v)));
    }
    if let Some(b) = d.as_bodies() {
        for i in 0..b.count() {
            let p = b.position(i).to_si();
            out.push(format!(
                "body {i} at {} {} {}, {}",
                bits(p.x),
                bits(p.y),
                bits(p.z),
                bits(b.value(i))
            ));
        }
    }
    out.push(format!("residual {}", bits(d.residual())));
    out.push(format!("limit {}", bits(d.max_stable_dt(now).to_si())));
    out.extend((e.probe)(d));
    out.push(format!("bus {bus:?}"));
    out
}

/// The first line two views differ in, for a message that says where rather than dumping both.
fn first_difference(a: &[String], b: &[String]) -> String {
    match a.iter().zip(b).position(|(x, y)| x != y) {
        Some(i) => format!("first at line {i}: {:?} against {:?}", a[i], b[i]),
        None => format!("{} lines against {}", a.len(), b.len()),
    }
}

const STEPS: usize = 3;

/// `STEPS` steps from `from`, with what was seen after each.
fn run<T: Domain>(
    ty: &str,
    e: &Entry<T>,
    d: &mut T,
    from: Time,
    dt: Time,
) -> Result<Vec<Vec<String>>, String> {
    let mut views = Vec::new();
    for n in 0..STEPS {
        let t = from + dt * n as f64;
        let mut bus = Exchange::new();
        (e.feed)(&mut bus);
        d.step(t, dt, &mut bus)
            .map_err(|v| format!("{ty}: step {n} was refused: {v:?}"))?;
        views.push(seen(e, d, t + dt, &bus));
    }
    Ok(views)
}

/// Both ways a coupling checkpoints: after some steps, and before the first, which is where an
/// advance from a fresh domain starts and where a step's one-off startup has not happened yet.
///
/// In that order, so each has a defect only it reports. The acoustic restores that set the
/// pressure a step before to the pressure now cost nothing before the first step, where the two
/// are equal; their missing startup offset costs nothing after it, where the offset has stopped
/// moving.
fn holds<T: Domain>(ty: &str, e: &Entry<T>) -> Result<(), String> {
    holds_after(ty, e, STEPS).map_err(|m| format!("{m} (checkpointed after {STEPS} steps)"))?;
    holds_after(ty, e, 0).map_err(|m| format!("{m} (checkpointed before the first step)"))
}

/// Stepped `warm` times, checkpointed, stepped, restored three times, then checkpointed again;
/// `Err` says what failed.
fn holds_after<T: Domain>(ty: &str, e: &Entry<T>, warm: usize) -> Result<(), String> {
    let mut d = (e.build)();
    if !d.supports_restore() {
        return Err(format!("{ty} is listed here and does not claim a restore"));
    }
    let start = Time::ZERO;
    let dt = match e.dt {
        Some(dt) => dt,
        None => {
            let limit = d.max_stable_dt(start);
            if !(limit.to_si().is_finite() && limit > Time::ZERO) {
                return Err(format!(
                    "{ty} has no finite limit, so its entry has to state a step"
                ));
            }
            limit * 0.5
        }
    };
    let mut views = Vec::new();
    for _ in 0..warm / STEPS {
        let from = start + dt * views.len() as f64;
        views.extend(run(ty, e, &mut d, from, dt)?);
    }
    let at = start + dt * views.len() as f64;
    d.checkpoint();
    let mark = seen(e, &d, at, &Exchange::new());

    let first = run(ty, e, &mut d, at, dt)?;
    // Each step against the one before it, the first against the checkpoint, leaving out the bus,
    // which a step that moved nothing could still have written to.
    let steps: Vec<&Vec<String>> = std::iter::once(&mark).chain(&first).collect();
    for n in 1..steps.len() {
        let (a, b) = (steps[n - 1], steps[n]);
        if a[..a.len() - 1] == b[..b.len() - 1] {
            return Err(format!(
                "{ty}: step {} moved nothing this test can see, so a restore could not be checked",
                n - 1
            ));
        }
    }
    for again in 1..=3 {
        d.restore();
        let back = seen(e, &d, at, &Exchange::new());
        if back != mark {
            return Err(format!(
                "{ty}: restore {again} did not give back the checkpointed state, {}",
                first_difference(&back, &mark)
            ));
        }
        let rerun = run(ty, e, &mut d, at, dt)?;
        if let Some(n) = (0..STEPS).find(|&n| rerun[n] != first[n]) {
            return Err(format!(
                "{ty}: step {n} after restore {again} is not step {n} after the checkpoint, {}",
                first_difference(&rerun[n], &first[n])
            ));
        }
    }
    // A second checkpoint replaces the first.
    let later = at + dt * STEPS as f64;
    d.checkpoint();
    let mark2 = seen(e, &d, later, &Exchange::new());
    if mark2 == mark {
        return Err(format!(
            "{ty}: three steps on looks like the first checkpoint"
        ));
    }
    run(ty, e, &mut d, later, dt)?;
    d.restore();
    let back = seen(e, &d, later, &Exchange::new());
    if back != mark2 {
        return Err(format!(
            "{ty}: a restore after a second checkpoint did not go back to the second, {}",
            first_difference(&back, &mark2)
        ));
    }
    Ok(())
}

#[test]
fn every_restorable_domain_can_be_restored_again() {
    let failed: Vec<String> = claimants()
        .iter()
        .filter_map(|c| (c.check)().err())
        .collect();
    assert!(failed.is_empty(), "{}", failed.join("\n"));
}

/// **A claim the list has not got fails here.**
///
/// Every `impl Domain for` at the start of a line in `crates/*/src` — which leaves out the test
/// modules' own, indented — and the body its `supports_restore` returns. A body that is neither
/// `true` nor `false` is a claim this cannot read, and fails rather than being guessed at.
#[cfg(not(target_family = "wasm"))]
#[test]
fn every_claim_in_the_library_is_in_the_list() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .expect("crates/pantometry has two ancestors")
        .to_path_buf();
    let mut files = Vec::new();
    let mut dirs: Vec<std::path::PathBuf> = std::fs::read_dir(root.join("crates"))
        .expect("crates/ is readable")
        .map(|e| e.expect("a directory entry").path().join("src"))
        .filter(|p| p.is_dir())
        .collect();
    while let Some(dir) = dirs.pop() {
        for e in std::fs::read_dir(&dir).expect("a source directory is readable") {
            let p = e.expect("a directory entry").path();
            if p.is_dir() {
                dirs.push(p);
            } else if p.extension().is_some_and(|x| x == "rs") {
                files.push(p);
            }
        }
    }

    let mut impls = 0;
    let mut claims = Vec::new();
    let mut unread = Vec::new();
    for f in &files {
        let krate = f
            .strip_prefix(root.join("crates"))
            .expect("under crates/")
            .components()
            .next()
            .expect("a crate directory")
            .as_os_str()
            .to_string_lossy()
            .into_owned();
        let text = std::fs::read_to_string(f).expect("a source file is readable");
        let lines: Vec<&str> = text.lines().collect();
        for (i, line) in lines.iter().enumerate() {
            // Any header at the start of a line that implements the trait, in whatever shape: one
            // this cannot read is a failure rather than a skip, so a generic impl, a `where`
            // clause or a path to the trait cannot take a claim past it.
            if !(line.starts_with("impl") && line.contains("Domain for ")) {
                continue;
            }
            // The kernel's forwarding impl, which claims what the box holds.
            if *line == "impl Domain for Box<dyn Domain> {" {
                continue;
            }
            let name = line
                .strip_prefix("impl Domain for ")
                .and_then(|n| n.strip_suffix(" {"))
                .filter(|n| n.chars().all(|ch| ch.is_alphanumeric() || ch == '_'));
            let Some(name) = name else {
                unread.push(format!("{krate}: a header this cannot read: {line:?}"));
                continue;
            };
            impls += 1;
            let body = lines[i + 1..].iter().take_while(|l| **l != "}");
            let mut it = body.skip_while(|l| !l.trim().starts_with("fn supports_restore"));
            match it.next().map(|l| l.trim()) {
                None => continue,
                Some("fn supports_restore(&self) -> bool {") => {}
                Some(other) => {
                    unread.push(format!(
                        "{krate}::{name}: a signature this cannot read: {other:?}"
                    ));
                    continue;
                }
            }
            match it.next().map(|l| l.trim()) {
                Some("true") => claims.push(format!("{krate}::{name}")),
                Some("false") => {}
                other => unread.push(format!("{krate}::{name} returns {other:?}")),
            }
        }
    }
    // What it could not read first, because a header it cannot read is also one it did not count.
    assert!(unread.is_empty(), "a claim this cannot read: {unread:?}");
    // **A walk that found nothing would pass for ever.** Twenty-four is how many the library has.
    assert!(
        impls >= 24,
        "only {impls} `impl Domain for` found under crates/"
    );

    claims.sort();
    let mut listed: Vec<String> = claimants().into_iter().map(|c| c.ty).collect();
    listed.sort();
    let missing: Vec<&String> = claims.iter().filter(|c| !listed.contains(c)).collect();
    let extra: Vec<&String> = listed.iter().filter(|l| !claims.contains(l)).collect();
    assert!(
        missing.is_empty() && extra.is_empty(),
        "claimed and not tested here: {missing:?}; listed and not claimed: {extra:?}"
    );
    assert_eq!(claims.len(), 18, "the count this file's prose states");
}
