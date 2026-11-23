//! A bar adds the heat it takes off the bus after its sweep, so the step that brings it does not
//! also conduct it.
//!
//! A forward step is `T′ = T + F·∇²T + s`: the sweep reads the state the step began in and what
//! arrives is added as the step's source. `Bar1D` put its heat in first, "so it spreads on the
//! same step it arrives", and that read the bar a sweep ahead of its heat: from a uniform start
//! its state after every step was this order's with one more sweep taken and nothing arriving, to
//! the bit, so a settled bar read every cell one step's arrival short, less its mean. Nothing in
//! the books could tell — the heat stays in the bar, and only moves along it. What moved was the
//! shape the bar settles to, which then depended on its step. Under constant power into the first
//! cell of an insulated bar of 16 cells, the first gap carried 0.893, 0.733 and 0.467 of its share
//! at Fourier numbers of 0.1, 0.25 and 0.5 — `1 − F·n/(n − 1)` exactly — while every other gap
//! carried all of its own.
//!
//! `ThermalNetwork` met this on its first joint and puts its heat into the same right-hand side
//! as its fluxes for it; `Solid3D`'s bus heat met it in a film, and goes in after the sweep.
//!
//! # What is checked against what
//!
//! **Constant power, settled.** An insulated bar under constant power settles to a shape that
//! rises at one rate everywhere, and Fourier's law across each gap then has to carry what the
//! cells beyond it are taking: `kA·(Tᵢ − Tᵢ₊₁)/dx = P₍≤i₎ − (i + 1)·P/n`, where `P₍≤i₎` is the
//! power arriving at or before cell `i`. No step enters that, so the march must land on it at
//! every step up to the limit — for heat with no place, which the bar puts in its first cell, and
//! for heat that knows where it landed, on a profile with no symmetry to hide behind. And under
//! the scheduler, which hands a window's heat to the bar's substeps a share at a time.
//!
//! **A slug.** One step that brings heat is the same step without it with the rise added, in
//! every cell, at the limit: the sweep reads what both steps began in.

use pantometry_core::conserved::quantity;
use pantometry_core::substance::ThermalProps;
use pantometry_core::units::{
    Area, Density, Length, SpecificHeat, Temperature, ThermalConductivity, ThermalExpansion, Time,
};
use pantometry_core::{
    Domain, Exchange, Flux, Kind, Ledger, Schedule, Simulation, Substance, Violation,
};
use pantometry_thermal::{Bar1D, HEAT};

const RHO: f64 = 1000.0;
const C: f64 = 4000.0;
const K: f64 = 50.0;
const DX: f64 = 1e-3;
const AREA: f64 = 1e-4;
/// One cell's heat capacity, from the substance and the geometry rather than from the bar.
const CELL: f64 = RHO * C * AREA * DX;

fn stuff() -> Substance {
    Substance::bulk("stuff", Density::from_si(RHO)).with_thermal(ThermalProps {
        conductivity: ThermalConductivity::w_per_m_k(K),
        specific_heat: SpecificHeat::j_per_kg_k(C),
        expansion: ThermalExpansion::ppm_per_k(50.0),
        emissivity: 0.0,
    })
}

/// A bar of `n` cells at 300 K, exposed so that heat can arrive on its faces.
fn bar(n: usize) -> Bar1D {
    Bar1D::new(
        "bar",
        stuff(),
        n,
        Length::from_si(DX),
        Area::from_si(AREA),
        Temperature::from_si(300.0),
    )
    .exposing("face", Area::from_si(1e-5))
}

/// The step at `fraction` of the explicit limit `dx²/2α`, from the literals.
fn step_at(fraction: f64) -> f64 {
    fraction * DX * DX * RHO * C / (2.0 * K)
}

/// Every gap's settled difference, `Tᵢ − Tᵢ₊₁`, as Fourier's law has it: what crosses gap `i` is
/// what arrives at or before cell `i`, less what those `i + 1` cells keep at the common rate.
fn by_fouriers_law(watts: &[f64]) -> Vec<f64> {
    let n = watts.len();
    let total: f64 = watts.iter().sum();
    let mut before = 0.0;
    (0..n - 1)
        .map(|i| {
            before += watts[i];
            (before - total * (i + 1) as f64 / n as f64) * DX / (K * AREA)
        })
        .collect()
}

/// How many steps at a Fourier number `f` settle a bar of `n` cells: enough for the slowest mode
/// to decay by `e^-118`. It goes as `1 − 4f·sin²(π/2n)` a step, and the highest, at the limit, as
/// `−cos(π/n)`, alternating, so `12n²/f` steps take either past it.
fn settling_steps(n: usize, f: f64) -> usize {
    (12.0 * (n * n) as f64 / f).ceil() as usize
}

/// March a bar under constant power until its shape has settled, and return each gap's
/// difference. `lumped` publishes the total on the plain channel, which the bar puts in its first
/// cell; otherwise each face gets its own. The mean rises the whole time, and that is the point:
/// the shape is what is checked.
fn settle(watts: &[f64], fraction: f64, lumped: bool) -> Vec<f64> {
    let n = watts.len();
    let dt = step_at(fraction);
    let mut bar = bar(n);
    let boundary = bar.boundary().expect("it was exposed").clone();
    let total: f64 = watts.iter().sum();
    let mut bus = Exchange::new();
    for _ in 0..settling_steps(n, fraction / 2.0) {
        if lumped {
            bus.publish(HEAT, total * dt);
        } else {
            let joules: Vec<f64> = watts.iter().map(|w| w * dt).collect();
            bus.publish_on(&boundary, HEAT, &Flux::from_faces(joules))
                .expect("one face per cell");
        }
        bar.step(Time::ZERO, Time::from_si(dt), &mut bus)
            .expect("inside the limit");
    }
    assert!(
        bus.unclaimed().next().is_none(),
        "the bar took everything it was offered"
    );
    gaps(&bar)
}

fn gaps(bar: &Bar1D) -> Vec<f64> {
    (0..bar.cell_count() - 1)
        .map(|i| bar.temperature_at(i).to_si() - bar.temperature_at(i + 1).to_si())
        .collect()
}

/// How far a settled shape is from Fourier's law, at its worst gap — and not a number if any gap
/// is not, where `f64::max` would have passed over it and let a bar of NaN settle at zero.
fn worst(got: &[f64], want: &[f64]) -> f64 {
    got.iter()
        .zip(want)
        .map(|(g, w)| (g - w).abs())
        .fold(0.0, |m, d| if d.is_nan() || d > m { d } else { m })
}

/// The largest residual a settled march leaves, in kelvin.
///
/// A march at 300 K and above holds its temperatures to about `6e-14` and its shape is a balance of
/// differences between them, so what is left is rounding and not the scheme: 8.8e-12 K at the
/// worst of these cases, measured, 33 cells at a tenth of the limit, where the march is longest.
/// The defect this holds out moved a gap by the difference of one step's rise across it: for heat
/// into the first cell, 5e-3 K at a tenth of the limit and 5e-2 K at the limit; on the faces,
/// 5e-4 K and 5e-3 K.
const SETTLED: f64 = 1e-10;

/// **Heat with no place.** The bar puts it in its first cell, and every gap then carries what the
/// cells beyond it take — at a tenth, a quarter, a half and all of the limit, the same shape.
#[test]
fn heat_into_the_first_cell_settles_every_gap_at_its_share() {
    for n in [16usize, 33] {
        let mut watts = vec![0.0; n];
        watts[0] = 0.5;
        let want = by_fouriers_law(&watts);
        for fraction in [0.1, 0.25, 0.5, 1.0] {
            let got = settle(&watts, fraction, true);
            let miss = worst(&got, &want);
            assert!(
                miss < SETTLED,
                "{n} cells at {fraction} of the limit: a gap is {miss:e} K off Fourier's law; the \
                 first carries {:.6} of its share",
                got[0] / want[0]
            );
        }
    }
}

/// **Heat that knows where it landed**, on a profile with no symmetry, so that no gap carries a
/// share a mirror image could have handed it.
#[test]
fn heat_on_the_faces_settles_every_gap_at_its_share() {
    let n = 21;
    let watts: Vec<f64> = (0..n).map(|i| 0.01 * (1 + (3 * i) % 8) as f64).collect();
    let want = by_fouriers_law(&watts);
    // Not a profile any gap of which Fourier's law settles flat: they run from 6.7e-4 K to
    // 7.3e-3 K in size and change sign ten times, against a residual held to 1e-10. A profile
    // that repeats a whole number of times along the bar has gaps carrying nothing, and the first
    // one tried here did.
    assert!(want.iter().all(|w| w.abs() > 1e-4), "{want:?}");
    for fraction in [0.1, 0.5, 1.0] {
        let got = settle(&watts, fraction, false);
        let miss = worst(&got, &want);
        assert!(
            miss < SETTLED,
            "at {fraction} of the limit a gap is {miss:e} K off Fourier's law"
        );
    }
}

/// A heater paying constant power onto the plain channel, a window's worth once a window.
struct Heater {
    watts: f64,
    paid: f64,
}

impl Domain for Heater {
    fn name(&self) -> &str {
        "heater"
    }
    fn kind(&self) -> Kind {
        Kind::QuasiStatic
    }
    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let joules = self.watts * dt.to_si();
        bus.publish(HEAT, joules);
        self.paid += joules;
        Ok(())
    }
    fn ledger(&self) -> Ledger {
        Ledger::new().with(quantity::ENERGY, -self.paid)
    }
    fn checkpoint(&mut self) {}
    fn restore(&mut self) {}
    fn supports_restore(&self) -> bool {
        true
    }
}

/// **Under the scheduler.** The heater publishes a window's heat at once and the bar takes the
/// window in four substeps, each of which must take its share and add it after its own sweep: the
/// first gap then settles where Fourier's law puts it, as it does stepped by hand. A bar taking the
/// whole window on its first substep would conduct three sweeps of it with nothing arriving.
#[test]
fn under_the_scheduler_each_substep_settles_with_its_share() {
    let n = 16;
    let mut watts = vec![0.0; n];
    watts[0] = 0.5;
    let want = by_fouriers_law(&watts);
    // Three and a half of the bar's limit, which it divides into four substeps of 0.875 of it.
    let window = Time::from_si(3.5 * step_at(1.0));
    let mut sim = Simulation::new(Schedule::Multirate)
        .with(Heater {
            watts: 0.5,
            paid: 0.0,
        })
        .with(bar(n));
    for _ in 0..settling_steps(n, 0.875 / 2.0).div_ceil(4) {
        sim.advance(window).expect("the books balance");
    }
    let report = sim.advance(window).expect("the books balance");
    assert!(
        report.substeps.contains(&("bar".to_string(), 4)),
        "the bar was to take the window in four substeps: {:?}",
        report.substeps
    );
    let bar: &Bar1D = sim.domain_as("bar").expect("the bar is there");
    let got = gaps(bar);
    let miss = worst(&got, &want);
    assert!(
        miss < SETTLED,
        "under the scheduler a gap is {miss:e} K off Fourier's law; the first carries {:.6} of its \
         share",
        got[0] / want[0]
    );
}

/// **A slug.** One step at the limit that brings heat, on its plain channel and on its faces, is
/// the step without it with the rise added in every cell.
///
/// From a shape with gradients in it, so that the step without the slug moves every cell: two bars
/// are marched identically under an uneven feed first, and then one is given the slug.
#[test]
fn a_slug_is_the_step_without_it_with_the_rise_added() {
    let n = 21;
    let dt = step_at(1.0);
    let make = || {
        let mut bar = bar(n);
        let boundary = bar.boundary().unwrap().clone();
        let mut bus = Exchange::new();
        for _ in 0..40 {
            let joules: Vec<f64> = (0..n).map(|i| 0.05 * ((5 * i) % 11) as f64).collect();
            bus.publish_on(&boundary, HEAT, &Flux::from_faces(joules))
                .unwrap();
            bar.step(Time::ZERO, Time::from_si(dt), &mut bus).unwrap();
        }
        (bar, boundary)
    };
    let (mut fed, boundary) = make();
    let (mut unfed, _) = make();
    let before: Vec<f64> = (0..n).map(|i| fed.temperature_at(i).to_si()).collect();

    let lumped = 0.4;
    let faces: Vec<f64> = (0..n).map(|i| 0.04 * (1 + (3 * i) % 7) as f64).collect();
    let mut bus = Exchange::new();
    bus.publish(HEAT, lumped);
    bus.publish_on(&boundary, HEAT, &Flux::from_faces(faces.clone()))
        .unwrap();
    fed.step(Time::ZERO, Time::from_si(dt), &mut bus).unwrap();
    unfed
        .step(Time::ZERO, Time::from_si(dt), &mut Exchange::new())
        .unwrap();

    let mut moved = 0.0f64;
    for i in 0..n {
        let rise = (faces[i] + if i == 0 { lumped } else { 0.0 }) / CELL;
        let got = fed.temperature_at(i).to_si() - unfed.temperature_at(i).to_si();
        // The fed cell is the unfed one's value with the rise added — once, or at cell 0 twice,
        // its plain share and then its face's — so the two differ from the rise by those
        // additions' rounding, at most an ulp at 325 K, 5.7e-14; 2.3e-14 is the most measured.
        // The defect conducted `F·∇²` of the rise as well, 0.35 K here at cell 0 and 0.5 K at
        // cell 1.
        assert!(
            (got - rise).abs() < 1e-13,
            "cell {i}: the slug raised it {got} K against {rise} K"
        );
        moved = moved.max((unfed.temperature_at(i).to_si() - before[i]).abs());
    }
    // And the step without the slug did something, so the comparison is not of two bars at rest.
    assert!(
        moved > 0.1,
        "the unfed step moved no cell more than {moved} K"
    );
}
