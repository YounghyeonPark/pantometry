//! The GPU against the domain it accelerates.
//!
//! Not "does it run" — every wrong stencil runs. The question is how far a single-precision port
//! lands from the `f64` domain that is the reference, and whether that distance is the shape of
//! rounding or the shape of a bug.
//!
//! # These tests skip when there is no GPU, and say so
//!
//! A CI runner usually has no adapter. Skipping is the honest outcome — the alternative is a
//! software rasteriser, which would be checking a different implementation than anyone runs. What
//! is **not** acceptable is skipping quietly, so every skip prints why.

use pantometry_core::units::{Energy, Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_gpu::GpuSolid;
use pantometry_thermal::Solid3D;

const N: usize = 16;
const DX: f64 = 1e-3;

/// Build the pair. `None` when this machine has no GPU, with a reason on stdout.
fn pair() -> Option<(Solid3D, GpuSolid)> {
    let cpu = Solid3D::new(
        "cpu",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    match GpuSolid::new(
        "gpu",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    ) {
        Ok(gpu) => Some((cpu, gpu)),
        Err(why) => {
            println!("skipped: {why}. Nothing here can run without one, and a software adapter");
            println!("         would be checking a different implementation than anyone uses.");
            None
        }
    }
}

/// Deposit the same joules in the same cell on both.
fn seed(cpu: &mut Solid3D, gpu: &mut GpuSolid) {
    cpu.deposit(N / 2, N / 2, N / 2, Energy::from_si(2.0));
    gpu.deposit(N / 2, N / 2, N / 2, Energy::from_si(2.0));
}

fn run(cpu: &mut Solid3D, gpu: &mut GpuSolid, steps: usize, dt: Time) {
    let mut bus = Exchange::new();
    let mut t = 0.0;
    for _ in 0..steps {
        cpu.step(Time::from_si(t), dt, &mut bus).expect("stable");
        gpu.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }
}

/// The largest relative difference between the two grids, measured against the range the field
/// actually spans rather than against each cell — a cell at ambient is 293 K and dividing by it
/// would report every real difference as tiny.
fn divergence(cpu: &Solid3D, gpu: &mut GpuSolid) -> f64 {
    let ambient = Temperature::celsius(20.0).to_si();
    let mut worst = 0.0f64;
    let mut scale: f64 = 1e-30;
    let cells = gpu.cells();
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                scale = scale.max((a - ambient).abs());
            }
        }
    }
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                let b = cells[i + N * (j + N * k)];
                worst = worst.max((a - b).abs() / scale);
            }
        }
    }
    worst
}

/// **The GPU reproduces the reference to single precision, and the figure is reported.**
///
/// The claim is not that they agree. WGSL has no `f64`, so they cannot: the port is a
/// lower-precision arithmetic and the only useful question is how much lower.
///
/// The tolerance is earned rather than tried. `f32` has about 7 decimal digits, so one update —
/// a sum of seven terms and two multiplies — loses of order `1e-7` relative. Over `k` steps the
/// error walks as roughly `√k`, so 200 steps is about `1.4e-6`. Anything at that order is
/// rounding; anything above it by orders is a different stencil.
#[test]
fn the_gpu_reproduces_the_cpu_to_single_precision() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    // Sixty, not two hundred. At half the stability limit the diffusion length after `k` steps is
    // `dx·√(k/3)`, so 200 steps reaches 8 mm on a 16 mm block — the spot has hit the walls and
    // levelled, and comparing two nearly-uniform grids is agreement for free. Sixty puts it at
    // 4.5 mm, a quarter of the block, where there is still structure to disagree about.
    let steps = 60;
    run(&mut cpu, &mut gpu, steps, dt);
    let worst = divergence(&cpu, &mut gpu);
    let expected = 1e-7 * (steps as f64).sqrt();

    println!("  after {steps} steps: worst relative difference {worst:.3e}");
    println!("  single precision over that many steps predicts about {expected:.3e}");
    assert!(
        worst < 20.0 * expected,
        "the two have diverged by more than rounding: {worst:.3e} against {expected:.3e}"
    );
    // And they are genuinely both computing something, or agreeing is free. Measured against the
    // deposit rather than against a round number: 2 J into one cell of this block is an 827 K
    // rise, and after 200 steps it has spread but is nowhere near level.
    let spread = cpu.peak_temperature().to_si() - cpu.coldest_temperature().to_si();
    let uniform = 2.0
        / Substance::aluminium_6061()
            .heat_capacity(cpu.volume())
            .unwrap()
            .to_si();
    println!(
        "  the spot is still {:.0}x the levelled rise",
        spread / uniform
    );
    assert!(
        spread > 5.0 * uniform,
        "the spot must still be a spot: {spread:.4} K against {uniform:.4} K levelled"
    );
}

/// **The divergence is the shape of rounding, not of a bug.**
///
/// A wrong stencil — a missing arm, a swapped axis, a mirror that is a zero — does not drift; it
/// is wrong immediately and stays wrong. Rounding grows slowly. So the check is the *shape*: after
/// ten steps the two must be far closer than after two hundred, and both far below what a real
/// disagreement would give.
#[test]
fn the_difference_grows_like_rounding_rather_than_appearing_at_once() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    run(&mut cpu, &mut gpu, 10, dt);
    let early = divergence(&cpu, &mut gpu);
    run(&mut cpu, &mut gpu, 50, dt);
    let late = divergence(&cpu, &mut gpu);

    println!("  10 steps {early:.3e}   60 steps {late:.3e}");
    assert!(
        late > early,
        "rounding accumulates; a wrong stencil would be wrong at step one"
    );
    assert!(
        early < 1e-6,
        "ten steps of f32 should be near-exact, was {early:.3e}"
    );
}

/// **Energy is conserved on the GPU too, to what single precision can hold.**
///
/// The faces are insulated and nothing is on the bus, so the mean is fixed. In `f64` the CPU holds
/// that to about `1e-12`; `f32` cannot, and the gap is the honest cost of the port.
///
/// This is why `GpuSolid` declines `books_balance` and why a scene using it needs
/// `conservation_tolerance_for(ENERGY, ..)` loosened: `Simulation`'s default `1e-9` is below what
/// the arithmetic can deliver, and a run would be refused for being single precision rather than
/// for being wrong.
#[test]
fn the_gpu_conserves_to_what_f32_can_hold() {
    let Some((mut cpu, mut gpu)) = pair() else {
        return;
    };
    seed(&mut cpu, &mut gpu);
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);

    let before_cpu = cpu.mean_temperature().to_si();
    let before_gpu = gpu.mean_temperature().to_si();
    run(&mut cpu, &mut gpu, 60, dt);
    let cpu_drift = (cpu.mean_temperature().to_si() - before_cpu).abs() / before_cpu;
    let gpu_drift = (gpu.mean_temperature().to_si() - before_gpu).abs() / before_gpu;

    println!("  mean drift: cpu {cpu_drift:.3e}   gpu {gpu_drift:.3e}");
    assert!(
        cpu_drift < 1e-12,
        "the f64 reference is exact: {cpu_drift:.3e}"
    );
    assert!(
        gpu_drift < 1e-6,
        "even f32 should conserve to a millionth: {gpu_drift:.3e}"
    );
    // The point of the pair: the GPU is measurably looser, and by how much is the finding.
    println!(
        "  so the accelerator is about {:.0}x looser than the domain it accelerates",
        (gpu_drift / cpu_drift.max(1e-16)).max(1.0)
    );
}

/// **The stability limit is refused on the GPU exactly as on the CPU.**
#[test]
fn past_the_limit_is_refused() {
    let Some((cpu, mut gpu)) = pair() else {
        return;
    };
    let limit = cpu.max_stable_dt(Time::from_si(0.0));
    assert!(
        (gpu.max_stable_dt(Time::from_si(0.0)).to_si() / limit.to_si() - 1.0).abs() < 1e-12,
        "the two report the same limit; the scheme is the same and only the precision differs"
    );
    let err = gpu
        .step(
            Time::from_si(0.0),
            Time::from_si(limit.to_si() * 1.05),
            &mut Exchange::new(),
        )
        .expect_err("5% past the limit must be refused");
    assert_eq!(err.quantity, "Fourier number");
}

/// **It is actually faster, and the README's table is this test's output.**
///
/// Reported rather than asserted. A timing threshold on somebody else's machine fails for reasons
/// that have nothing to do with the code, and a GPU slower than a CPU on a small grid is a true fact
/// about small grids rather than a defect.
///
/// # One size per process, and that is the whole finding
///
/// The README claimed **191× at 64³** and `ARCHITECTURE.md` claimed it twice, and nothing in the
/// repository measured it: this test ran at one hard-coded size and printed one row.
///
/// Replacing it with a sweep over seven sizes did not work, and the way it failed is worth keeping.
/// This machine slows under sustained load — **both columns together**, so it is neither the device
/// nor the allocator — and the drift is seconds fast. Whatever size went last was penalised, by a
/// factor of two to three: run backwards, 128³ read **42×** where forwards it read **21×**, and 16³
/// went from `2.2e-5` to `4.7e-5` s a step for identical work. Best-of-three did not help, because
/// all three reps of a row sit at the same place in the sweep. Round-robining the reps did not help
/// either; it only made the spread column honest about how large the drift is — 100% to 400%.
///
/// What *is* reproducible is the first thing measured in a fresh process: to ±5% across runs, every
/// time. So the size comes from `PANTOMETRY_SWEEP` and the README loops the shell over it. Each
/// number is then a first measurement, and rows compare because none of them paid for the others.
#[test]
fn how_much_faster() {
    let n: usize = std::env::var("PANTOMETRY_SWEEP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(64);
    // Fewer steps on the big grids: 128³ at 400 is minutes of wall clock to sharpen a ratio that is
    // clear at 100.
    let steps = if n > 64 { 100 } else { 400 };

    // Best of three of *this* size, back to back. Same-size reps sit close enough together that the
    // drift inside them is a few per cent, which the spread column reports either way.
    const REPS: usize = 3;
    let mut cpu = Spread::new();
    let mut gpu = Spread::new();
    let mut where_it_ran = String::new();
    for _ in 0..REPS {
        let Some((cpu_time, gpu_time, on)) = speed_at(n, steps) else {
            return;
        };
        cpu.saw(cpu_time);
        gpu.saw(gpu_time);
        where_it_ran = on;
    }

    println!(
        "  | {n}³ | {steps} | {} | {:.3e} | {:.3e} | **{}** | {:.0}% / {:.0}% |",
        // **What the CPU column actually used**, not what the machine has: `threads_for` is
        // `√(cells/39000)` capped at the cores, so 48³ and under run on one thread however many
        // the machine has. A "32 threads" claim was in the first draft of the README's table.
        pantometry_core::sweep::threads_for(n * n * n),
        cpu.best / steps as f64,
        gpu.best / steps as f64,
        // **Two significant figures below ten.** `{:.0}` printed the 16³ ratio of 0.69 as `1×`,
        // which is a table rounding a loss into a tie.
        ratio(cpu.best / gpu.best.max(1e-9)),
        cpu.spread_percent(),
        gpu.spread_percent(),
    );
    println!("  on {where_it_ran}");
    println!("  PANTOMETRY_SWEEP picks the grid; a readback is one transfer, so a run that audits");
    println!("  every step pays for it");
}

/// Wall clock for `steps` on both, at `n³`, and what the device calls itself. `None` when this
/// machine has no adapter.
fn speed_at(n: usize, steps: usize) -> Option<(f64, f64, String)> {
    let mut cpu = Solid3D::new(
        "cpu",
        Substance::aluminium_6061(),
        (n, n, n),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    );
    let mut gpu = match GpuSolid::new(
        "gpu",
        Substance::aluminium_6061(),
        (n, n, n),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    ) {
        Ok(gpu) => gpu,
        Err(why) => {
            println!("  skipped at {n}³: {why}");
            return None;
        }
    };
    cpu.deposit(n / 2, n / 2, n / 2, Energy::from_si(2.0));
    gpu.deposit(n / 2, n / 2, n / 2, Energy::from_si(2.0));

    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();

    let start = std::time::Instant::now();
    for _ in 0..steps {
        cpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
    }
    let cpu_time = start.elapsed().as_secs_f64();

    let start = std::time::Instant::now();
    for _ in 0..steps {
        gpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
    }
    // Force the queue to drain, or this times the submission and not the work.
    let _ = gpu.mean_temperature();
    let gpu_time = start.elapsed().as_secs_f64();
    Some((cpu_time, gpu_time, gpu.device_name().to_string()))
}

/// Best and worst of a few timings of the same thing.
struct Spread {
    best: f64,
    worst: f64,
}

impl Spread {
    fn new() -> Spread {
        Spread {
            best: f64::MAX,
            worst: 0.0,
        }
    }

    fn saw(&mut self, t: f64) {
        self.best = self.best.min(t);
        self.worst = self.worst.max(t);
    }

    /// How much the slowest run exceeded the fastest. Printed because a column with a 100% spread
    /// is a column whose single-shot number means nothing — which is how a 96³ row once came out
    /// *below* the 64³ one.
    fn spread_percent(&self) -> f64 {
        100.0 * (self.worst - self.best) / self.best.max(1e-12)
    }
}

/// Two significant figures below ten, none above: `0.69×`, `4.8×`, `87×`.
fn ratio(r: f64) -> String {
    if r < 10.0 {
        format!("{r:.2}×")
    } else {
        format!("{r:.0}×")
    }
}
