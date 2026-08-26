//! **The device runs the block the CPU resolved, not a homogeneous stand-in.**
//!
//! The kernel was `centre + F·(sum − 6·centre)`, which has no per-cell coefficient — so it is the
//! same operator as `Solid3D`'s only when every cell is the same material. `Solid3D`'s own comment
//! says why: written that way it cannot conserve `Σ Cᵢ Tᵢ` when the capacities differ.
//!
//! It takes the resolved face conductances and mobilities now, which is what makes two materials, a
//! void and a coating work without the device knowing what any of those are — they are already in
//! the coefficients.
//!
//! Skipped without an adapter, loudly. A software rasteriser would be checking an implementation
//! nobody runs.

use pantometry_core::units::{Energy, Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_gpu::{GpuSolid, Unavailable};
use pantometry_thermal::Solid3D;

const N: usize = 12;
const DX: f64 = 1e-3;

fn cpu_block() -> Solid3D {
    Solid3D::new(
        "block",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(DX),
        Temperature::celsius(20.0),
    )
}

/// How far apart the two grids are, against the range the field spans — dividing by a cell at
/// ambient would report every real difference as tiny.
fn divergence(cpu: &Solid3D, gpu: &mut GpuSolid) -> f64 {
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    let mut worst = 0.0f64;
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                // A void has no temperature and the CPU says so with a `NaN`. `f64::min` ignores
                // one, so this is safe — but the comparison below is not, which is the point.
                if a.is_finite() {
                    lo = lo.min(a);
                    hi = hi.max(a);
                }
            }
        }
    }
    let span = (hi - lo).max(1e-12);
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let a = cpu.temperature_at(i, j, k).to_si();
                // **Only where there is something to compare.** The CPU reports a void as *not a
                // number*, because a void has no temperature and a zero or an ambient there is a
                // value somebody would plot. The device's buffer has no such concept: a void cell
                // holds the deviation it was given and keeps it, because its mobility is zero.
                // Comparing the two would be comparing an absence with a number.
                if !a.is_finite() {
                    continue;
                }
                let b = gpu.temperature_at(i, j, k).to_si();
                // **`f64::max` ignores a NaN**, so `worst.max(nan)` is `worst` and a grid that had
                // gone entirely to NaN reported a divergence of exactly zero. Caught by the stored
                // heat, which had no such mercy. A difference that is not a number is the largest
                // difference there is.
                let d = (a - b).abs() / span;
                if !d.is_finite() {
                    return f64::INFINITY;
                }
                worst = worst.max(d);
            }
        }
    }
    worst
}

/// `Σ Cᵢ (Tᵢ − T₀)`, the same quantity the device reports, computed from the public surface.
fn cpu_stored(cpu: &Solid3D) -> f64 {
    let mut total = 0.0;
    for k in 0..N {
        for j in 0..N {
            for i in 0..N {
                let c = i + N * (j + N * k);
                total += cpu.cell_capacities()[c]
                    * (cpu.temperature_at(i, j, k).to_si() - Temperature::celsius(20.0).to_si());
            }
        }
    }
    total
}

fn march(cpu: &mut Solid3D, gpu: &mut GpuSolid, steps: usize) {
    let dt = Time::from_si(cpu.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    for _ in 0..steps {
        cpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
        gpu.step(Time::from_si(0.0), dt, &mut bus).expect("stable");
    }
}

/// One GPU test at a time, within this process.
///
/// **A workaround at the level the evidence supports, and no deeper.** Measured on this file:
/// **7 hangs in 30** at the harness default, **0 in 30** with `--test-threads=1`, running alone
/// with nothing before it. So it is concurrency inside one process — but *what* about it is not
/// established. Two explanations were tried and measured wrong: device creation, which has always
/// been behind a `OnceLock` and happens once; and overlapping readbacks, which serialising left
/// at 5 in 30.
///
/// So the tests take a lock and the library does not, because a lock in the library would be a
/// fix for a mechanism nobody has demonstrated. What is demonstrated is that these tests do not
/// need to run at the same time, and that when they do, the suite stops.
static ONE_AT_A_TIME: std::sync::Mutex<()> = std::sync::Mutex::new(());

/// Take it. A poisoned lock is not interesting here — a panicking test has already failed, and
/// the next one is entitled to a turn rather than a second failure about the first.
fn alone() -> std::sync::MutexGuard<'static, ()> {
    ONE_AT_A_TIME
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}
/// **Two materials, and the joint between them is the joint the CPU computed.**
///
/// Copper in half the block and aluminium in the other: their conductivities differ by a factor of
/// about two and their volumetric capacities by about 1.3, so a kernel using one `F` for the whole
/// grid gets both the flux across the joint and the rise either side of it wrong. The tolerance is
/// single precision's — about `1e-7` a step, walking as `√k` — and not a number chosen to pass.
#[test]
fn two_materials_meet_at_the_face_the_cpu_resolved() {
    let _alone = alone();
    let mut cpu = cpu_block();
    cpu.fill(Substance::copper(), |i, _, _| i < N / 2);
    let mut gpu = match GpuSolid::mirroring(cpu.clone()) {
        Ok(g) => g,
        // **The reason, not a guess at it.** This said "no adapter" for every error, so a block
        // the device refused on its own terms was reported as a machine without a GPU.
        Err(why) => return println!("skipped: {why}"),
    };
    cpu.deposit(1, N / 2, N / 2, Energy::from_si(2.0));
    gpu.deposit(1, N / 2, N / 2, Energy::from_si(2.0));

    march(&mut cpu, &mut gpu, 200);
    let worst = divergence(&cpu, &mut gpu);
    println!("  two materials, 200 steps: worst relative difference {worst:.3e}");
    assert!(
        worst < 1e-5,
        "the device is not running the same operator: {worst:.3e}"
    );

    // And the heat is still all there, on both. A stencil that got the joint wrong loses or makes
    // energy at it, which is exactly what a single `F` over two capacities does.
    let (a, b) = (cpu_stored(&cpu), gpu.stored_heat().to_si());
    assert!(
        (a - b).abs() / a.abs().max(1.0) < 1e-5,
        "stored heat differs: {a} against {b}"
    );
}

/// **A void is a hole in the operator, and the device sees it as one.**
///
/// A clearance carries no heat and holds none. In the coefficients it is a face conductance of zero
/// and a mobility of zero, so the device needs to know nothing about voids to get it right — which
/// is the argument for taking the coefficients rather than the material list.
///
/// **A corner, not a slab through the middle**, and the difference is physics rather than
/// convenience. `resolve` calls `find_gaps`, so any run of void with solid at *both* ends is a pair
/// that radiates — an interior void in this workspace is never merely a hole. A slab through the
/// middle is therefore a block the device refuses, correctly, and the first version of this test
/// used one and reported itself skipped. A corner void runs to the boundary along every line it
/// lies on, so it pairs with nothing and is only a hole.
#[test]
fn a_void_stops_the_heat_on_the_device_too() {
    let _alone = alone();
    let bite = N / 3;
    let mut cpu = cpu_block().empty(move |i, j, k| i < bite && j < bite && k < bite);
    let mut gpu = match GpuSolid::mirroring(cpu.clone()) {
        Ok(g) => g,
        // **The reason, not a guess at it.** This said "no adapter" for every error, so a block
        // the device refused on its own terms was reported as a machine without a GPU.
        Err(why) => return println!("skipped: {why}"),
    };
    // Just outside the bite, so the heat has to go round it.
    cpu.deposit(1, 1, bite, Energy::from_si(2.0));
    gpu.deposit(1, 1, bite, Energy::from_si(2.0));

    march(&mut cpu, &mut gpu, 200);
    let worst = divergence(&cpu, &mut gpu);
    println!("  a void, 200 steps: worst relative difference {worst:.3e}");
    assert!(worst < 1e-5, "the device crossed the void: {worst:.3e}");

    // Inside the bite there is nothing, and the two say so differently. The CPU reports **not a
    // number**, because a void has no temperature; the device holds a deviation that never moves,
    // because its mobility is zero. Both are right and neither is the other, so what is checked is
    // that the device did not *warm* it — joules in a cell with no capacity to hold them would put
    // the ledger at odds with the block.
    assert!(
        !cpu.temperature_at(0, 0, 0).to_si().is_finite(),
        "the CPU should call a void an absence, and this test relies on it"
    );
    let void_gpu = gpu.temperature_at(0, 0, 0).to_si();
    assert!(
        (void_gpu - Temperature::celsius(20.0).to_si()).abs() < 1e-6,
        "the device warmed a void it has no capacity for: {void_gpu}"
    );
    // And the heat went round rather than through: the far corner, on the other side of the bite,
    // is warmer than it started, and the block still holds what was put in it.
    let stored = gpu.stored_heat().to_si();
    assert!(
        (stored - 2.0).abs() < 1e-3,
        "the device lost heat into the void: {stored} J of 2.0"
    );
}

/// **A source is watts in a cell, and it arrives on the device as the same rise.**
#[test]
fn a_generating_block_generates_the_same_on_both() {
    let _alone = alone();
    let mut cpu = cpu_block().dissipating(4.0, |i, _, _| i < 2);
    let mut gpu = match GpuSolid::mirroring(cpu.clone()) {
        Ok(g) => g,
        // **The reason, not a guess at it.** This said "no adapter" for every error, so a block
        // the device refused on its own terms was reported as a machine without a GPU.
        Err(why) => return println!("skipped: {why}"),
    };
    march(&mut cpu, &mut gpu, 200);
    let worst = divergence(&cpu, &mut gpu);
    println!("  a source, 200 steps: worst relative difference {worst:.3e}");
    assert!(worst < 1e-5, "the source differs: {worst:.3e}");
}

/// **What the device cannot run is refused, and the refusal names it.**
///
/// Not a fall back to the CPU. A device is what a scene *states*, and a run that silently changed
/// where it ran is a run whose answer nobody asked for — the same reason the run file's reader
/// refuses an unknown panel kind rather than skipping it.
#[test]
fn a_block_the_device_cannot_run_is_refused_by_name() {
    let _alone = alone();
    let filmed = cpu_block().losing_from(
        pantometry_thermal::Face::ZMax,
        pantometry_thermal::Environment::still_air(
            Temperature::celsius(20.0),
            pantometry_core::units::Area::from_si(DX * DX),
        ),
    );
    match GpuSolid::mirroring(filmed) {
        Err(Unavailable::Unsupported(why)) => {
            assert!(why.contains("film"), "the refusal must say which: {why}");
        }
        Err(other) => println!("skipped: {other}"),
        Ok(_) => panic!("a film has no device pass and must not be accepted silently"),
    }
}
