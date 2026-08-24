//! **Rule 3, for the field stencils: the answer does not depend on the core count.**
//!
//! `CLAUDE.md` says results are bit-for-bit identical across platforms, optimisation levels,
//! WebAssembly and thread counts. Until now the last clause was free — the stencils had no
//! threads. It is not free any more, and this is what holds it.
//!
//! Compared **bit for bit**, not to a tolerance. A tolerance would pass a version that reordered
//! the arithmetic, and reordering is exactly the failure a parallel sweep can introduce: the
//! stencil is safe because each cell reads a snapshot, and a *sum* is not, which is why `lost` is
//! folded in a sequential pass afterwards.
//!
//! The thread count comes from `PANTOMETRY_THREADS`. It exists for this test: the only way to
//! check "the answer is the same at any count" on one machine is to run it at two counts.

use pantometry_core::units::{Energy, Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Substance};
use pantometry_thermal::Solid3D;

/// Big enough that `sweep::fill` really splits — the floor is 8,192 cells a thread, and 40³ is
/// 64,000. A grid that stayed sequential whatever the setting would pass this test by not
/// testing it.
const N: usize = 40;

/// Every number the run produces: the grid, and the scalars that are sums over it.
fn run_at(threads: usize, steps: usize) -> (Vec<u64>, Vec<u64>) {
    std::env::set_var("PANTOMETRY_THREADS", threads.to_string());

    let mut block = Solid3D::new(
        "block",
        Substance::aluminium_6061(),
        (N, N, N),
        Length::from_si(1e-3),
        Temperature::celsius(20.0),
    );
    // Off centre and asymmetric, so a sweep that transposed an axis or mirrored a chunk would
    // produce a different field rather than the same one by symmetry.
    block.deposit(N / 3, N / 2 + 1, 2 * N / 5, Energy::from_si(4.0));
    block.deposit(1, 1, 1, Energy::from_si(1.5));

    let dt = Time::from_si(block.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
    let mut bus = Exchange::new();
    let mut t = 0.0;
    for _ in 0..steps {
        block.step(Time::from_si(t), dt, &mut bus).expect("stable");
        t += dt.to_si();
    }

    let field = (0..N * N * N)
        .map(|c| {
            let (i, j, k) = (c % N, (c / N) % N, c / (N * N));
            block.temperature_at(i, j, k).to_si().to_bits()
        })
        .collect();
    let readings = block.readings().iter().map(|r| r.value.to_bits()).collect();
    (field, readings)
}

/// **One thread and many produce the same bits, in the grid and in the scalars.**
#[test]
fn a_threaded_sweep_is_the_sequential_one() {
    let (one_field, one_readings) = run_at(1, 60);

    for threads in [2usize, 3, 5, 8, 16] {
        let (field, readings) = run_at(threads, 60);
        // Report where, not just that: a single differing cell is a chunk boundary and a
        // scattered thousand is a race.
        let differing: Vec<usize> = (0..field.len())
            .filter(|&c| field[c] != one_field[c])
            .collect();
        assert!(
            differing.is_empty(),
            "{threads} threads changed {} of {} cells; the first is {:?}",
            differing.len(),
            field.len(),
            differing
                .first()
                .map(|&c| (c % N, (c / N) % N, c / (N * N)))
        );
        assert_eq!(
            readings, one_readings,
            "{threads} threads changed a scalar — a sum folded out of index order"
        );
    }

    std::env::remove_var("PANTOMETRY_THREADS");
}

/// **And a block that sheds to air**, which is the case where the sweep produces a *second*
/// per-cell number that a sequential pass has to total.
///
/// `lost` is the accumulator: joules the surface film took out. Summed inside the parallel pass it
/// would depend on how the work divided; summed afterwards in index order it does not. This is the
/// test that would fail if somebody moved that addition back into the sweep.
#[test]
fn the_shed_joules_total_the_same_however_it_was_split() {
    let build = |threads: usize| {
        std::env::set_var("PANTOMETRY_THREADS", threads.to_string());
        let mut block = Solid3D::new(
            "block",
            Substance::aluminium_6061(),
            (N, N, N),
            Length::from_si(1e-3),
            Temperature::celsius(200.0),
        )
        .losing_from(
            pantometry_thermal::Face::ZMax,
            pantometry_thermal::Environment::still_air(
                Temperature::celsius(20.0),
                pantometry_core::units::Area::from_si(1e-3 * 1e-3),
            ),
        );
        let dt = Time::from_si(block.max_stable_dt(Time::from_si(0.0)).to_si() * 0.5);
        let mut bus = Exchange::new();
        for _ in 0..60 {
            block
                .step(Time::from_si(0.0), dt, &mut bus)
                .expect("stable");
        }
        block
            .readings()
            .iter()
            .map(|r| (r.label.clone(), r.value.to_bits()))
            .collect::<Vec<_>>()
    };

    let one = build(1);
    for threads in [2usize, 4, 9] {
        assert_eq!(
            build(threads),
            one,
            "{threads} threads changed what a cooled block reports"
        );
    }
    std::env::remove_var("PANTOMETRY_THREADS");
}
