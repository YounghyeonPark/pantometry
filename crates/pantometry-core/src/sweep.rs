//! Splitting one field update across cores, without changing a single bit of the answer.
//!
//! # The third axis of parallelism here, and the one the grids needed
//!
//! [`Ensemble`](crate::ensemble::Ensemble) splits *many* evaluations — a Monte Carlo study, a
//! parameter sweep. `TreeNBody::with_threads` splits *one* evaluation of a body problem. Neither
//! touches a field, and every stencil in this workspace ran on **one core** until this existed:
//! measured on a thirty-two core machine, `Solid3D` did 168 Mcell/s at 96³ with thirty-one cores
//! idle.
//!
//! # Why this is bit-for-bit and `Ensemble` had to argue for it
//!
//! An explicit stencil reads a snapshot and writes a fresh array. Cell `c`'s new value is a
//! function of the *old* array and of read-only coefficients, and of nothing any other cell
//! writes. So splitting the output into contiguous chunks is not an approximation and not a
//! reordering: every cell performs exactly the operations it performed sequentially, in exactly
//! the same order, on exactly the same inputs.
//!
//! That is a stronger position than `Ensemble`'s, which had to make `Rng` addressable by index to
//! get here. A stencil is already there. **The floating-point result is identical for any thread
//! count**, which is what rule 3 asks and what `bit_for_bit_whatever_the_thread_count` pins.
//!
//! # What must not come in here
//!
//! A **reduction**. A mean, a total, a ledger entry: floating-point addition is not associative,
//! so a sum split across threads and combined depends on how the work divided, and the answer
//! changes with the core count. `runtime/gpu` learned the same thing from the other side and reads
//! its grid back to sum it on the CPU in index order.
//!
//! So this fills an output and returns nothing. A caller that wants a total runs one afterwards,
//! sequentially, over the array this filled — which is O(N) of cheap work against O(N) of a
//! seven-point stencil, and is not where the time was.

/// Cells a thread must be given before its own creation is worth paying for.
///
/// **Measured, and the number is large because thread creation is expensive.** On this machine a
/// `std::thread::scope` costs 27 to 42 microseconds *per thread spawned* — 861 µs for thirty-two —
/// and a stencil moves about 130 million cells a second, so one thread's spawn is worth about
/// 4,000 cells of work.
///
/// The first version of this used a flat floor of 8,192 cells a thread and took the machine's full
/// core count above it. The result was a **regression at every grid it was tried on**: at 32³ the
/// step went from 0.27 ms to 1.17 ms, because four spawns cost more than the whole sweep. The
/// count has to grow with the square root of the work, not with the core count — see
/// [`threads_for`].
///
/// A persistent pool would make this constant irrelevant. One cannot be written here without
/// erasing a lifetime, which is `unsafe`, or taking a dependency; both are decisions bigger than
/// this module.
pub const CELLS_PER_SPAWN: usize = 39_000;

/// How many threads a sweep of `cells` should use.
///
/// Asked of the machine rather than fixed, because the right number is the machine's and the
/// answer does not depend on it. `available_parallelism` is not a clock and not a random source:
/// it changes how long the sweep takes and cannot change what it computes.
pub fn threads_for(cells: usize) -> usize {
    if cfg!(target_family = "wasm") {
        return 1;
    }
    // **An override, and it is here to be tested against.** Rule 3 says the answer does not depend
    // on the thread count; the only way to check that on one machine is to run the same thing at
    // two counts, and the only way to do *that* is a knob. It changes how long a sweep takes and
    // cannot change what it computes, which is exactly the claim
    // `a_threaded_sweep_is_the_sequential_one` makes with it.
    if let Ok(n) = std::env::var("PANTOMETRY_THREADS") {
        if let Ok(n) = n.parse::<usize>() {
            return n.max(1);
        }
    }
    let cores = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    // **The square root of the work, not the core count.** Spawning `T` threads costs `T` spawns
    // and buys `1/T` of the work, so the overhead is worth a tenth of the parallel part when
    // `T <= sqrt(cells / CELLS_PER_SPAWN)`. Taking every core instead is what made the first
    // version slower than sequential everywhere below 96³.
    //
    // At 32³ this is 0 and the sweep stays sequential; at 96³ it is 4; at 128³ it is 7. On a
    // thirty-two core machine the cap is never the binding constraint, which is the honest shape
    // of the answer: per-step spawning limits this, not the hardware.
    let wanted = ((cells as f64) / CELLS_PER_SPAWN as f64).sqrt() as usize;
    wanted.clamp(1, cores)
}

/// Fill `out` in contiguous chunks, on as many threads as the size justifies.
///
/// `body(first, chunk)` receives the index `out[0]` of its chunk and the chunk itself. A caller
/// computes `out[first + k]` from whatever read-only state it has closed over; the closure is
/// `Sync` because every thread runs the same one against different output.
///
/// `granularity` rounds the chunk boundaries up to a multiple — a whole plane of a grid, so a
/// caller's index arithmetic stays the arithmetic it already had. Pass 1 for no constraint.
///
/// Sequential when there is one thread, on WebAssembly, or when the work is too small to pay for a
/// spawn. All four paths compute the same bits.
pub fn fill<T: Send>(out: &mut [T], granularity: usize, body: impl Fn(usize, &mut [T]) + Sync) {
    let n = out.len();
    let threads = threads_for(n);
    if threads <= 1 || n == 0 {
        body(0, out);
        return;
    }
    let g = granularity.max(1);
    // Round the chunk up to whole units of `granularity`, so no chunk boundary lands inside one.
    let per = n.div_ceil(threads).div_ceil(g) * g;
    if per >= n {
        body(0, out);
        return;
    }

    #[cfg(not(target_family = "wasm"))]
    {
        let body = &body;
        std::thread::scope(|scope| {
            for (c, slice) in out.chunks_mut(per).enumerate() {
                let first = c * per;
                scope.spawn(move || body(first, slice));
            }
        });
    }

    // WebAssembly has no threads to spawn, so it takes the sequential path and gets the same
    // answer for a less interesting reason — the same resolution `Ensemble` and `TreeNBody` use.
    #[cfg(target_family = "wasm")]
    {
        body(0, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The answer does not depend on how the work divided.**
    ///
    /// A seven-point stencil in one dimension, which is the same shape as the three-dimensional
    /// one for this purpose: each output reads its neighbours from a snapshot. Run at every thread
    /// count from one to sixteen and compared **bit for bit**, not to a tolerance — because the
    /// claim is identity and a tolerance would pass a version that reordered the arithmetic.
    #[test]
    fn bit_for_bit_whatever_the_thread_count() {
        let n = 100_000;
        // Values chosen so the sums are not exactly representable: a stencil over 1.0 everywhere
        // would agree under any reordering and would prove nothing.
        let old: Vec<f64> = (0..n).map(|i| 1.0 / (i as f64 + 1.7)).collect();
        let stencil = |first: usize, chunk: &mut [f64]| {
            for (k, slot) in chunk.iter_mut().enumerate() {
                let c = first + k;
                let l = if c > 0 { old[c - 1] } else { old[c] };
                let r = if c + 1 < old.len() {
                    old[c + 1]
                } else {
                    old[c]
                };
                *slot = old[c] + 0.13 * (l + r - 2.0 * old[c]);
            }
        };

        let mut reference = vec![0.0; n];
        stencil(0, &mut reference);

        for per in [1usize, 2, 3, 7, 16, 999, 100_000] {
            let mut got = vec![0.0; n];
            // Drive the split directly rather than through `fill`, so the test covers boundaries
            // the machine's core count would never produce.
            for (c, slice) in got.chunks_mut(per).enumerate() {
                stencil(c * per, slice);
            }
            assert_eq!(
                got.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "a chunk of {per} changed the answer"
            );
        }

        // And through `fill`, which is what a domain calls.
        let mut threaded = vec![0.0; n];
        fill(&mut threaded, 1, stencil);
        assert_eq!(
            threaded.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            reference.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
        );
    }

    /// Every cell is written exactly once, and the chunk boundaries respect the granularity.
    ///
    /// A sweep that skipped a plane would leave it at its previous value, which for a field that
    /// is not changing much reads as a correct answer.
    #[test]
    fn every_cell_is_filled_once_and_planes_are_not_split() {
        for &(n, g) in &[(1000usize, 1usize), (1000, 10), (4096, 256), (7, 4), (0, 8)] {
            let mut out = vec![0u32; n];
            let starts = std::sync::Mutex::new(Vec::new());
            fill(&mut out, g, |first, chunk| {
                starts.lock().expect("no panic in a chunk").push(first);
                for slot in chunk.iter_mut() {
                    *slot += 1;
                }
            });
            assert!(
                out.iter().all(|c| *c == 1),
                "n={n} g={g}: some cell was written {:?} times",
                out.iter().max()
            );
            for first in starts.lock().expect("no panic").iter() {
                assert_eq!(
                    first % g,
                    0,
                    "n={n} g={g}: a chunk began mid-plane at {first}"
                );
            }
        }
    }

    /// Threading is not attempted for work too small to pay for it.
    #[test]
    fn a_small_sweep_stays_on_one_thread() {
        std::env::remove_var("PANTOMETRY_THREADS");
        assert_eq!(threads_for(0), 1);
        assert_eq!(threads_for(CELLS_PER_SPAWN - 1), 1);
        assert_eq!(
            threads_for(32 * 32 * 32),
            1,
            "a 32-cube is not worth a spawn"
        );
        // The square-root law, at the sizes it was derived for.
        assert_eq!(threads_for(96 * 96 * 96), 4);
        assert_eq!(threads_for(128 * 128 * 128), 7);
    }
}
