//! Work done on several threads and summed on one, in one order: the same bits on any number of
//! threads.
//!
//! A floating-point sum is not associative, so a sum split across threads and combined depends on
//! how the work was divided ([`pantometry_core::sweep`] says the same of a reduction). What is
//! split here is never a sum. Each item — a span of atoms — is turned on some thread into what it
//! will add, each pair's terms, and the calling thread adds every item's in the items' order,
//! whatever order the threads finished them in. Each pair's terms are a function of the positions
//! and the parameters alone, so they are the same bits on any thread, and the additions are the
//! one-thread evaluation's, in its order. Or each item owns what it writes — a run of atoms'
//! forces, a slab of a grid, some of its columns ([`chunks`]) — and takes its additions there in
//! the one-thread order, and nothing is added across items at all.
//!
//! That assumes every thread rounds alike. A new thread starts in the default floating-point
//! environment, round to nearest with subnormals kept, and Rust offers no way to leave it; a host
//! that set flush-to-zero on the calling thread through other code would see terms with a
//! subnormal in them computed differently there.
//!
//! What it costs is the terms' round trip through memory, against the one-thread evaluation's
//! adding them where they are computed; and a spawn a thread a round, since a persistent pool
//! could not be written here without `unsafe` or a dependency ([`pantometry_core::sweep`]'s reason
//! too). A round's first thread runs about 80 µs after it starts and each further one about 35 µs
//! later ([`FIRST_THREAD_NS`]), so a part of an evaluation takes threads only where they shorten
//! it ([`threads_for_work`]). WebAssembly has no threads to spawn and takes one, the calling thread,
//! and so does a machine that refuses to spawn one.

use std::ops::Range;

/// How many threads this machine gives, read once: `PANTOMETRY_THREADS` when it is set to a
/// number, and the available parallelism otherwise; one on WebAssembly. Here it is a ceiling — an
/// evaluation too small to pay for a spawn still takes one — where
/// [`pantometry_core::sweep::threads_for`] takes it as the count. It is neither a clock nor a
/// random source, and it decides how long an evaluation takes, never what it computes.
pub(crate) fn available() -> usize {
    if cfg!(target_family = "wasm") {
        return 1;
    }
    static MACHINE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *MACHINE.get_or_init(|| {
        if let Ok(n) = std::env::var("PANTOMETRY_THREADS") {
            if let Ok(n) = n.parse::<usize>() {
                return n.max(1);
            }
        }
        std::thread::available_parallelism().map_or(1, |n| n.get())
    })
}

/// `0..n` cut into at most `pieces` runs of about equal length, in order.
pub(crate) fn ranges(n: usize, pieces: usize) -> Vec<Range<usize>> {
    let per = n.div_ceil(pieces.max(1)).max(1);
    (0..n.div_ceil(per))
        .map(|k| k * per..((k + 1) * per).min(n))
        .collect()
}

/// `slice` cut as [`ranges`] cuts its indices, each run with its range: disjoint, so that each
/// can go to a thread of its own.
pub(crate) fn chunks<T>(slice: &mut [T], pieces: usize) -> Vec<(Range<usize>, &mut [T])> {
    let mut rest = slice;
    let mut out = Vec::new();
    for r in ranges(rest.len(), pieces) {
        let (head, tail) = std::mem::take(&mut rest).split_at_mut(r.len());
        out.push((r, head));
        rest = tail;
    }
    out
}

/// How long a round of scoped spawns takes to have its first thread running, nanoseconds, and how
/// much later each further thread starts: measured as rounds shaped like [`in_order`]'s, each
/// after 300 µs of the calling thread's work alone as an evaluation's parts are, release,
/// `x86_64-pc-windows-gnu`, an i9-13900HX, medians of 200 rounds — the first thread running at 64
/// to 86 µs, and the last of two, four and eight at 117, 184 and 343 µs, 31 to 40 µs a thread.
/// Threads that do nothing, spawned and joined back to back, are 27 µs a thread, which a first
/// version took for the whole cost and then threaded parts that one thread finished sooner.
pub(crate) const FIRST_THREAD_NS: f64 = 80_000.0;

/// See [`FIRST_THREAD_NS`].
pub(crate) const SPAWN_NS: f64 = 35_000.0;

/// The threads a part of an evaluation that takes about `ns` nanoseconds on one thread takes:
/// `asked` when it was asked for, and otherwise the `T` that makes `T` spawns and `1/T` of the
/// work least, `⌊√(ns / SPAWN_NS)⌋` at most the machine's ([`available`]) — or one, when that
/// would take longer than one thread does, [`FIRST_THREAD_NS`] counted. One on WebAssembly.
pub(crate) fn threads_for_work(asked: Option<usize>, ns: f64) -> usize {
    if cfg!(target_family = "wasm") {
        return 1;
    }
    if let Some(t) = asked {
        return t.max(1);
    }
    by_work(ns, available())
}

/// [`threads_for_work`]'s default for a machine of `available` threads.
fn by_work(ns: f64, available: usize) -> usize {
    let t = ((ns / SPAWN_NS).sqrt() as usize).clamp(1, available.max(1));
    let threaded = FIRST_THREAD_NS + (t - 1) as f64 * SPAWN_NS + ns / t as f64;
    if t > 1 && threaded < ns {
        t
    } else {
        1
    }
}

/// The fewest and the most threads the parts of an evaluation spawned, which the tests read: an
/// evaluation that asked for several is to have spawned more than one in every part.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Spawned {
    /// The fewest any part spawned, `None` before the first.
    pub fewest: Option<usize>,
    /// The most any part spawned.
    pub most: usize,
    /// How many parts there were.
    pub parts: usize,
}

impl Spawned {
    /// Counts a part that spawned `threads`.
    pub(crate) fn add(&mut self, threads: usize) {
        self.fewest = Some(self.fewest.map_or(threads, |f| f.min(threads)));
        self.most = self.most.max(threads);
        self.parts += 1;
    }
}

#[cfg(all(test, not(target_family = "wasm")))]
thread_local! {
    /// Set by a test to stand for a machine that refuses every spawn.
    static REFUSED: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Whether a spawn is to be taken as refused: only ever in a test.
#[cfg(not(target_family = "wasm"))]
fn refused() -> bool {
    #[cfg(test)]
    {
        REFUSED.with(std::cell::Cell::get)
    }
    #[cfg(not(test))]
    {
        false
    }
}

/// Calls `produce(item, out)` for every item, on up to `threads` threads, and `consume(out)` on the
/// calling thread for each item's `out` **in the items' order**, and returns how many threads it
/// spawned. `out` is a buffer handed from item to item: `consume` leaves it as `produce` expects to
/// find it, empty. The buffers are taken from `spare` and given back to it, at most one more than
/// the threads, so a caller that keeps `spare` between calls allocates them about once.
///
/// On one thread, or one item, or WebAssembly, each item is produced and consumed in turn on the
/// calling thread, and none is spawned; and so they are on a machine that refuses the first spawn.
///
/// # Panics
///
/// If `produce` panics on any thread.
pub(crate) fn in_order<W, O, P, C>(
    threads: usize,
    items: Vec<W>,
    spare: &mut Vec<O>,
    produce: P,
    mut consume: C,
) -> usize
where
    W: Send,
    O: Default + Send,
    P: Fn(W, &mut O) + Sync,
    C: FnMut(&mut O),
{
    let count = items.len();
    if threads <= 1 || count <= 1 || cfg!(target_family = "wasm") {
        let mut out = spare.pop().unwrap_or_default();
        for item in items {
            produce(item, &mut out);
            consume(&mut out);
        }
        spare.push(out);
        return 0;
    }

    #[cfg(not(target_family = "wasm"))]
    {
        use std::sync::{mpsc, Mutex, MutexGuard, PoisonError};
        fn lock<T>(m: &Mutex<T>) -> MutexGuard<'_, T> {
            m.lock().unwrap_or_else(PoisonError::into_inner)
        }
        // The items are handed out in order, so the one the calling thread waits for has always
        // been taken before any later one: it is being produced, or it is done.
        let queue = Mutex::new(items.into_iter().enumerate());
        let pool: Mutex<Vec<O>> = Mutex::new(std::mem::take(spare));
        let (done, finished) = mpsc::channel::<(usize, O)>();
        let produce = &produce;
        let workers = std::thread::scope(|scope| {
            let mut workers = 0;
            for _ in 0..threads.min(count) {
                if refused() {
                    break;
                }
                let done = done.clone();
                let (queue, pool) = (&queue, &pool);
                let spawned = std::thread::Builder::new().spawn_scoped(scope, move || loop {
                    let next = lock(queue).next();
                    let Some((k, item)) = next else { break };
                    let mut out = lock(pool).pop().unwrap_or_default();
                    produce(item, &mut out);
                    if done.send((k, out)).is_err() {
                        break;
                    }
                });
                if spawned.is_err() {
                    break;
                }
                workers += 1;
            }
            // Only the threads hold a sender now, so a thread that panics ends the wait below
            // rather than leaving it waiting for an item nobody will finish.
            drop(done);
            if workers == 0 {
                // No thread was given: the calling thread takes every item, in order.
                let mut out = lock(&pool).pop().unwrap_or_default();
                loop {
                    let next = lock(&queue).next();
                    let Some((_, item)) = next else { break };
                    produce(item, &mut out);
                    consume(&mut out);
                }
                lock(&pool).push(out);
                return 0;
            }
            let mut waiting: Vec<Option<O>> = (0..count).map(|_| None).collect();
            let mut next = 0;
            while next < count {
                let (k, out) = finished
                    .recv()
                    .expect("a thread of an ordered evaluation panicked");
                waiting[k] = Some(out);
                while next < count {
                    let Some(mut out) = waiting[next].take() else {
                        break;
                    };
                    consume(&mut out);
                    lock(&pool).push(out);
                    next += 1;
                }
            }
            workers
        });
        let mut kept = pool.into_inner().unwrap_or_else(PoisonError::into_inner);
        kept.truncate(workers + 1);
        *spare = kept;
        workers
    }

    #[cfg(target_family = "wasm")]
    {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    /// Every third term negative: `% 3`, as `is_multiple_of` is Rust 1.87 and this is 1.78.
    fn term(k: usize) -> f64 {
        1.0 / (k as f64 + 0.37) * [-1.0, 1.0, 1.0][k % 3]
    }

    /// **The sum is the one-thread sum, to the bit, on any number of threads and any cut**: a sum
    /// of terms whose order changes its bits — summed backwards, or in five pieces each summed and
    /// then combined, they differ — on 1, 2, 3, 7 and 16 threads cut into 1 to 20 000 items. Every
    /// call with more than one thread and more than one item is asserted to have spawned its
    /// threads and to have produced on them rather than on the calling thread; WebAssembly spawns
    /// none and is asserted to.
    #[test]
    fn the_sum_is_the_one_threads_whatever_the_threads() {
        let n = 20_000;
        let mut reference = 0.0f64;
        for k in 0..n {
            reference += term(k);
        }
        let mut backwards = 0.0f64;
        for k in (0..n).rev() {
            backwards += term(k);
        }
        assert_ne!(reference.to_bits(), backwards.to_bits());
        let combined: f64 = ranges(n, 5)
            .into_iter()
            .map(|r| {
                let mut s = 0.0f64;
                for k in r {
                    s += term(k);
                }
                s
            })
            .fold(0.0, |a, b| a + b);
        assert_ne!(reference.to_bits(), combined.to_bits());

        let caller = std::thread::current().id();
        for threads in [1usize, 2, 3, 7, 16] {
            for pieces in [1usize, 2, 5, 64, 999, n] {
                let items = ranges(n, pieces);
                let count = items.len();
                let elsewhere = AtomicBool::new(false);
                let mut sum = 0.0f64;
                let mut seen = 0;
                let mut spare = Vec::new();
                let spawned = in_order(
                    threads,
                    items,
                    &mut spare,
                    |range: Range<usize>, out: &mut Vec<f64>| {
                        if std::thread::current().id() != caller {
                            elsewhere.store(true, Ordering::Relaxed);
                        }
                        // Uneven work, so that later items often finish first.
                        if range.start & 1 == 0 {
                            std::thread::yield_now();
                        }
                        out.extend(range.map(term));
                    },
                    |out: &mut Vec<f64>| {
                        for t in out.iter() {
                            sum += t;
                            seen += 1;
                        }
                        out.clear();
                    },
                );
                let label = format!("{threads} threads, {count} items");
                let split = threads > 1 && count > 1 && !cfg!(target_family = "wasm");
                let expected = if split { threads.min(count) } else { 0 };
                assert_eq!(spawned, expected, "{label}");
                assert_eq!(elsewhere.into_inner(), split, "{label}");
                // The buffers kept for the next call: at most one more than the threads.
                assert!(
                    (1..=spawned + 1).contains(&spare.len()),
                    "{label}: {}",
                    spare.len()
                );
                assert_eq!(seen, n, "{label}");
                assert_eq!(
                    sum.to_bits(),
                    reference.to_bits(),
                    "{label}: {sum} against {reference}"
                );
            }
        }
    }

    /// **An item finished before the one ahead of it waits its turn**, by construction rather
    /// than by the scheduler: on two threads, item 0 is not finished until item 2 has been
    /// started, which the other thread does only after it has sent item 1. So item 1 always
    /// arrives first, and the sum is still item 0's terms, then 1's, then 2's — which every other
    /// order of the three would change.
    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn an_item_finished_first_waits_its_turn() {
        let items: Vec<Range<usize>> = vec![0..700, 700..1400, 1400..2100];
        let sum_in = |order: [usize; 3]| {
            let mut s = 0.0f64;
            for i in order {
                for k in items[i].clone() {
                    s += term(k);
                }
            }
            s.to_bits()
        };
        let reference = sum_in([0, 1, 2]);
        for other in [[1, 0, 2], [1, 2, 0], [0, 2, 1], [2, 0, 1], [2, 1, 0]] {
            assert_ne!(sum_in(other), reference, "{other:?} sums alike");
        }
        let started = AtomicBool::new(false);
        let mut sum = 0.0f64;
        let spawned = in_order(
            2,
            items.clone(),
            &mut Vec::new(),
            |range: Range<usize>, out: &mut Vec<f64>| {
                if range.start == 1400 {
                    started.store(true, Ordering::Release);
                }
                if range.start == 0 {
                    let since = std::time::Instant::now();
                    while !started.load(Ordering::Acquire) {
                        assert!(since.elapsed().as_secs() < 60, "item 2 never started");
                        std::thread::yield_now();
                    }
                }
                out.extend(range.map(term));
            },
            |out: &mut Vec<f64>| {
                for t in out.iter() {
                    sum += t;
                }
                out.clear();
            },
        );
        assert_eq!(spawned, 2);
        assert_eq!(sum.to_bits(), reference);
    }

    /// A machine that refuses to spawn a thread gets the calling thread's evaluation: every item,
    /// in order, and nothing spawned.
    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn a_machine_that_gives_no_thread_takes_the_calling_one() {
        let n = 5_000;
        let mut reference = 0.0f64;
        for k in 0..n {
            reference += term(k);
        }
        REFUSED.with(|r| r.set(true));
        let caller = std::thread::current().id();
        let elsewhere = AtomicBool::new(false);
        let mut sum = 0.0f64;
        let mut spare = Vec::new();
        let spawned = in_order(
            8,
            ranges(n, 40),
            &mut spare,
            |range: Range<usize>, out: &mut Vec<f64>| {
                if std::thread::current().id() != caller {
                    elsewhere.store(true, Ordering::Relaxed);
                }
                out.extend(range.map(term));
            },
            |out: &mut Vec<f64>| {
                for t in out.iter() {
                    sum += t;
                }
                out.clear();
            },
        );
        REFUSED.with(|r| r.set(false));
        assert_eq!(spawned, 0);
        assert!(!elsewhere.into_inner());
        assert_eq!(sum.to_bits(), reference.to_bits());
        assert_eq!(spare.len(), 1);
    }

    /// **A part takes threads only where they pay**: as many as asked when asked, one for nothing to
    /// do, and by default `⌊√(ns / SPAWN_NS)⌋` of the machine's where that beats one thread with the
    /// first thread's start counted, and one where it does not. Two threads pay from 230 µs exactly,
    /// `80 + 35 + 230/2`: the octanol box's transform along the first axis, 229 µs by the model
    /// (`1.4 ns × 32³ × 5`, measured 0.23 ms), stays on one, and its splines, 850 µs
    /// (`167 ns × 5 089`), take four. The rule is checked exactly for a machine of 32 and of three,
    /// whatever this one is, and on this one as [`threads_for_work`] reads it. One on WebAssembly,
    /// asked or not.
    #[test]
    fn a_part_takes_threads_only_where_they_pay() {
        let wasm = cfg!(target_family = "wasm");
        for t in [1usize, 2, 7, 64] {
            let want = if wasm { 1 } else { t };
            assert_eq!(threads_for_work(Some(t), 0.0), want);
            assert_eq!(threads_for_work(Some(t), 1e12), want);
        }
        assert_eq!(threads_for_work(Some(0), 1e12), 1);
        for (ns, on_32, on_3) in [
            (0.0, 1, 1),
            (1.4 * 32_768.0 * 5.0, 1, 1),
            (230_000.0, 1, 1),
            (230_001.0, 2, 2),
            (167.0 * 5_089.0, 4, 3),
            (1e12, 32, 3),
        ] {
            assert_eq!(by_work(ns, 32), on_32, "{ns} ns on 32");
            assert_eq!(by_work(ns, 3), on_3, "{ns} ns on 3");
            assert_eq!(by_work(ns, 1), 1, "{ns} ns on 1");
            let here = if wasm { 1 } else { by_work(ns, available()) };
            assert_eq!(threads_for_work(None, ns), here, "{ns} ns here");
        }
    }

    /// `chunks` cuts a slice where `ranges` cuts its indices, each piece its own range's elements.
    #[test]
    fn the_chunks_are_the_runs() {
        for n in [0usize, 1, 7, 100, 5089] {
            for pieces in [1usize, 3, 8, 64, 10_000] {
                let mut v: Vec<usize> = (0..n).collect();
                let runs = ranges(n, pieces);
                let cut = chunks(&mut v, pieces);
                assert_eq!(cut.len(), runs.len(), "{n} in {pieces}");
                for ((range, chunk), run) in cut.into_iter().zip(&runs) {
                    assert_eq!(&range, run);
                    assert_eq!(chunk.to_vec(), range.collect::<Vec<_>>());
                }
            }
        }
    }

    /// The runs cover `0..n` once each, in order, and are no more than asked for.
    #[test]
    fn the_runs_cover_every_index_once_in_order() {
        for n in [0usize, 1, 2, 7, 100, 5089] {
            for pieces in [0usize, 1, 2, 3, 8, 64, 10_000] {
                let r = ranges(n, pieces);
                assert!(
                    r.len() <= pieces.max(1),
                    "{n} cut into {} of {pieces}",
                    r.len()
                );
                let mut next = 0;
                for run in &r {
                    assert_eq!(run.start, next, "{n} in {pieces}: a gap or an overlap");
                    assert!(run.end > run.start, "{n} in {pieces}: an empty run");
                    next = run.end;
                }
                assert_eq!(next, n, "{n} in {pieces}: not every index");
            }
        }
    }

    /// **A thread that panics ends the evaluation with a panic, not a wait** for an item that will
    /// never come: run on a thread of its own and waited for at most a minute, so that the wait
    /// fails this test rather than holding the run, and the panic is the calling thread's, from
    /// the threads' path — a panic on the calling thread would carry the item's own message.
    #[test]
    #[cfg(not(target_family = "wasm"))]
    fn a_panic_on_a_thread_is_a_panic_here() {
        let (tell, told) = std::sync::mpsc::channel();
        std::thread::spawn(move || {
            let caught = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                in_order(
                    4,
                    ranges(100, 10),
                    &mut Vec::new(),
                    |range: Range<usize>, out: &mut Vec<usize>| {
                        assert!(range.start != 30, "the item that fails");
                        out.extend(range);
                    },
                    |out: &mut Vec<usize>| out.clear(),
                )
            }));
            let message = caught.err().map(|e| {
                e.downcast_ref::<String>()
                    .cloned()
                    .or_else(|| e.downcast_ref::<&str>().map(|s| s.to_string()))
                    .unwrap_or_default()
            });
            let _ = tell.send(message);
        });
        let message = told
            .recv_timeout(std::time::Duration::from_secs(60))
            .expect("the evaluation waited rather than panicking")
            .expect("a thread's panic is a panic here");
        assert!(
            message.contains("a thread of an ordered evaluation panicked"),
            "{message}"
        );
    }
}
