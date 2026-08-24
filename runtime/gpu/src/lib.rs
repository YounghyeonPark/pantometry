//! An accelerator for pantometry's explicit stencils, with the CPU domain as the reference.
//!
//! # The one rule
//!
//! **`Solid3D` is the answer and this is a cache of it.** Where the two disagree, the CPU is
//! right. That is not modesty — it is the only arrangement under which the library's promises
//! survive a GPU at all, and it is what makes the comparison in `tests/` the point of the crate
//! rather than a nicety.
//!
//! # Why a GPU port is a different computation, not a faster one
//!
//! WGSL has no `f64`. Every number below is `f32`, against the domain's `f64`, so this is not the
//! same arithmetic run quicker — it is a lower-precision arithmetic. `tests/against_the_cpu.rs`
//! measures how far apart they land rather than asserting they agree, and the figure it reports
//! is the honest cost of the acceleration.
//!
//! That has a consequence worth stating plainly: `Simulation`'s conservation audit defaults to a
//! relative `1e-9`, and single precision cannot hold that. A scene using this needs
//! `conservation_tolerance_for(quantity::ENERGY, ..)` set to something `f32` can meet, and
//! choosing that number is choosing what the run is allowed to lose.
//!
//! # What is on the GPU and what is deliberately not
//!
//! The **stencil** is: each cell reads six neighbours and writes itself, in no particular order,
//! which is exactly the shape a GPU is for and exactly the shape that has no ordering problem.
//!
//! The **reductions are not.** A mean or a total summed with atomics depends on the order the
//! workgroups happened to finish, and floating-point addition is not associative — so the answer
//! would change between runs on one machine. [`GpuSolid::ledger`] reads the grid back and sums it
//! on the CPU in index order. That costs a transfer per audited step and buys back the thing the
//! whole workspace is built on.
//!
//! `Ensemble` solved the same problem on threads with fixed-size blocks, and the same discipline
//! would work here. It is not done yet because a readback is simpler and correct, and a faster
//! deterministic reduction is worth writing when somebody's grid is large enough to need it.
//!
//! # The buffer holds the **deviation**, and that is not a detail
//!
//! Cells carry `T − T₀`, not `T`. The first version stored absolute kelvin and diverged from the
//! reference by `1.4e-3` after two hundred steps — a thousand times what accumulation predicts —
//! and the cause was not accumulation at all.
//!
//! The update is `centre + F·(sum − 6·centre)`. On absolute temperatures near 293 K that sum is
//! about 1759, where `f32`'s resolution is `1.2e-4`. The difference being extracted from it is of
//! order `1e-3` K. **Subtracting two numbers that agree to five digits keeps less than one digit
//! of the answer**, every step, forever.
//!
//! On deviations the same numbers are near 1 K, where the resolution is `1.2e-7`, and the
//! subtraction keeps about four digits. The stencil is linear, so subtracting a constant commutes
//! with it exactly and the fix costs nothing — which is the useful part: single precision was not
//! the problem, spending it on an offset nobody needed was.

#![deny(missing_docs)]
#![forbid(unsafe_code)]

use pantometry_core::conserved::quantity;
use pantometry_core::units::{Length, Temperature, Time};
use pantometry_core::{Domain, Exchange, Ledger, Reading, Substance, Violation};
use pantometry_thermal::{Solid3D, STABLE_FOURIER_3D};

/// The heat channel, the same name the thermal crate publishes on.
pub const HEAT: &str = quantity::ENERGY;

/// A block of conducting material whose stencil runs on the GPU.
///
/// Mirrors `pantometry_thermal::Solid3D`: cubic cells, a seven-point stencil, insulated faces. What
/// differs is the arithmetic — `f32` here against `f64` there — and that difference is measured
/// rather than assumed.
pub struct GpuSolid {
    name: String,
    counts: (usize, usize, usize),
    dx: f64,
    alpha: f64,
    /// Heat capacity per cell, in J/K. **Per cell, not one number**: a block of two materials has
    /// two, and a void has none. `deposit` and `stored_heat` both divide by it.
    capacity: Vec<f64>,
    reference: f64,
    absorbed: f64,
    gpu: Context,
    /// The last state read back, as **deviations** from `reference`. See the crate docs.
    mirror: Vec<f32>,
    mirror_valid: bool,
}

/// The device, the pipeline and the two buffers the stencil ping-pongs between.
/// The part of a device context that does not depend on the block: the adapter, the queue, the
/// compiled shader and its layout.
///
/// **One per process, behind a `OnceLock`.** Every `GpuSolid::new` used to create its own
/// `wgpu::Instance`, request its own adapter and compile the shader again. Two costs, and the second
/// is the one that mattered: creating them **concurrently blocks**, so `cargo test --release` — the
/// command this crate's README gives — hung at the first row of the timing sweep while four other
/// tests held devices of their own. It ran in 32 s with `--test-threads=1` and not at all in ten
/// minutes without. Sharing is also simply right: the shader is identical every time.
struct Shared {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    /// What this actually ran on, from the adapter itself.
    ///
    /// A table of timings that does not name the device is prose nobody can check. `how_much_faster`
    /// prints it beside the numbers so the README's table has a machine attached to it.
    adapter: String,
}

/// One block's buffers, on the process's shared device.
struct Context {
    shared: std::sync::Arc<Shared>,
    cells: [wgpu::Buffer; 2],
    /// The resolved operator, uploaded once: face conductances, mobility, and the per-cell source
    /// rise. They do not change unless the block is rebuilt, so they cost one transfer and not one
    /// a step.
    coefficients: [wgpu::Buffer; 5],
    uniforms: wgpu::Buffer,
    readback: wgpu::Buffer,
    /// One per direction of the ping-pong, built once.
    ///
    /// **This used to be built per step.** Nothing in it changes except which of the two cell
    /// buffers is read, so there are exactly two and they are known at construction. Six thousand
    /// bind groups over a sweep is waste; removing it did not move the timings measurably, and the
    /// claim here is only that the work was never needed.
    binds: [wgpu::BindGroup; 2],
    /// Which of `cells` currently holds the state.
    front: usize,
    count: usize,
}

/// What went wrong before there was anything to compute.
#[derive(Clone, Debug)]
pub enum Unavailable {
    /// The block asks for something the device kernel has no pass for, and the reasons.
    ///
    /// Not a fall back to the CPU: a device is what a scene *states*, and a run that silently
    /// changed where it ran is a run whose answer nobody asked for.
    Unsupported(String),
    /// No GPU this process can reach.
    NoAdapter,
    /// One was found and would not give a device.
    NoDevice(String),
    /// The substance has no diffusivity, so there is no stencil to run.
    NotConducting,
}

impl std::fmt::Display for Unavailable {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unavailable::Unsupported(why) => {
                write!(f, "this block cannot run on a device: {why}")
            }
            Unavailable::NoAdapter => write!(f, "no GPU adapter"),
            Unavailable::NoDevice(e) => write!(f, "no GPU device: {e}"),
            Unavailable::NotConducting => write!(f, "the substance has no diffusivity"),
        }
    }
}

impl std::error::Error for Unavailable {}

impl GpuSolid {
    /// A block of `counts` cubic cells of side `dx`, all starting at `initial`.
    ///
    /// Fails rather than falling back to the CPU. A silent fallback is the worst outcome here:
    /// a caller who asked for a GPU and got a slow CPU run with no message has no way to tell,
    /// and the measurement they are about to take is of the wrong thing.
    pub fn new(
        name: impl Into<String>,
        substance: Substance,
        counts: (usize, usize, usize),
        dx: Length,
        initial: Temperature,
    ) -> Result<GpuSolid, Unavailable> {
        GpuSolid::mirroring(Solid3D::new(name, substance, counts, dx, initial))
    }

    /// Run **this** block's operator on the device.
    ///
    /// The CPU domain resolves the stencil — face conductances, mobility, sources — and this
    /// uploads what it resolved. Everything difficult about a real block is already in those
    /// arrays: two materials meeting are a face conductance that is neither one's, a void is a
    /// conductance and a mobility of zero, a coating is a thin row with its own `k`. Reading the
    /// substances and rebuilding the operator here would be a second implementation of `resolve`,
    /// and the first defect it had would be a physics difference nobody could see.
    ///
    /// # What it refuses, and why refusing is right
    ///
    /// A film, a gap exchange or a phase change is not a stencil, and the device has no pass for
    /// one. Asking for the device with any of those present is an **error** — see
    /// [`Unsupported`](Unavailable::Unsupported) — rather than a quiet fall back to the CPU. A run
    /// that silently changed where it ran is a run whose answer nobody asked for, and the device is
    /// something a scene *states*, not something a heuristic picks.
    pub fn mirroring(cpu: Solid3D) -> Result<GpuSolid, Unavailable> {
        let why = cpu.unsupported_on_a_device();
        if !why.is_empty() {
            return Err(Unavailable::Unsupported(why.join("; ")));
        }
        let counts = cpu.counts();
        let n = counts.0 * counts.1 * counts.2;
        let dx = cpu.spacing().to_si();
        let alpha = STABLE_FOURIER_3D * dx * dx / cpu.max_stable_dt(Time::from_si(0.0)).to_si();

        let c = cpu.coefficients();
        let (kx, ky, kz, mobility, source) = (c.kx, c.ky, c.kz, c.mobility, c.source);
        let down = |v: &[f64]| v.iter().map(|x| *x as f32).collect::<Vec<f32>>();
        // The source is watts; the kernel adds a rise, so the conversion happens once here rather
        // than every step on the device. `mobility` is `dx/C`, so `S·dt/C` is `S·dt·mobility/dx`.
        let rise: Vec<f32> = source
            .iter()
            .zip(mobility)
            .map(|(w, m)| (w * m / dx) as f32)
            .collect();
        let (fx, fy, fz, fm) = (down(kx), down(ky), down(kz), down(mobility));

        let reference = cpu.mean_temperature().to_si();
        let capacity: Vec<f64> = (0..n)
            .map(|c| {
                let (nx, ny, _) = counts;
                let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
                let _ = (i, j, k);
                cpu.cell_capacities()[c]
            })
            .collect();

        // Zero, because the buffer holds `T - T0` and `T` starts uniform at `T0`. A block whose
        // cells already differ is uploaded below.
        let mut gpu = Context::new(n, 0.0, [&fx, &fy, &fz, &fm, &rise])?;
        let mut mirror = vec![0.0f32; n];
        let (nx, ny, _) = counts;
        for (c, slot) in mirror.iter_mut().enumerate() {
            let (i, j, k) = (c % nx, (c / nx) % ny, c / (nx * ny));
            let t = cpu.temperature_at(i, j, k).to_si();
            // **A void arrives as zero, not as the absence the CPU reports.**
            //
            // `Solid3D::temperature_at` answers `NaN` for a void, deliberately: a void has no
            // temperature and a zero or an ambient there is a value somebody would plot. A device
            // buffer has no way to say that — and worse, `0.0 * NaN` is `NaN`, so a single absent
            // cell uploaded as one poisoned the whole grid within a few steps *even though its
            // face conductances are zero*. The first version of this did exactly that and every
            // cell came back `NaN`.
            //
            // Zero is safe because the cell cannot move: its mobility is zero, so the kernel
            // writes back what it read. What it holds means nothing, and nothing reads it — the
            // CPU is the reference and answers the question about voids.
            *slot = if t.is_finite() {
                (t - reference) as f32
            } else {
                0.0
            };
        }
        gpu.upload(&mirror);

        Ok(GpuSolid {
            name: pantometry_core::Domain::name(&cpu).to_string(),
            counts,
            dx,
            alpha,
            capacity,
            reference,
            absorbed: 0.0,
            gpu,
            mirror,
            mirror_valid: true,
        })
    }

    /// Cell counts along x, y and z.
    pub fn counts(&self) -> (usize, usize, usize) {
        self.counts
    }

    /// `α·dt/dx²`. Must stay at or under `STABLE_FOURIER_3D`.
    pub fn fourier_number(&self, dt: Time) -> f64 {
        self.alpha * dt.to_si() / (self.dx * self.dx)
    }

    /// The temperature of one cell, in kelvin.
    ///
    /// Reads the grid back if anything has changed since the last read, so a caller looping over
    /// every cell pays one transfer rather than one per cell.
    pub fn temperature_at(&mut self, i: usize, j: usize, k: usize) -> Temperature {
        self.sync();
        let (nx, ny, nz) = self.counts;
        let idx = i.min(nx - 1) + nx * (j.min(ny - 1) + ny * k.min(nz - 1));
        Temperature::from_si(self.reference + self.mirror[idx] as f64)
    }

    /// Set one cell, for an initial condition.
    pub fn set_temperature(&mut self, i: usize, j: usize, k: usize, t: Temperature) {
        self.sync();
        let (nx, ny, nz) = self.counts;
        if i >= nx || j >= ny || k >= nz {
            return;
        }
        self.mirror[i + nx * (j + ny * k)] = (t.to_si() - self.reference) as f32;
        self.gpu.upload(&self.mirror);
    }

    /// Put joules into one cell, as heat that arrived there.
    pub fn deposit(
        &mut self,
        i: usize,
        j: usize,
        k: usize,
        joules: pantometry_core::units::Energy,
    ) {
        let (nx, ny, nz) = self.counts;
        if i >= nx || j >= ny || k >= nz {
            return;
        }
        self.sync();
        let idx = i + nx * (j + ny * k);
        self.mirror[idx] += (joules.to_si() / self.capacity[idx]) as f32;
        self.absorbed += joules.to_si();
        self.gpu.upload(&self.mirror);
    }

    /// Every cell as a **deviation** from the starting temperature, in index order. One transfer.
    pub fn deviations(&mut self) -> &[f32] {
        self.sync();
        &self.mirror
    }

    /// Every cell as an absolute temperature in kelvin, in index order.
    pub fn cells(&mut self) -> Vec<f64> {
        self.sync();
        let reference = self.reference;
        self.mirror.iter().map(|v| reference + *v as f64).collect()
    }

    /// What this block is running on, as the adapter names itself.
    ///
    /// For a measurement to say where it was taken. A `191×` in this repository's prose had no
    /// machine, no build profile and no test behind it, and turned out to be nearer `48×`.
    pub fn device_name(&self) -> &str {
        &self.gpu.shared.adapter
    }

    /// Mean over every cell, summed **on the CPU in index order**.
    pub fn mean_temperature(&mut self) -> Temperature {
        self.sync();
        let n = self.mirror.len() as f64;
        Temperature::from_si(
            self.reference + self.mirror.iter().map(|v| *v as f64).sum::<f64>() / n,
        )
    }

    /// The hottest cell.
    pub fn peak_temperature(&mut self) -> Temperature {
        self.sync();
        Temperature::from_si(
            self.reference + self.mirror.iter().fold(f32::MIN, |m, v| m.max(*v)) as f64,
        )
    }

    fn sync(&mut self) {
        if !self.mirror_valid {
            self.gpu.download(&mut self.mirror);
            self.mirror_valid = true;
        }
    }

    /// Heat held, measured from the temperature it started at — the same reference `Solid3D` uses,
    /// and for the same reason: against absolute zero the interesting millijoule is lost in the
    /// last bits of a kilojoule.
    ///
    /// The buffer already holds deviations, which is exactly what a stored heat is measured from —
    /// the same fact that makes the `f32` arithmetic work at all.
    pub fn stored_heat(&mut self) -> pantometry_core::units::Energy {
        self.sync();
        pantometry_core::units::Energy::from_si(
            self.capacity
                .iter()
                .zip(&self.mirror)
                .map(|(c, d)| c * *d as f64)
                .sum::<f64>(),
        )
    }

    /// Heat taken from the bus over the run.
    pub fn absorbed_energy(&self) -> pantometry_core::units::Energy {
        pantometry_core::units::Energy::from_si(self.absorbed)
    }
}

impl Domain for GpuSolid {
    fn name(&self) -> &str {
        &self.name
    }

    /// `dx²/(6α)`, the same limit the CPU domain reports — the scheme is the same, only the
    /// precision differs.
    fn max_stable_dt(&self, _now: Time) -> Time {
        Time::from_si(STABLE_FOURIER_3D * self.dx * self.dx / self.alpha)
    }

    fn step(&mut self, _t: Time, dt: Time, bus: &mut Exchange) -> Result<(), Violation> {
        let f = self.fourier_number(dt);
        if f > STABLE_FOURIER_3D + 1e-12 {
            return Err(Violation {
                quantity: "Fourier number".to_string(),
                site: format!("{} (explicit 3D conduction, on the GPU)", self.name),
                before: STABLE_FOURIER_3D,
                after: f,
                scale: STABLE_FOURIER_3D,
                tolerance: 1e-12,
            });
        }

        // Placeless heat spreads evenly, exactly as the CPU domain does it — a block has six faces
        // and no distinguished cell, so choosing one would invent a location the bus never carried.
        let gained = bus.take_share(HEAT, dt);
        if gained != 0.0 {
            self.absorbed += gained;
            self.sync();
            // Spread to a uniform **rise**, in proportion to each cell's capacity rather than in
            // equal joules — equal joules would warm the low-capacity material more and so would
            // say where the heat landed, which the bus never carried. The CPU domain's rule.
            let total: f64 = self.capacity.iter().sum();
            let per_cell = if total > 0.0 {
                (gained / total) as f32
            } else {
                0.0
            };
            for cell in self.mirror.iter_mut() {
                *cell += per_cell;
            }
            self.gpu.upload(&self.mirror);
        }

        let (nx, ny, nz) = self.counts;
        // **Seconds, not the Fourier number.** The uniform-coefficient kernel folded `α dt/dx²`
        // into one scalar; the conductance form needs the step itself, because the coefficients
        // carry the rest. Passing the old scalar made the rise about eight hundred times too
        // large, and two hundred steps of that is a grid of `NaN`.
        self.gpu
            .dispatch(dt.to_si() as f32, [nx as u32, ny as u32, nz as u32]);
        let _ = f;
        self.mirror_valid = false;
        Ok(())
    }

    /// Heat gained since the start, summed on the CPU.
    fn ledger(&self) -> Ledger {
        // `ledger` takes `&self`, so the mirror cannot be refreshed here — which is exactly why
        // the audit reads a value that is one dispatch behind unless something has already synced.
        // `Simulation` calls `ledger` after every domain has stepped, and `readings` below syncs,
        // so in practice this is current. Stated because a stale ledger is an audit that passes.
        Ledger::new().with(
            quantity::ENERGY,
            self.capacity
                .iter()
                .zip(&self.mirror)
                .map(|(c, d)| c * *d as f64)
                .sum::<f64>(),
        )
    }

    fn books_balance(&self) -> bool {
        // Not claimed. The books would balance in `f64`; in `f32` the sum drifts from what the bus
        // moved by more than the per-domain check allows, and claiming otherwise would make a
        // correct domain fail. See the crate docs on what single precision costs.
        false
    }

    fn readings(&self) -> Vec<Reading> {
        let n = self.mirror.len().max(1) as f64;
        let mean = self.reference + self.mirror.iter().map(|v| *v as f64).sum::<f64>() / n;
        let peak = self.reference + self.mirror.iter().fold(f32::MIN, |m, v| m.max(*v)) as f64;
        let coldest = self.reference + self.mirror.iter().fold(f32::MAX, |m, v| m.min(*v)) as f64;
        vec![
            Reading::new(&self.name, "peak", peak - 273.15, "C"),
            Reading::new(&self.name, "mean", mean - 273.15, "C"),
            // **`coldest`, which was missing.** The CPU reports four scalars and this reported
            // three, so a scene run on the device produced a CSV with a column absent and nothing
            // to say a column was absent. Found by comparing the two runs of one scene rather than
            // by reading either.
            Reading::new(&self.name, "coldest", coldest - 273.15, "C"),
            Reading::new(&self.name, "absorbed", self.absorbed, "J"),
        ]
    }

    fn as_any(&self) -> Option<&dyn std::any::Any> {
        Some(self)
    }

    fn as_any_mut(&mut self) -> Option<&mut dyn std::any::Any> {
        Some(self)
    }
}

impl Drop for Context {
    /// **Give the buffers back to the driver now, not whenever.**
    ///
    /// Every `GpuSolid` used to own a whole `wgpu::Device`, so dropping one freed its allocations as
    /// a side effect of the device dying. Sharing the device removed that, and nine buffers a block
    /// at 8 MB each on a 128³ grid is 64 MB left to a lazy reclaim. An application that rebuilds a
    /// scene repeatedly — an editor, a batch of runs — is the case that cares.
    ///
    /// **This is hygiene, not a measured fix.** It was written to explain a 5.4× slowdown in the
    /// timing sweep and it did not explain it: the slowdown followed *position in the sweep* rather
    /// than grid size, dragged the CPU column along with it, and turned out to be the machine
    /// drifting under load. Kept because it is right, and recorded as unmeasured because the first
    /// version of this comment claimed a cause it had not established.
    fn drop(&mut self) {
        for b in self
            .cells
            .iter()
            .chain(self.coefficients.iter())
            .chain([&self.uniforms, &self.readback])
        {
            b.destroy();
        }
        // **And no poll here.** `destroy` is the part that is safe to do from a block's own drop.
        // Two attempts to make the device reclaim eagerly both wedged the suite past ten minutes:
        // `Maintain::Poll` after every `submit`, and `Maintain::Wait` in this drop — the latter
        // because `Wait` waits for the whole *device*, and three other tests were still submitting
        // to it. A shared device is shared: one block does not get to stop it.
    }
}

impl Shared {
    /// The process's device, made on first use and reused after.
    ///
    /// A failure is cached too: a machine with no adapter should not pay a probe per block, and the
    /// answer cannot change within a run.
    fn get() -> Result<std::sync::Arc<Shared>, Unavailable> {
        static SHARED: std::sync::OnceLock<Result<std::sync::Arc<Shared>, Unavailable>> =
            std::sync::OnceLock::new();
        SHARED.get_or_init(Shared::make).clone()
    }

    fn make() -> Result<std::sync::Arc<Shared>, Unavailable> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or(Unavailable::NoAdapter)?;
        let info = adapter.get_info();
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("pantometry stencil"),
                required_features: wgpu::Features::empty(),
                // **Seven storage buffers, and the downlevel default is four.** The uniform-
                // coefficient kernel needed two — a source and a destination — and the conductance
                // form needs five more: three face conductances, the mobility and the source.
                // Asked for explicitly so an adapter that cannot provide them fails here, with a
                // reason, rather than at `create_bind_group_layout` with a validation panic.
                required_limits: wgpu::Limits {
                    max_storage_buffers_per_shader_stage: 8,
                    ..wgpu::Limits::downlevel_defaults()
                },
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .map_err(|e| Unavailable::NoDevice(e.to_string()))?;

        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("stencil"),
            entries: &[
                entry(0, true),
                entry(1, false),
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                entry(3, true),
                entry(4, true),
                entry(5, true),
                entry(6, true),
                entry(7, true),
            ],
        });
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("stencil"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[&layout],
            push_constant_ranges: &[],
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("stencil"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: "sweep",
            compilation_options: Default::default(),
            cache: None,
        });

        Ok(std::sync::Arc::new(Shared {
            device,
            queue,
            pipeline,
            layout,
            adapter: format!("{} ({:?}, {:?})", info.name, info.device_type, info.backend),
        }))
    }
}

impl Context {
    fn new(count: usize, initial: f32, operator: [&[f32]; 5]) -> Result<Context, Unavailable> {
        let shared = Shared::get()?;
        let device = &shared.device;
        let queue = &shared.queue;

        let bytes = (count * 4) as u64;
        let make = |label: &str, usage: wgpu::BufferUsages| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(label),
                size: bytes.max(4),
                usage,
                mapped_at_creation: false,
            })
        };
        let storage = wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC;
        let cells = [make("cells a", storage), make("cells b", storage)];
        // One buffer each, sized to its own array: the face conductances are not cell-sized.
        let labels = ["kx", "ky", "kz", "mobility", "source"];
        let coefficients: [wgpu::Buffer; 5] = std::array::from_fn(|a| {
            let b = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some(labels[a]),
                size: ((operator[a].len() * 4) as u64).max(4),
                usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
            queue.write_buffer(
                &b,
                0,
                &operator[a]
                    .iter()
                    .flat_map(|f| f.to_le_bytes())
                    .collect::<Vec<u8>>(),
            );
            b
        });
        queue.write_buffer(
            &cells[0],
            0,
            &vec![initial; count]
                .iter()
                .flat_map(|f| f.to_le_bytes())
                .collect::<Vec<u8>>(),
        );

        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("stencil uniforms"),
            size: 32,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: bytes.max(4),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let binds = std::array::from_fn(|front| {
            let back = 1 - front;
            shared.device.create_bind_group(&wgpu::BindGroupDescriptor {
                label: Some(if front == 0 {
                    "sweep a to b"
                } else {
                    "sweep b to a"
                }),
                layout: &shared.layout,
                entries: &[
                    wgpu::BindGroupEntry {
                        binding: 0,
                        resource: cells[front].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 1,
                        resource: cells[back].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 2,
                        resource: uniforms.as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 3,
                        resource: coefficients[0].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 4,
                        resource: coefficients[1].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 5,
                        resource: coefficients[2].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 6,
                        resource: coefficients[3].as_entire_binding(),
                    },
                    wgpu::BindGroupEntry {
                        binding: 7,
                        resource: coefficients[4].as_entire_binding(),
                    },
                ],
            })
        });

        Ok(Context {
            shared,
            cells,
            coefficients,
            uniforms,
            readback,
            binds,
            front: 0,
            count,
        })
    }

    fn upload(&mut self, cells: &[f32]) {
        let bytes: Vec<u8> = cells.iter().flat_map(|f| f.to_le_bytes()).collect();
        self.shared
            .queue
            .write_buffer(&self.cells[self.front], 0, &bytes);
    }

    fn download(&mut self, into: &mut [f32]) {
        let mut encoder = self
            .shared
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.copy_buffer_to_buffer(
            &self.cells[self.front],
            0,
            &self.readback,
            0,
            (self.count * 4) as u64,
        );
        self.shared.queue.submit(Some(encoder.finish()));

        let slice = self.readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.shared.device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .expect("the map completes")
            .expect("a mapped buffer");
        {
            let mapped = slice.get_mapped_range();
            for (out, chunk) in into.iter_mut().zip(mapped.chunks_exact(4)) {
                *out = f32::from_le_bytes([chunk[0], chunk[1], chunk[2], chunk[3]]);
            }
        }
        self.readback.unmap();
    }

    /// One sweep: read the front buffer, write the back one, then swap.
    ///
    /// Two buffers rather than one, because a stencil that reads and writes the same array reads
    /// some neighbours already updated and some not — a Gauss-Seidel sweep pretending to be
    /// Jacobi, which is a different scheme with a different stability limit and no ordering
    /// anybody chose.
    fn dispatch(&mut self, dt: f32, counts: [u32; 3]) {
        let mut uniforms = Vec::with_capacity(32);
        uniforms.extend(counts[0].to_le_bytes());
        uniforms.extend(counts[1].to_le_bytes());
        uniforms.extend(counts[2].to_le_bytes());
        uniforms.extend((self.count as u32).to_le_bytes());
        uniforms.extend(dt.to_le_bytes());
        uniforms.extend([0u8; 12]);
        self.shared.queue.write_buffer(&self.uniforms, 0, &uniforms);

        let back = 1 - self.front;
        let mut encoder = self
            .shared
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("sweep"),
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.shared.pipeline);
            pass.set_bind_group(0, &self.binds[self.front], &[]);
            pass.dispatch_workgroups((self.count as u32).div_ceil(64), 1, 1);
        }
        self.shared.queue.submit(Some(encoder.finish()));
        self.front = back;
    }
}

fn entry(binding: u32, read_only: bool) -> wgpu::BindGroupLayoutEntry {
    wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    }
}

/// The seven-point stencil **in conductance form** — the same update `Solid3D` does, in `f32`.
///
/// `Cᵢ ΔTᵢ = dt Σ_f G_f (T_f − Tᵢ)`, with the per-face conductances and the per-cell mobility the
/// CPU domain resolved. It was `centre + F·(sum − 6·centre)`, which has no per-cell coefficient and
/// so is only the same operator when every cell is the same material — a limit the CPU's own
/// comment names, and the reason this port could run a homogeneous block and nothing else.
///
/// What the coefficients already carry, so the kernel does not have to: two materials meeting are a
/// face conductance that is neither one's; a **void** is a face conductance of zero and a mobility
/// of zero; a coating is a thin row with its own `k`. An outer face contributes no term at all,
/// which is the mirror boundary — a zero neighbour would be a wall at absolute zero and would drain
/// the block.
const SHADER: &str = r#"
struct Grid {
    nx: u32, ny: u32, nz: u32, total: u32,
    dt: f32, _pad0: f32, _pad1: f32, _pad2: f32,
};

@group(0) @binding(0) var<storage, read>       src: array<f32>;
@group(0) @binding(1) var<storage, read_write> dst: array<f32>;
@group(0) @binding(2) var<uniform>             g: Grid;
@group(0) @binding(3) var<storage, read>       kx: array<f32>;
@group(0) @binding(4) var<storage, read>       ky: array<f32>;
@group(0) @binding(5) var<storage, read>       kz: array<f32>;
@group(0) @binding(6) var<storage, read>       mobility: array<f32>;
@group(0) @binding(7) var<storage, read>       source: array<f32>;

@compute @workgroup_size(64)
fn sweep(@builtin(global_invocation_id) id: vec3<u32>) {
    let n = id.x;
    if (n >= g.total) { return; }

    let i = n % g.nx;
    let j = (n / g.nx) % g.ny;
    let k = n / (g.nx * g.ny);

    let t = src[n];
    let xr = (g.nx + 1u) * (j + g.ny * k);
    let yr = g.nx * (j + (g.ny + 1u) * k);
    var flux = 0.0;
    if (i > 0u)          { flux = flux + kx[i + xr]      * (src[n - 1u] - t); }
    if (i + 1u < g.nx)   { flux = flux + kx[i + 1u + xr] * (src[n + 1u] - t); }
    if (j > 0u)          { flux = flux + ky[i + yr]      * (src[n - g.nx] - t); }
    if (j + 1u < g.ny)   { flux = flux + ky[i + g.nx + yr] * (src[n + g.nx] - t); }
    if (k > 0u)          { flux = flux + kz[n]           * (src[n - g.nx * g.ny] - t); }
    if (k + 1u < g.nz)   { flux = flux + kz[n + g.nx * g.ny] * (src[n + g.nx * g.ny] - t); }

    // A source is a constant over the step and commutes with the stencil exactly, so adding it
    // after costs nothing in order — unlike a flux proportional to `T`, which would be Lie
    // splitting. `mobility` is `dx/C`, and a source is watts, so the rise is `S·dt/C`: the
    // division is folded into the coefficient the same way the CPU folds it.
    dst[n] = t + g.dt * (mobility[n] * flux + source[n]);
}
"#;

/// The [`Accelerator`](pantometry_world::Accelerator) a scene's `"device": "gpu"` is honoured by.
///
/// # Why this lives here and not in the library
///
/// `pantometry-world` is in the library's workspace: thirteen external crates, every one
/// licence-gated by `deny.toml`, all of them compiling to `wasm32` and to Rust 1.78. A wgpu stack is
/// eighty-six crates and none of those three things. So the scene format *carries* the request and
/// an application honours it, which is what this is for:
///
/// ```no_run
/// use pantometry_gpu::OnTheGpu;
/// use pantometry_world::{OnDisk, Scene, World};
///
/// # fn main() -> Result<(), String> {
/// let scene: Scene = serde_json::from_str(
///     r#"{ "title": "on the device", "duration_s": 0.02, "frames": 3,
///          "conservation_tolerance": 1e-4,
///          "domains": [ { "kind": "block", "name": "part", "cells": [64, 64, 64],
///                        "cell_mm": 1.0, "material": "aluminium", "initial_c": 20.0,
///                        "device": "gpu" } ] }"#,
/// )
/// .map_err(|e| e.to_string())?;
///
/// let mut world = World::build_with_accelerator(scene, &OnDisk, &OnTheGpu)?;
///
/// // `run` fails with a `Violation` — the conservation audit, which is why the scene above
/// // loosens it: `f32` cannot hold the default `1e-9`.
/// let frames = world.run().map_err(|v| v.to_string())?;
/// # let _ = frames;
/// # Ok(())
/// # }
/// ```
///
/// `no_run` rather than `ignore`: it was `ignore`, which is a snippet nobody compiles, and the point
/// of an example on a trait implementation is that the signature still fits.
///
/// # What it refuses
///
/// Everything `GpuSolid::mirroring` refuses, with the same reasons — a film, a gap exchange, a phase
/// change — plus a domain kind that has no device port at all. Never a quiet fall back to the CPU:
/// the scene said where to run and an answer from somewhere else is not the answer it asked for.
pub struct OnTheGpu;

impl pantometry_world::Accelerator for OnTheGpu {
    fn take(
        &self,
        spec: &pantometry_world::DomainSpec,
        device: pantometry_world::Device,
        cpu: Box<dyn Domain>,
    ) -> Result<Box<dyn Domain>, String> {
        if device != pantometry_world::Device::Gpu {
            return Ok(cpu);
        }
        // Downcast rather than rebuild from the spec: the block the library built has its materials,
        // voids, coatings and sources already resolved, and those resolved coefficients are the
        // whole reason the device can run a real block at all.
        let block = cpu
            .as_any()
            .and_then(|a| a.downcast_ref::<Solid3D>())
            .ok_or_else(|| {
                format!(
                    "{}: only a block runs on the gpu, and this domain is not one",
                    spec.name()
                )
            })?;
        match GpuSolid::mirroring(block.clone()) {
            Ok(gpu) => Ok(Box::new(gpu)),
            Err(why) => Err(format!("{}: {why}", spec.name())),
        }
    }
}
