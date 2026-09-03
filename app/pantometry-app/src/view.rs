//! A native window for a pantometry run.
//!
//! ```text
//! cargo run --release -- run.json
//! ```
//!
//! Drag to rotate, scroll to zoom, space to play, left and right to scrub.
//!
//! ```text
//! cargo run --release -- run.json --snapshot out.ppm
//! ```
//!
//! renders one frame **without opening a window** and writes it out. That mode is not a
//! convenience: it is the only way anything can check that this file draws. A window nobody can
//! photograph proves it did not panic, which is a much weaker claim than it looks — every wrong
//! projection, every empty vertex buffer and every silently-failed pipeline also does not panic.
//!
//! # What is here and what is deliberately not
//!
//! This is the shell. Everything that could be got wrong twice — the colour scale across a run,
//! the framing, the projection — is in `viewer-core`, which has no GPU dependency and is tested
//! against real run files. What is left here is a surface, a line pipeline and an event loop.
//!
//! That split is not tidiness. A renderer is the one place where a wrong answer looks like a
//! picture, so the arithmetic lives where a test can reach it and this file only draws what it is
//! handed.
//!
//! # It does not depend on `pantometry`
//!
//! Deliberately. It reads the JSON a run wrote and nothing else, so "the wire format carries
//! enough to draw a run" is demonstrated rather than asserted. If this needed to link the library
//! for something the file did not have, the format would be the thing to fix.

use std::sync::Arc;

use image::ImageEncoder;
use viewer_core::{segments, Camera, Framing, Run};
use winit::application::ApplicationHandler;
use winit::event::{ElementState, MouseScrollDelta, WindowEvent};
use winit::event_loop::{ActiveEventLoop, ControlFlow, EventLoop};
use winit::keyboard::{KeyCode, PhysicalKey};
use winit::window::{Window, WindowId};

/// One end of a line, as the shader wants it.
#[repr(C)]
#[derive(Clone, Copy)]
struct Vertex {
    position: [f32; 2],
    colour: [f32; 3],
}

/// A stand-in for `bytemuck`, which is one dependency for one trait.
mod bytemuck_lite {
    /// Marker for a type that can be reinterpreted as bytes.
    ///
    /// # Safety
    ///
    /// Implementors must be `#[repr(C)]`, contain no padding and no references. Both types here
    /// are arrays of `f32`, which satisfies all three.
    pub unsafe trait Pod: Copy {}

    /// The bytes of a slice.
    pub fn cast_slice<T: Pod>(v: &[T]) -> &[u8] {
        // SAFETY: `T: Pod` promises a plain-data layout, and the length is derived from the
        // slice's own, so the range is exactly the memory the slice owns.
        unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
    }
}

// SAFETY: `#[repr(C)]`, five `f32` with no padding and no references.
unsafe impl bytemuck_lite::Pod for Vertex {}

/// `pantometry view <run.json> [--snapshot out.ppm]`.
pub fn run(args: &[String]) -> i32 {
    let path = match args.first() {
        Some(p) => p.clone(),
        None => {
            eprintln!("usage: pantometry view <run.json> [--snapshot out.ppm]");
            eprintln!("  produced by `pantometry run <scene> out.json`, or by any run that calls");
            eprintln!("  pantometry_view::to_json");
            return 2;
        }
    };
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("cannot read {path}: {e}");
            std::process::exit(1);
        }
    };
    let run = match Run::from_json(&text) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    };

    // Which panel to show: the first one there is.
    //
    // **It used to be the first one with `Paths`**, and it refused everything else — "fields and
    // point clouds are different pipelines and are not built yet". They are not different
    // pipelines: a body and a field sample are points, and a point is two short segments in
    // screen space, which is the pipeline that was already here. What the old rule cost is that
    // `PanelData::paths` is built by two examples and by tests and by no `Domain` at all, so
    // **none of the thirty shipped scenes could be drawn by this shell** — the refusal was the
    // only thing it ever said about a scene.
    let Some(panel) = run.panels().into_iter().next() else {
        eprintln!("{}: this run has no panels at all", run.title);
        std::process::exit(1);
    };

    println!(
        "{} — {} frames, showing {panel}",
        run.title,
        run.frames.len()
    );
    println!("  drag to rotate, scroll to zoom, space to play, left/right to scrub");

    // Headless: render one frame to a texture, read it back, write a PPM. No window, no display,
    // so a machine with neither can still check that the renderer puts something on the canvas.
    // From the dispatcher, past the run file. `std::env::args().skip(2)` was right when this was
    // its own binary and is off by one as a subcommand — it would have read the run file's path as
    // the flag.
    // **Found by scanning, like its siblings.** `--snapshot` had to be the argument straight after
    // the run file while `--frame` and `--thumbnail` are looked for anywhere, and the asymmetry
    // was silent both ways: `view run.json --thumbnail --snapshot out.png` fell past this branch
    // and opened a window, and `view run.json --snapshot --frame 5` took `--frame` as the output
    // path and wrote a 1100x720 PPM into a file called `--frame`.
    let rest: Vec<&str> = args.iter().skip(1).map(String::as_str).collect();
    if let Some(i) = rest.iter().position(|a| *a == "--snapshot") {
        let out = match rest.get(i + 1) {
            Some(p) if !p.starts_with("--") => (*p).to_string(),
            _ => {
                eprintln!(
                    "--snapshot takes a file to write; got {:?}",
                    rest.get(i + 1)
                );
                return 2;
            }
        };
        // Which frame, because the first one is often the least interesting: a run that fills up
        // over its length — a shot, a spreading spot — has nothing in it at `t = 0`, and a
        // snapshot of that is a picture of an empty box that still counts as "the renderer works".
        let at: usize = rest
            .iter()
            .position(|a| *a == "--frame")
            .and_then(|i| rest.get(i + 1))
            .and_then(|n| n.parse().ok())
            .unwrap_or(0)
            .min(run.frames.len().saturating_sub(1));
        // **A tile for the New-project chooser**: cropped to what was drawn and shrunk. The
        // cropping is the point — the camera fits the panel's bounds, and a one-dimensional bar
        // is a hairline in a frame that is otherwise background.
        let tile = rest.contains(&"--thumbnail");
        let mut app = App::new(run, panel);
        app.frame = at;
        println!("  snapshot of frame {at} of {}", app.run.frames.len());
        // Rendered three times the tile and averaged down, because a marker one pixel wide
        // vanishes under nearest-neighbour and a tile is not where to discover that.
        let (rw, rh) = if tile {
            (THUMB.0 * 3, THUMB.1 * 3)
        } else {
            (1100, 720)
        };
        match app.snapshot(rw, rh) {
            Ok(pixels) => {
                let (pw, ph, pixels) = if tile {
                    let (cw, ch, cropped) = crop_to_content(rw, rh, &pixels);
                    (THUMB.0, THUMB.1, shrink(cw, ch, &cropped, THUMB.0, THUMB.1))
                } else {
                    (rw, rh, pixels)
                };
                write_image(&out, pw, ph, &pixels);
                // **Compared against the corner pixel, not against a constant.** The target is
                // sRGB, so the clear colour is stored far brighter than the linear number the
                // pass was given — 56,66,77 rather than 10,14,19. A fixed threshold called every
                // pixel in the image a line and reported 100%, which is exactly the kind of
                // "measurement" a renderer check exists to avoid.
                let background = [pixels[0], pixels[1], pixels[2]];
                // `as_chunks` rather than `chunks_exact`: each pixel is a `[u8; 4]` by type, which
                // is what it is. This workspace has no MSRV to keep it off — see the note beside
                // the same call in `pantometry-gpu`.
                let (rgba, _) = pixels.as_chunks::<4>();
                let lit = rgba.iter().filter(|p| p[..3] != background).count();
                let mut shades: Vec<[u8; 3]> = Vec::new();
                for p in rgba {
                    let c = [p[0], p[1], p[2]];
                    if c != background && !shades.contains(&c) {
                        shades.push(c);
                    }
                }
                // **The image's own size, not the default one.** `lit` is counted over
                // `pixels`, which under `--thumbnail` is 240x156 and not 1100x720 — so a tile
                // lighting 173 of its 37 440 was reported as "173 of 792000 (0.02%)", a fifth of
                // a percent read as a fiftieth. The count was right and the denominator was a
                // constant from before there was a second size.
                let all = pw as u64 * ph as u64;
                println!(
                    "  wrote {out} — {lit} of {all} pixels carry a line ({:.2}%), in {} shades",
                    100.0 * lit as f64 / all as f64,
                    shades.len()
                );
                if lit == 0 {
                    eprintln!("  nothing was drawn");
                    std::process::exit(1);
                }
            }
            Err(e) => {
                eprintln!("  no GPU available for a snapshot: {e}");
                std::process::exit(3);
            }
        }
        return 0;
    }

    let event_loop = EventLoop::new().expect("an event loop");
    event_loop.set_control_flow(ControlFlow::Poll);
    let mut app = App::new(run, panel);
    event_loop.run_app(&mut app).expect("the window runs");
    0
}

/// Everything the window holds between frames.
struct App {
    run: Run,
    panel: String,
    framing: Framing,
    camera: Camera,
    /// The run-wide range the shading is measured against.
    span: (f64, f64),
    frame: usize,
    playing: bool,
    dragging: Option<(f64, f64)>,
    gpu: Option<Gpu>,
}

struct Gpu {
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    pipeline: wgpu::RenderPipeline,
}

impl App {
    fn new(run: Run, panel: String) -> App {
        let framing = Framing::of(
            run.framing_of(&panel)
                .unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
        );
        // Once, from the whole run. Re-fitting it per frame is what `Run::scale_of` exists to
        // stop, and for a while nothing called it.
        let span = run.scale_of(&panel).unwrap_or((0.0, 1.0));
        // Framed to what is actually there. The window can still be zoomed; `--snapshot` cannot,
        // and a fixed distance is a distance chosen for a cube.
        let mut camera = Camera::default();
        camera.fit(
            run.framing_of(&panel)
                .unwrap_or([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0]),
            &framing,
            16.0 / 9.0,
            0.85,
        );
        App {
            run,
            panel,
            framing,
            span,
            camera,
            frame: 0,
            playing: true,
            dragging: None,
            gpu: None,
        }
    }

    /// The line vertices for the current frame.
    fn vertices(&self, aspect: f64) -> Vec<Vertex> {
        let Some(panel) = self
            .run
            .frames
            .get(self.frame)
            .and_then(|f| f.panels.iter().find(|p| p.name() == self.panel))
        else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for s in segments(panel, &self.camera, &self.framing, aspect, self.span) {
            let colour = ramp(s.shade);
            out.push(Vertex {
                position: [s.from.x as f32, s.from.y as f32],
                colour,
            });
            out.push(Vertex {
                position: [s.to.x as f32, s.to.y as f32],
                colour,
            });
        }
        out
    }
}

impl ApplicationHandler for App {
    fn resumed(&mut self, event_loop: &ActiveEventLoop) {
        if self.gpu.is_some() {
            return;
        }
        let attrs = Window::default_attributes()
            .with_title(format!("pantometry — {}", self.run.title))
            .with_inner_size(winit::dpi::LogicalSize::new(1100.0, 720.0));
        let window = Arc::new(event_loop.create_window(attrs).expect("a window"));

        let instance = wgpu::Instance::default();
        let surface = instance
            .create_surface(window.clone())
            .expect("a surface on that window");
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: Some(&surface),
            force_fallback_adapter: false,
        }))
        .expect("a GPU adapter");
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("pantometry viewer"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .expect("a device");

        let size = window.inner_size();
        let caps = surface.get_capabilities(&adapter);
        let config = wgpu::SurfaceConfiguration {
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            format: caps.formats[0],
            width: size.width.max(1),
            height: size.height.max(1),
            present_mode: wgpu::PresentMode::Fifo,
            alpha_mode: caps.alpha_modes[0],
            view_formats: vec![],
            desired_maximum_frame_latency: 2,
        };
        surface.configure(&device, &config);

        let pipeline = line_pipeline(&device, config.format);

        self.gpu = Some(Gpu {
            window,
            surface,
            device,
            queue,
            config,
            pipeline,
        });
    }

    fn window_event(&mut self, event_loop: &ActiveEventLoop, _id: WindowId, event: WindowEvent) {
        let Some(gpu) = self.gpu.as_mut() else {
            return;
        };
        match event {
            WindowEvent::CloseRequested => event_loop.exit(),
            WindowEvent::Resized(size) => {
                gpu.config.width = size.width.max(1);
                gpu.config.height = size.height.max(1);
                gpu.surface.configure(&gpu.device, &gpu.config);
            }
            WindowEvent::MouseInput { state, .. } => {
                self.dragging = matches!(state, ElementState::Pressed).then_some((0.0, 0.0));
            }
            WindowEvent::CursorMoved { position, .. } => {
                if let Some((lx, ly)) = self.dragging {
                    let (x, y) = (position.x, position.y);
                    if lx != 0.0 || ly != 0.0 {
                        self.camera.turn((x - lx) * 0.008, (y - ly) * 0.006);
                    }
                    self.dragging = Some((x, y));
                }
            }
            WindowEvent::MouseWheel { delta, .. } => {
                let step = match delta {
                    MouseScrollDelta::LineDelta(_, y) => y as f64,
                    MouseScrollDelta::PixelDelta(p) => p.y / 60.0,
                };
                self.camera.zoom(1.0 - step * 0.1);
            }
            WindowEvent::KeyboardInput { event, .. } if event.state == ElementState::Pressed => {
                let n = self.run.frames.len().max(1);
                match event.physical_key {
                    PhysicalKey::Code(KeyCode::Space) => self.playing = !self.playing,
                    PhysicalKey::Code(KeyCode::ArrowRight) => {
                        self.playing = false;
                        self.frame = (self.frame + 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::ArrowLeft) => {
                        self.playing = false;
                        self.frame = (self.frame + n - 1) % n;
                    }
                    PhysicalKey::Code(KeyCode::Escape) => event_loop.exit(),
                    _ => {}
                }
            }
            WindowEvent::RedrawRequested => self.draw(),
            _ => {}
        }
    }

    fn about_to_wait(&mut self, _event_loop: &ActiveEventLoop) {
        if let Some(gpu) = &self.gpu {
            gpu.window.request_redraw();
        }
    }
}

impl App {
    fn draw(&mut self) {
        let Some(gpu) = self.gpu.as_ref() else {
            return;
        };
        let aspect = gpu.config.width as f64 / gpu.config.height.max(1) as f64;
        let verts = self.vertices(aspect);

        let Ok(surface_texture) = gpu.surface.get_current_texture() else {
            return;
        };
        let view = surface_texture
            .texture
            .create_view(&wgpu::TextureViewDescriptor::default());
        let buffer = if verts.is_empty() {
            None
        } else {
            Some(
                gpu.device
                    .create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)),
            )
        };

        let mut encoder = gpu
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("lines"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            if let Some(b) = &buffer {
                pass.set_pipeline(&gpu.pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        gpu.queue.submit(Some(encoder.finish()));
        surface_texture.present();

        if self.playing && !self.run.frames.is_empty() {
            self.frame = (self.frame + 1) % self.run.frames.len();
        }
    }
}

impl App {
    /// Render one frame with no window at all, and hand back the pixels as RGBA.
    ///
    /// The same camera, the same `viewer-core` segments and the same pipeline the window uses —
    /// only the target differs. A snapshot path that built its own vertices would be checking
    /// itself.
    fn snapshot(&self, width: u32, height: u32) -> Result<Vec<u8>, String> {
        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .ok_or("no adapter")?;
        let (device, queue) = pollster::block_on(adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: Some("snapshot"),
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
                memory_hints: wgpu::MemoryHints::default(),
            },
            None,
        ))
        .map_err(|e| e.to_string())?;

        let format = wgpu::TextureFormat::Rgba8UnormSrgb;
        let target = device.create_texture(&wgpu::TextureDescriptor {
            label: Some("snapshot"),
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
            view_formats: &[],
        });
        let view = target.create_view(&wgpu::TextureViewDescriptor::default());
        let pipeline = line_pipeline(&device, format);

        let verts = self.vertices(width as f64 / height as f64);
        let buffer = (!verts.is_empty())
            .then(|| device.create_buffer_init_lite(bytemuck_lite::cast_slice(&verts)));

        // Copies out of a texture want rows padded to 256 bytes, so the readback is wider than the
        // image and the padding is dropped below rather than left in the file as a smear.
        let unpadded = width * 4;
        let padded = unpadded.div_ceil(256) * 256;
        let readback = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback"),
            size: (padded * height) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });

        let mut encoder =
            device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("snapshot"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: None,
                timestamp_writes: None,
                occlusion_query_set: None,
            });
            if let Some(b) = &buffer {
                pass.set_pipeline(&pipeline);
                pass.set_vertex_buffer(0, b.slice(..));
                pass.draw(0..verts.len() as u32, 0..1);
            }
        }
        encoder.copy_texture_to_buffer(
            wgpu::ImageCopyTexture {
                texture: &target,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::ImageCopyBuffer {
                buffer: &readback,
                layout: wgpu::ImageDataLayout {
                    offset: 0,
                    bytes_per_row: Some(padded),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        queue.submit(Some(encoder.finish()));

        let slice = readback.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        device.poll(wgpu::Maintain::Wait);
        rx.recv()
            .map_err(|e| e.to_string())?
            .map_err(|e| e.to_string())?;

        let mapped = slice.get_mapped_range();
        let mut out = Vec::with_capacity((unpadded * height) as usize);
        for row in 0..height {
            let from = (row * padded) as usize;
            out.extend_from_slice(&mapped[from..from + unpadded as usize]);
        }
        drop(mapped);
        readback.unmap();
        Ok(out)
    }
}

/// A plain PPM, because a picture nothing can open is not evidence and an encoder is a dependency.
/// A thumbnail's size. Wide enough to tell a room's standing wave from a block's hot core in a
/// grid of thirty, small enough that thirty of them are under two hundred kilobytes.
pub const THUMB: (u32, u32) = (240, 156);

/// Trim the background around what was drawn, then pad back to `THUMB`'s aspect.
///
/// **The camera fits the panel's bounds, and a bar is a line across an otherwise empty frame.**
/// At thumbnail size that is a picture of nothing: measured on `05-beam-on-bar`, where the drawn
/// content was a diagonal hairline in a frame that was 95% background. Cropping to the content is
/// what makes a grid of these worth looking at.
///
/// The background is the corner pixel, which is what the snapshot's own "not background" count
/// already assumes: the shell follows no theme, but the clear colour is one place and this reads
/// it rather than naming it twice.
///
/// # The magnification is clamped, and the number is chosen rather than derived
///
/// A scene with one body has a content box one marker wide, and cropping to that fills the tile
/// with a single cross — a picture of the *marker* rather than of the run. Measured on
/// `07-bouncing-ball`, which is exactly one body. **Three** is a rendering taste; it is written
/// down as one so the next reader does not go looking for where it came from.
fn crop_to_content(w: u32, h: u32, rgba: &[u8]) -> (u32, u32, Vec<u8>) {
    const MAGNIFY: f64 = 3.0;
    const MARGIN: f64 = 0.06;
    let at = |x: u32, y: u32| ((y * w + x) * 4) as usize;
    let bg = [rgba[0], rgba[1], rgba[2]];
    let (mut x0, mut y0, mut x1, mut y1) = (w as i64, h as i64, -1i64, -1i64);
    for y in 0..h {
        for x in 0..w {
            let k = at(x, y);
            if (0..3).any(|c| rgba[k + c].abs_diff(bg[c]) > 6) {
                x0 = x0.min(x as i64);
                y0 = y0.min(y as i64);
                x1 = x1.max(x as i64);
                y1 = y1.max(y as i64);
            }
        }
    }
    if x1 < x0 || y1 < y0 {
        return (w, h, rgba.to_vec());
    }

    let (mut cw, mut ch) = ((x1 - x0 + 1) as f64, (y1 - y0 + 1) as f64);
    let (cx, cy) = (x0 as f64 + cw / 2.0, y0 as f64 + ch / 2.0);
    let pad = cw.max(ch) * MARGIN;
    cw += 2.0 * pad;
    ch += 2.0 * pad;
    cw = cw.max(w as f64 / MAGNIFY);
    ch = ch.max(h as f64 / MAGNIFY);

    // To the tile's aspect, so nothing is stretched.
    let want = THUMB.0 as f64 / THUMB.1 as f64;
    if cw / ch < want {
        cw = ch * want;
    } else {
        ch = cw / want;
    }
    let (x0, y0) = (cx - cw / 2.0, cy - ch / 2.0);

    let (ow, oh) = (cw.round() as u32, ch.round() as u32);
    let mut out = vec![0u8; (ow * oh * 4) as usize];
    for y in 0..oh {
        for x in 0..ow {
            let (sx, sy) = (x0 as i64 + x as i64, y0 as i64 + y as i64);
            let k = ((y * ow + x) * 4) as usize;
            if sx >= 0 && sy >= 0 && (sx as u32) < w && (sy as u32) < h {
                let j = at(sx as u32, sy as u32);
                out[k..k + 4].copy_from_slice(&rgba[j..j + 4]);
            } else {
                out[k..k + 3].copy_from_slice(&bg);
                out[k + 3] = 255;
            }
        }
    }
    (ow, oh, out)
}

/// Box-average down to `THUMB`.
///
/// Averaged rather than sampled, because a marker one pixel wide vanishes under
/// nearest-neighbour: five bodies at 1100 points across become nothing at 240.
fn shrink(w: u32, h: u32, rgba: &[u8], tw: u32, th: u32) -> Vec<u8> {
    let mut out = vec![255u8; (tw * th * 4) as usize];
    for y in 0..th {
        let y0 = y * h / th;
        let y1 = (((y + 1) * h / th).max(y0 + 1)).min(h);
        for x in 0..tw {
            let x0 = x * w / tw;
            let x1 = (((x + 1) * w / tw).max(x0 + 1)).min(w);
            let (mut r, mut g, mut b, mut n) = (0u64, 0u64, 0u64, 0u64);
            for yy in y0..y1 {
                for xx in x0..x1 {
                    let k = ((yy * w + xx) * 4) as usize;
                    r += rgba[k] as u64;
                    g += rgba[k + 1] as u64;
                    b += rgba[k + 2] as u64;
                    n += 1;
                }
            }
            let k = ((y * tw + x) * 4) as usize;
            out[k] = (r / n) as u8;
            out[k + 1] = (g / n) as u8;
            out[k + 2] = (b / n) as u8;
        }
    }
    out
}

/// Write the pixels, as whatever the path's extension asks for.
///
/// PPM is what this wrote and is a format you can produce in six lines and read in a hex editor,
/// which is why the snapshot started there. PNG is what a picture that goes anywhere else has to
/// be, and `image` is already in this binary's tree through `eframe` — declaring it added no crate
/// to the lockfile, which is the only reason it is here.
fn write_image(path: &str, width: u32, height: u32, rgba: &[u8]) {
    let png = std::path::Path::new(path)
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("png"));
    let wrote = if png {
        // **Three channels, not four.** The pass writes an opaque frame and the alpha is 255
        // everywhere; carrying it cost 536 kB across thirty tiles against 188 kB without, for a
        // channel that says nothing.
        let (pixels, _) = rgba.as_chunks::<4>();
        let rgb: Vec<u8> = pixels.iter().flat_map(|p| p[..3].to_vec()).collect();
        // **Best compression, chosen rather than defaulted.** These are written once and read
        // from a binary that embeds them, so encoder time is free and every kilobyte is carried
        // in the repository for as long as the scene exists: the default setting cost 520 kB
        // across twenty-seven tiles against 188 kB at this one.
        std::fs::File::create(path)
            .map_err(|e| e.to_string())
            .and_then(|f| {
                image::codecs::png::PngEncoder::new_with_quality(
                    std::io::BufWriter::new(f),
                    image::codecs::png::CompressionType::Best,
                    image::codecs::png::FilterType::Adaptive,
                )
                .write_image(&rgb, width, height, image::ExtendedColorType::Rgb8)
                .map_err(|e| e.to_string())
            })
    } else {
        let mut out = format!("P6\n{width} {height}\n255\n").into_bytes();
        let (pixels, _) = rgba.as_chunks::<4>();
        for p in pixels {
            out.extend_from_slice(&p[..3]);
        }
        std::fs::write(path, out).map_err(|e| e.to_string())
    };
    if let Err(e) = wrote {
        eprintln!("cannot write {path}: {e}");
        std::process::exit(1);
    }
}

/// The clear colour, shared by the window and the snapshot so the two agree.
const BACKGROUND: wgpu::Color = wgpu::Color {
    r: 0.039,
    g: 0.055,
    b: 0.075,
    a: 1.0,
};

/// The line pipeline, built once here so the window and the snapshot cannot drift apart.
fn line_pipeline(device: &wgpu::Device, format: wgpu::TextureFormat) -> wgpu::RenderPipeline {
    let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("lines"),
        source: wgpu::ShaderSource::Wgsl(SHADER.into()),
    });
    let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: None,
        bind_group_layouts: &[],
        push_constant_ranges: &[],
    });
    device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("lines"),
        layout: Some(&layout),
        vertex: wgpu::VertexState {
            module: &shader,
            entry_point: "vs",
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<Vertex>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &[
                    wgpu::VertexAttribute {
                        offset: 0,
                        shader_location: 0,
                        format: wgpu::VertexFormat::Float32x2,
                    },
                    wgpu::VertexAttribute {
                        offset: 8,
                        shader_location: 1,
                        format: wgpu::VertexFormat::Float32x3,
                    },
                ],
            }],
            compilation_options: Default::default(),
        },
        fragment: Some(wgpu::FragmentState {
            module: &shader,
            entry_point: "fs",
            targets: &[Some(format.into())],
            compilation_options: Default::default(),
        }),
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::LineList,
            ..Default::default()
        },
        depth_stencil: None,
        multisample: wgpu::MultisampleState::default(),
        multiview: None,
        cache: None,
    })
}

/// A vertex buffer from bytes, without pulling in `wgpu::util`.
trait BufferInit {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer;
}

impl BufferInit for wgpu::Device {
    fn create_buffer_init_lite(&self, bytes: &[u8]) -> wgpu::Buffer {
        let buffer = self.create_buffer(&wgpu::BufferDescriptor {
            label: Some("lines"),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::VERTEX,
            mapped_at_creation: true,
        });
        buffer
            .slice(..)
            .get_mapped_range_mut()
            .copy_from_slice(bytes);
        buffer.unmap();
        buffer
    }
}

/// The same ramp the HTML report uses, so the two views of one run agree about colour.
fn ramp(t: f64) -> [f32; 3] {
    let t = t.clamp(0.0, 1.0);
    let stops = [
        [24.0, 52.0, 110.0],
        [61.0, 139.0, 255.0],
        [255.0, 194.0, 71.0],
        [255.0, 82.0, 33.0],
    ];
    let (a, b, u) = if t < 0.34 {
        (stops[0], stops[1], t / 0.34)
    } else if t < 0.67 {
        (stops[1], stops[2], (t - 0.34) / 0.33)
    } else {
        (stops[2], stops[3], (t - 0.67) / 0.33)
    };
    [
        ((a[0] + (b[0] - a[0]) * u) / 255.0) as f32,
        ((a[1] + (b[1] - a[1]) * u) / 255.0) as f32,
        ((a[2] + (b[2] - a[2]) * u) / 255.0) as f32,
    ]
}

const SHADER: &str = r#"
struct In { @location(0) pos: vec2<f32>, @location(1) colour: vec3<f32> };
struct Out { @builtin(position) clip: vec4<f32>, @location(0) colour: vec3<f32> };

@vertex
fn vs(v: In) -> Out {
    var out: Out;
    out.clip = vec4<f32>(v.pos, 0.0, 1.0);
    out.colour = v.colour;
    return out;
}

@fragment
fn fs(v: Out) -> @location(0) vec4<f32> {
    return vec4<f32>(v.colour, 1.0);
}
"#;
