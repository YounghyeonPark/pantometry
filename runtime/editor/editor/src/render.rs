//! **Shaded surfaces with a depth buffer, in the viewport, on the GPU.**
//!
//! What this replaces is the reason it exists. The viewport drew a field as translucent circles
//! sorted far to near — a volume render by painter's algorithm — and a body as a filled circle. It
//! reads as a point cloud, which is the same criticism the HTML report's first draft earned, and it
//! is not what a reader comparing this against usdview or 3ds Max sees there.
//!
//! # Three things it is careful about
//!
//! **The geometry is `pantometry_view::mesh`.** The same function that writes glTF and USD builds
//! what is on screen here, so the editor cannot disagree with an export about how big a solid is or
//! where its boundary lies. That arithmetic has already been wrong once — a 40 mm cube exported 80
//! mm across, because a field sampled corner to corner means an end node owns half a cell — and the
//! fix belongs in one place.
//!
//! **The colour is `editor_core::Colouring`**, which is Planck's law where the field is hot enough
//! to glow and `pantometry_view::ramp` where it is not — the same decision the flat painter made,
//! now made once. It arrives here already in **linear** RGB, because a shader interpolates across a
//! triangle and multiplies by a light, and both of those are wrong in sRGB: the glTF exporter
//! shipped that defect and every export came out about 2.3× too bright in the midtones.
//!
//! **Vertices are framing-local.** `Framing::local` subtracts the centre in `f64` and only the
//! result narrows, so a 9 mm block sitting 200 mm from the origin arrives as numbers of order one.
//! Folding the centre into the matrix instead was measured disagreeing with `Camera::project` by
//! `4.4e-6` — `f32` keeping two digits of seven out of a subtraction.

use eframe::glow::{self, HasContext};

/// What the viewport owns on the GPU, shared with the paint callback.
///
/// `egui::PaintCallback` holds an `Arc<dyn Any + Send + Sync>`, so the meshes cannot live in `App`
/// and be borrowed into the closure — they live here behind one lock, and the closure and the frame
/// builder take turns. The GL context does not exist until a callback runs, which is also why the
/// meshes are built there rather than at start-up.
#[derive(Default)]
pub struct Shared {
    solid: Option<Mesh>,
    lines: Option<Mesh>,
    /// Uploaded on the next callback, then dropped. `None` means what is on the GPU is current.
    pub pending_solid: Option<Batch>,
    /// As above, for the line pass.
    pub pending_lines: Option<Batch>,
    /// Why there is no shaded viewport, if there is not.
    ///
    /// A machine whose driver will not give a 3.3 core context should lose the surfaces and keep the
    /// editor — and should be **told which of those happened**, in the viewport, rather than shown
    /// an empty rectangle that looks exactly like a scene with nothing in it.
    pub error: Option<String>,
    /// How many triangles and lines the last paint actually drew, and how many paints there have
    /// been.
    ///
    /// A statistics readout, which every DCC viewport has — and the only way to tell a viewport
    /// drawing nothing from a viewport that was never asked to. The first version of this pass drew
    /// nothing and looked identical to a callback that had not run; there was no logger installed,
    /// so `egui_glow`'s own warning about an unrecognised callback went nowhere.
    pub drawn: (usize, usize, u64),
}

impl Shared {
    /// Upload anything pending and draw both passes.
    ///
    /// `clip` is column-major from [`viewer_core::Camera::matrix`], and `pixels` is the callback's
    /// own viewport, which the depth clear is scissored to.
    pub fn paint(&mut self, gl: &glow::Context, clip: &[f32; 16], pixels: [i32; 4]) {
        if self.error.is_some() {
            return;
        }
        for (slot, pending, lines) in [
            (0usize, self.pending_solid.take(), false),
            (1, self.pending_lines.take(), true),
        ] {
            let held = if slot == 0 {
                &mut self.solid
            } else {
                &mut self.lines
            };
            if held.is_none() {
                match Mesh::new(gl, lines) {
                    Ok(m) => *held = Some(m),
                    Err(why) => {
                        // On stderr as well as in the struct: losing the 3D view of a 3D editor is
                        // not something to mention only in a status line somebody may not read.
                        eprintln!("viewport: {why}");
                        self.error = Some(why);
                        return;
                    }
                }
            }
            if let (Some(mesh), Some(batch)) = (held.as_mut(), pending) {
                mesh.upload(gl, &batch);
            }
        }

        // **The depth buffer is the whole point, and egui does not want one.** egui paints in order
        // with depth testing off, so it has to be turned on here and turned back off, or every
        // widget painted after this frame fails its test against geometry it knows nothing about.
        // Only the depth is cleared, and only inside this rect: the colour here is whatever egui
        // has already painted underneath, and `glClear` obeys the scissor rather than the viewport.
        unsafe {
            gl.enable(glow::SCISSOR_TEST);
            gl.scissor(pixels[0], pixels[1], pixels[2], pixels[3]);
            gl.enable(glow::DEPTH_TEST);
            gl.depth_func(glow::LESS);
            gl.depth_mask(true);
            gl.clear_depth_f64(1.0);
            gl.clear(glow::DEPTH_BUFFER_BIT);
            // Off, because a slice through a voxel boundary shows its inside faces and a solid with
            // its far side culled reads as a shell of the wrong thickness.
            gl.disable(glow::CULL_FACE);
        }

        for mesh in [self.solid.as_ref(), self.lines.as_ref()]
            .into_iter()
            .flatten()
        {
            mesh.draw(gl, clip);
        }
        self.drawn = (
            self.solid.as_ref().map_or(0, |m| m.count as usize / 3),
            self.lines.as_ref().map_or(0, |m| m.count as usize / 2),
            self.drawn.2 + 1,
        );
        // **A way to ask the pass what it did, from outside.** `PANTOMETRY_VIEWPORT` prints the
        // counts and the viewport rect on every paint. Written because a viewport that drew nothing
        // and a viewport that was never asked to run look identical from a screenshot, and there is
        // no logger in this binary for `egui_glow`'s own warnings to reach.
        if std::env::var_os("PANTOMETRY_VIEWPORT").is_some() {
            eprintln!(
                "viewport: {} triangles, {} lines, paint {}, scissor {:?}",
                self.drawn.0, self.drawn.1, self.drawn.2, pixels
            );
        }
        if let Some(path) = std::env::var_os("PANTOMETRY_VIEWPORT_SHOT") {
            match self.snapshot(gl, pixels, std::path::Path::new(&path)) {
                Ok(lit) => eprintln!("viewport: wrote {path:?}, {lit} pixels not background"),
                Err(why) => eprintln!("viewport: no snapshot: {why}"),
            }
        }

        unsafe {
            gl.disable(glow::DEPTH_TEST);
            gl.depth_mask(false);
        }
    }

    /// Read the viewport back out of the framebuffer and write it as a PPM.
    ///
    /// **The point of this is that it is not a screenshot.** Capturing the window from outside gets
    /// whatever the desktop compositor feels like handing over — a DPI-virtualised `PrintWindow`
    /// returned a plausible image of the panels with the 3D content simply absent, and a screen grab
    /// got a different application's window, because a background process cannot raise a window on
    /// Windows. This reads the pixels the pass itself just wrote, which is the only thing that
    /// answers "did the geometry land where the matrix said".
    ///
    /// A PPM because it is eleven bytes of header and no dependency. Returns how many pixels differ
    /// from the corner pixel, which is a one-number answer to "is anything there".
    fn snapshot(
        &self,
        gl: &glow::Context,
        pixels: [i32; 4],
        path: &std::path::Path,
    ) -> Result<usize, String> {
        let (w, h) = (pixels[2].max(1) as usize, pixels[3].max(1) as usize);
        let mut buf = vec![0u8; w * h * 4];
        unsafe {
            gl.read_pixels(
                pixels[0],
                pixels[1],
                w as i32,
                h as i32,
                glow::RGBA,
                glow::UNSIGNED_BYTE,
                glow::PixelPackData::Slice(Some(&mut buf)),
            );
        }

        // GL reads bottom-up; a PPM is written top-down.
        let mut out = format!("P6\n{w} {h}\n255\n").into_bytes();
        let corner = [buf[0], buf[1], buf[2]];
        let mut lit = 0usize;
        for row in (0..h).rev() {
            for col in 0..w {
                let i = 4 * (row * w + col);
                let rgb = [buf[i], buf[i + 1], buf[i + 2]];
                // "Not background" against the corner pixel rather than against white: the editor
                // follows the system theme and the background is whatever egui painted.
                if rgb.iter().zip(&corner).any(|(a, b)| a.abs_diff(*b) > 6) {
                    lit += 1;
                }
                out.extend(rgb);
            }
        }
        std::fs::write(path, out).map_err(|e| e.to_string())?;
        Ok(lit)
    }

    /// Give everything back. Called from `eframe::App::on_exit`, which is where the context still
    /// exists — dropping these without a context leaks the driver's objects for the process's life.
    pub fn destroy(&self, gl: &glow::Context) {
        for mesh in [self.solid.as_ref(), self.lines.as_ref()]
            .into_iter()
            .flatten()
        {
            mesh.destroy(gl);
        }
    }
}

/// One mesh on the GPU: triangles, or lines.
struct Mesh {
    program: glow::Program,
    vao: glow::VertexArray,
    buffers: [glow::Buffer; 4],
    /// How many indices to draw. Zero is a valid state — a run with nothing present — and draws
    /// nothing rather than being an error.
    count: i32,
    lines: bool,
}

/// What to draw, with a colour already resolved for every vertex.
///
/// Deliberately not `pantometry_view::mesh::Surface`: lines arrive here too, and the caller has
/// already placed every value on the run-wide scale. Keeping the colour out would mean the shader
/// knowing about scales, which is what `editor_core::Colouring` exists to own.
pub struct Batch {
    /// Framing-local positions, three floats a vertex.
    pub positions: Vec<f32>,
    /// Unit normals, three floats a vertex. Empty for a line batch.
    pub normals: Vec<f32>,
    /// Linear RGB, three floats a vertex.
    pub colours: Vec<f32>,
    /// Triangles, or vertex pairs for a line batch.
    pub indices: Vec<u32>,
}

impl Batch {
    /// An empty batch, which draws nothing.
    pub fn new() -> Batch {
        Batch {
            positions: Vec::new(),
            normals: Vec::new(),
            colours: Vec::new(),
            indices: Vec::new(),
        }
    }

    /// How many vertices are in it, which is what the next index has to be offset by.
    pub fn vertices(&self) -> u32 {
        (self.positions.len() / 3) as u32
    }

    /// Add one vertex.
    pub fn push(&mut self, at: [f32; 3], normal: [f32; 3], colour: [f32; 3]) {
        self.positions.extend(at);
        self.normals.extend(normal);
        self.colours.extend(colour);
    }
}

impl Default for Batch {
    fn default() -> Batch {
        Batch::new()
    }
}

const SHADER_VERTEX: &str = r#"#version 330 core
layout(location = 0) in vec3 a_position;
layout(location = 1) in vec3 a_normal;
layout(location = 2) in vec3 a_colour;

uniform mat4 u_clip;

out vec3 v_normal;
out vec3 v_colour;

void main() {
    // The normal is rotated by the same matrix's linear part. A framing is one uniform scale on all
    // three axes -- `Framing` keeps a single span precisely so that a cube is drawn as a cube -- so
    // the inverse transpose is that matrix up to a factor, and normalising absorbs the factor.
    v_normal = mat3(u_clip) * a_normal;
    v_colour = a_colour;
    gl_Position = u_clip * vec4(a_position, 1.0);
}
"#;

/// Two lights and an ambient, which is the smallest arrangement that shows a shape: a key over the
/// viewer's shoulder and a fill from the other side and below, so a face turned away is dark rather
/// than black. One light leaves half of every object unreadable; none is the flat-coloured blob the
/// splats already were.
const SHADER_FRAGMENT: &str = r#"#version 330 core
in vec3 v_normal;
in vec3 v_colour;

uniform int u_lit;

out vec4 f_colour;

void main() {
    vec3 shade = vec3(1.0);
    if (u_lit != 0) {
        // Two-sided. A voxel boundary is closed and outward-facing, but a slice through one is not,
        // and a surface lit from behind should read as turned away rather than as absent -- an unlit
        // black face is indistinguishable from a hole, and a hole is a physical claim here.
        vec3 n = normalize(v_normal);
        if (!gl_FrontFacing) {
            n = -n;
        }
        vec3 key = normalize(vec3(-0.35, 0.55, 1.0));
        vec3 fill = normalize(vec3(0.6, -0.4, 0.4));
        float l = 0.30 + 0.62 * max(dot(n, key), 0.0) + 0.22 * max(dot(n, fill), 0.0);
        shade = vec3(min(l, 1.25));
    }
    vec3 linear = clamp(v_colour * shade, 0.0, 1.0);
    // egui's glow painter draws into a framebuffer it treats as sRGB, so the encode belongs here.
    // Leaving it out is the glTF exporter's defect in the other direction: linear written where
    // sRGB is expected comes out dark and muddy rather than too bright.
    f_colour = vec4(pow(linear, vec3(1.0 / 2.2)), 1.0);
}
"#;

impl Mesh {
    /// Compile the program and allocate the buffers, or say why not.
    fn new(gl: &glow::Context, lines: bool) -> Result<Mesh, String> {
        unsafe {
            let program = gl
                .create_program()
                .map_err(|e| format!("no shader program: {e}"))?;
            let mut shaders = Vec::new();
            for (kind, source) in [
                (glow::VERTEX_SHADER, SHADER_VERTEX),
                (glow::FRAGMENT_SHADER, SHADER_FRAGMENT),
            ] {
                let s = gl
                    .create_shader(kind)
                    .map_err(|e| format!("no shader: {e}"))?;
                gl.shader_source(s, source);
                gl.compile_shader(s);
                if !gl.get_shader_compile_status(s) {
                    return Err(format!(
                        "the viewport shader did not compile: {}",
                        gl.get_shader_info_log(s)
                    ));
                }
                gl.attach_shader(program, s);
                shaders.push(s);
            }
            gl.link_program(program);
            if !gl.get_program_link_status(program) {
                return Err(format!(
                    "the viewport shader did not link: {}",
                    gl.get_program_info_log(program)
                ));
            }
            for s in shaders {
                gl.detach_shader(program, s);
                gl.delete_shader(s);
            }

            let vao = gl
                .create_vertex_array()
                .map_err(|e| format!("no vertex array: {e}"))?;
            // Four separate calls rather than an array default: a `glow::Buffer` is the driver's
            // name for an object and there is no such thing as a blank one.
            let mut made = Vec::with_capacity(4);
            for _ in 0..4 {
                made.push(gl.create_buffer().map_err(|e| format!("no buffer: {e}"))?);
            }
            let buffers: [glow::Buffer; 4] = [made[0], made[1], made[2], made[3]];
            Ok(Mesh {
                program,
                vao,
                buffers,
                count: 0,
                lines,
            })
        }
    }

    /// Replace what is on the GPU with `batch`.
    ///
    /// Called when the frame, the visibility set or the scene changes — not every paint. A drag
    /// redraws sixty times a second and re-uploading a boundary each time is the difference between
    /// a viewport you can aim and one you fight.
    fn upload(&mut self, gl: &glow::Context, batch: &Batch) {
        // A normal per vertex either way: an unbound attribute reads whatever the last upload left
        // there, and for a line batch that would be somebody else's geometry.
        let zeros = vec![0.0f32; batch.positions.len()];
        let normals = if batch.normals.len() == batch.positions.len() {
            &batch.normals
        } else {
            &zeros
        };
        unsafe {
            gl.bind_vertex_array(Some(self.vao));
            for (index, data) in [(0u32, &batch.positions), (1, normals), (2, &batch.colours)] {
                gl.bind_buffer(glow::ARRAY_BUFFER, Some(self.buffers[index as usize]));
                gl.buffer_data_u8_slice(glow::ARRAY_BUFFER, as_bytes(data), glow::DYNAMIC_DRAW);
                gl.enable_vertex_attrib_array(index);
                gl.vertex_attrib_pointer_f32(index, 3, glow::FLOAT, false, 12, 0);
            }
            gl.bind_buffer(glow::ELEMENT_ARRAY_BUFFER, Some(self.buffers[3]));
            gl.buffer_data_u8_slice(
                glow::ELEMENT_ARRAY_BUFFER,
                as_bytes(&batch.indices),
                glow::DYNAMIC_DRAW,
            );
            gl.bind_vertex_array(None);
        }
        self.count = batch.indices.len() as i32;
    }

    /// Draw it. Depth state is the caller's, because both passes share one buffer.
    fn draw(&self, gl: &glow::Context, clip: &[f32; 16]) {
        if self.count == 0 {
            return;
        }
        unsafe {
            gl.use_program(Some(self.program));
            if let Some(u) = gl.get_uniform_location(self.program, "u_clip") {
                gl.uniform_matrix_4_f32_slice(Some(&u), false, clip);
            }
            if let Some(u) = gl.get_uniform_location(self.program, "u_lit") {
                gl.uniform_1_i32(Some(&u), i32::from(!self.lines));
            }
            gl.bind_vertex_array(Some(self.vao));
            let mode = if self.lines {
                glow::LINES
            } else {
                glow::TRIANGLES
            };
            gl.draw_elements(mode, self.count, glow::UNSIGNED_INT, 0);
            gl.bind_vertex_array(None);
        }
    }

    /// Give the buffers and the program back.
    fn destroy(&self, gl: &glow::Context) {
        unsafe {
            gl.delete_program(self.program);
            gl.delete_vertex_array(self.vao);
            for b in self.buffers {
                gl.delete_buffer(b);
            }
        }
    }
}

/// A slice of plain numbers as the bytes the driver wants, in native order.
fn as_bytes<T: Copy>(v: &[T]) -> &[u8] {
    // SAFETY: called only with `f32` and `u32`, which have no padding and no invalid bit patterns.
    // The result borrows `v`, so it cannot outlive it, and the length is exact.
    unsafe { std::slice::from_raw_parts(v.as_ptr() as *const u8, std::mem::size_of_val(v)) }
}
