//! A frame as glTF 2.0, so somebody else's renderer can draw it.
//!
//! # Why export rather than render
//!
//! Blender, three.js, Omniverse, Isaac Sim, macOS Quick Look and every USD pipeline read glTF.
//! Competing with those on rendering is a losing use of effort; what this workspace has that they
//! do not is physics that is audited, deterministic and checked against closed forms. So the move
//! is to make a result *reachable from* them.
//!
//! # It costs nothing
//!
//! glTF is JSON with the binary data base64'd into a `data:` URI, and this crate already writes
//! JSON by hand on purpose. So there is no new dependency, no encoder, and no build step — the
//! same reason SVG was chosen over a raster format. `pantometry-view` still has exactly one
//! dependency, and it is `pantometry-scene`.
//!
//! # Surfaces, not points
//!
//! This exported point clouds. A field became one vertex per cell, and a point cloud in Blender
//! has no silhouette, takes no light and casts no shadow: it is a picture of the *sampling*. A
//! solid is a solid, and a solid has a surface.
//!
//! A three-dimensional field is now **the surface of its material**, coloured by the field on it.
//! One quad per cell face whose neighbour is absent or outside the grid, so the interior is not
//! emitted and a void inside the block produces a real interior surface rather than a hole in the
//! colours. Measured on a solid 9×9×9 block: 486 quads against the 4,374 a per-cell box would
//! write, because 89% of the faces of a solid are between two cells and cannot be seen.
//!
//! Flat-shaded on purpose — four vertices per quad carrying that face's normal — because a voxel
//! surface *is* faceted and smoothing it would draw a rounded object the simulation does not
//! have.
//!
//! **A surface hides the inside, and that is a real loss.** A hot spot in the middle of a block
//! exports as a block whose faces are all at ambient, because they are: the point cloud showed the
//! interior and this does not. Rendered from Blender, scene 15 is a uniformly cold cube, which is
//! a true picture of its surface and says nothing about the 353 K cell at its centre. The
//! interior is what the HTML report's raycast and slice montage are for, and the glTF is for the
//! object. Neither substitutes and both are written from the same run.
//!
//! Where a void makes the interior *be* a surface, this shows it: scene 23 is a hot part and a
//! cooled lid with a real gap between them, and the two solids come out as two solids.
//!
//! One thing to expect on import: a millimetre-scale part sits inside Blender's default camera
//! clip of 0.1 m, so the first render of a 9 mm block was an empty frame with the geometry present
//! and correct. Set the near plane from the model, not from the default.
//!
//! # What each shape becomes
//!
//! | panel | glTF |
//! | --- | --- |
//! | a 3D field | `TRIANGLES` — the surface of the present cells, with normals |
//! | a 2D field | `TRIANGLES` — a plane of quads over the extent it was sampled on |
//! | a 1D field | nothing — a row of samples along a line is a graph, not geometry |
//! | points | `TRIANGLES` — a sphere each, up to [`MAX_SPHERES`], then `POINTS` |
//! | paths | `LINES`, one segment per adjacent pair, coloured per path |
//!
//! The 1D row is a refusal rather than an omission, and [`gltf`] reports what it left out. The 2D
//! row used to be one too, with the same reasoning — and the reasoning had expired: a plane of
//! samples was a graph only because the run did not record the box it was sampled over, so there
//! was no size to give it. `PanelData::Field` carries `extent_m` now, and a 2D field over a real
//! box is a plate, a floor or a wall.
//!
//! # One frame, and why not the run
//!
//! [`gltf`] takes a single [`Frame`]. glTF animates **node transforms and morph targets** — a
//! thing moving or a mesh deforming between fixed vertex counts. A field whose values change, a
//! ray bundle that retraces, a body count that is constant only by luck: none of those is either.
//! Encoding a run as an animation would mean choosing a lie about what is moving.
//!
//! Exporting the frame a caller asks for is the honest shape. For a whole run there is
//! [`to_json`](crate::to_json), which carries every frame and is what the native viewer reads.
//!
//! One consequence to know about: the colour scale here spans **this frame**, because one frame is
//! all there is. Every other view in this crate spans the run, and two frames exported separately
//! are therefore not on the same scale. [`Exported::notes`] says so on every export that has a
//! colour in it.

use crate::ramp;
use pantometry_scene::{Frame, Panel, PanelData};

/// The most bodies that become spheres before the exporter falls back to a point cloud.
///
/// A sphere is 96 quads, so 256 bodies is about 49,000 triangles — a comfortable Blender scene.
/// Ten thousand bodies would be two million, which is a worse picture than points and a much
/// slower one.
pub const MAX_SPHERES: usize = 256;

/// The most cell faces one field may write.
///
/// Past this the field is subsampled with a stride and [`Exported::notes`] says so, because a
/// silently decimated surface is a picture of a coarser simulation than the one that ran. A
/// 100³ field has a 60,000-face surface, which is fine; a 100³ field *full of holes* can have
/// far more, which is what this is really for.
pub const MAX_FACES: usize = 200_000;

/// One frame as a glTF 2.0 document, and what was left out of it.
pub struct Exported {
    /// The document. Write it to a `.gltf` file; the binary data is embedded, so there is no
    /// sidecar to lose.
    pub document: String,
    /// Panels that produced no geometry, and why — so a caller is never left wondering whether an
    /// empty scene is the physics or the exporter.
    pub skipped: Vec<String>,
    /// Choices this export made that a reader should know about: a size that is not in the data,
    /// a subsampling, the scale being this frame's.
    ///
    /// Separate from [`Exported::skipped`], which is about what is *not* in the file. These are
    /// about what is, and on what terms.
    pub notes: Vec<String>,
}

/// A colour for a value in `0..=1`, as glTF wants it.
///
/// **Linear, not sRGB, and that was a real defect.** glTF 2.0 specifies `COLOR_0` and
/// `baseColorFactor` in *linear* space; only textures are sRGB-encoded. This wrote `byte / 255`
/// straight from an sRGB ramp, so every colour that ever left this workspace was decoded as
/// though it were already linear — which lightens and desaturates it. A mid-grey `#808080` is
/// 0.5 in sRGB and **0.216** in linear: the export was handing renderers a value 2.3x too bright
/// in the midtones, uniformly, and it looked like a plausible picture the whole time.
///
/// The ramp itself is [`crate::ramp`], so the export, the HTML report and both editor shells are
/// one scale. What was here was a fifth copy of the four-stop gradient whose lightness folds back
/// on itself at 0.67.
fn colour(t: f64, signed: bool) -> [f32; 4] {
    let rgb = if signed {
        ramp::diverging(t)
    } else {
        ramp::sequential(t)
    };
    let linear = |b: u8| {
        let c = b as f64 / 255.0;
        let v = if c <= 0.040_45 {
            c / 12.92
        } else {
            ((c + 0.055) / 1.055).powf(2.4)
        };
        v as f32
    };
    [linear(rgb[0]), linear(rgb[1]), linear(rgb[2]), 1.0]
}

/// One mesh under construction.
struct Mesh {
    name: String,
    positions: Vec<[f32; 3]>,
    /// Empty for a point cloud and for lines, which take no light.
    normals: Vec<[f32; 3]>,
    colours: Vec<[f32; 4]>,
    /// Empty for a point cloud.
    indices: Vec<u32>,
    /// glTF primitive mode: 0 points, 1 lines, 4 triangles.
    mode: u32,
}

impl Mesh {
    fn new(name: &str, mode: u32) -> Mesh {
        Mesh {
            name: name.to_string(),
            positions: Vec::new(),
            normals: Vec::new(),
            colours: Vec::new(),
            indices: Vec::new(),
            mode,
        }
    }

    /// Add one flat-shaded quad: four corners in winding order, one normal, one colour.
    ///
    /// Its own four vertices rather than shared ones, because the normal belongs to the face. A
    /// voxel surface sharing corner vertices would average three perpendicular normals into a
    /// direction no face points, and shade a cube as a blob.
    fn quad(&mut self, corners: [[f32; 3]; 4], normal: [f32; 3], colour: [f32; 4]) {
        let base = self.positions.len() as u32;
        for c in corners {
            self.positions.push(c);
            self.normals.push(normal);
            self.colours.push(colour);
        }
        self.indices
            .extend([base, base + 1, base + 2, base, base + 2, base + 3]);
    }
}

/// Export one frame.
///
/// Every panel that is geometry becomes a node under one scene, named after its domain, so a
/// reader opening the file sees the same names the run reported.
pub fn gltf(title: &str, frame: &Frame) -> Exported {
    let mut meshes = Vec::new();
    let mut skipped = Vec::new();
    let mut notes = Vec::new();
    let mut coloured = false;

    for panel in &frame.panels {
        match &panel.data {
            PanelData::Paths {
                vertices,
                starts,
                values,
                ..
            } => {
                let (lo, hi, signed) = span(values);
                coloured = true;
                let mut mesh = Mesh::new(&panel.name, 1);
                mesh.positions.reserve(vertices.len());
                for (k, start) in starts.iter().enumerate() {
                    let from = *start;
                    let to = starts.get(k + 1).copied().unwrap_or(vertices.len());
                    let c = colour(
                        place(values.get(k).copied().unwrap_or(lo), lo, hi, signed),
                        signed,
                    );
                    let base = mesh.positions.len() as u32;
                    for v in &vertices[from..to] {
                        mesh.positions.push([v[0] as f32, v[1] as f32, v[2] as f32]);
                        mesh.colours.push(c);
                    }
                    // One segment per adjacent pair. Indices rather than duplicated vertices, so a
                    // ray of n points costs n positions and not 2(n-1).
                    for i in 0..(to - from).saturating_sub(1) as u32 {
                        mesh.indices.push(base + i);
                        mesh.indices.push(base + i + 1);
                    }
                }
                if !mesh.indices.is_empty() {
                    meshes.push(mesh);
                }
            }
            PanelData::Points {
                positions,
                values,
                bounds,
                ..
            } => {
                let (lo, hi, signed) = span(values);
                coloured = true;
                if positions.len() <= MAX_SPHERES {
                    let r = body_radius(positions, bounds);
                    notes.push(format!(
                        "{}: bodies drawn as spheres of radius {r:.6} m — a **size this run does \
                         not carry**. A body set records positions and a value, not an extent, so \
                         the radius is a quarter of the median distance to the nearest neighbour, \
                         which keeps them apart and is a drawing convention rather than a \
                         measurement",
                        panel.name
                    ));
                    let mut mesh = Mesh::new(&panel.name, 4);
                    for (i, p) in positions.iter().enumerate() {
                        let c = colour(
                            place(values.get(i).copied().unwrap_or(lo), lo, hi, signed),
                            signed,
                        );
                        sphere(&mut mesh, [p[0], p[1], p[2]], r, c);
                    }
                    meshes.push(mesh);
                } else {
                    notes.push(format!(
                        "{}: {} bodies is over the {MAX_SPHERES} a sphere each is worth, so this \
                         is a point cloud — which takes no light and casts no shadow",
                        panel.name,
                        positions.len()
                    ));
                    let mut mesh = Mesh::new(&panel.name, 0);
                    for (i, p) in positions.iter().enumerate() {
                        mesh.positions.push([p[0] as f32, p[1] as f32, p[2] as f32]);
                        mesh.colours.push(colour(
                            place(values.get(i).copied().unwrap_or(lo), lo, hi, signed),
                            signed,
                        ));
                    }
                    meshes.push(mesh);
                }
            }
            PanelData::Field {
                nx,
                ny,
                nz,
                extent_m,
                values,
            } => {
                let dims = [*nx, *ny, *nz].iter().filter(|&&n| n > 1).count();
                if dims < 2 {
                    skipped.push(format!(
                        "{} is a {nx}x{ny}x{nz} field: a row of samples along a line is a graph, \
                         not geometry, so it is not in the scene",
                        panel.name
                    ));
                    continue;
                }
                let (lo, hi, signed) = span(values);
                coloured = true;
                let mut mesh = Mesh::new(&panel.name, 4);
                let stride = surface(
                    &mut mesh,
                    (*nx, *ny, *nz),
                    *extent_m,
                    values,
                    lo,
                    hi,
                    signed,
                );
                if stride > 1 {
                    notes.push(format!(
                        "{}: subsampled every {stride} cells — the full surface is over \
                         {MAX_FACES} faces",
                        panel.name
                    ));
                }
                if mesh.indices.is_empty() {
                    skipped.push(format!(
                        "{} is a {nx}x{ny}x{nz} field with no cell that holds a value, so it has \
                         no surface",
                        panel.name
                    ));
                } else {
                    meshes.push(mesh);
                }
            }
        }
    }

    if coloured {
        notes.push(
            "the colour scale spans **this frame**, because one frame is all a glTF export \
             carries — two frames exported separately are not on the same scale"
                .to_string(),
        );
    }

    Exported {
        document: document(title, &meshes),
        skipped,
        notes,
    }
}

/// Where a value sits on its colour scale, in `0..=1`.
///
/// A signed range is placed about **zero** rather than about the middle of the range, so the
/// neutral colour lands on the value zero. For −100 to +300 those are a quarter of the scale
/// apart.
fn place(v: f64, lo: f64, hi: f64, signed: bool) -> f64 {
    if signed {
        let reach = lo.abs().max(hi.abs()).max(f64::MIN_POSITIVE);
        0.5 + 0.5 * v / reach
    } else {
        (v - lo) / (hi - lo)
    }
}

/// The surface of the cells that hold a value, coloured by the field on it.
///
/// Returns the stride used. A face is written when the neighbour across it is absent or outside
/// the grid, which is what makes this the surface of the material rather than a box per cell: on a
/// solid block 89% of the faces are between two present cells and no camera can see them.
///
/// A collapsed axis — one sample — is a plane, and gets its two faces so it is visible from both
/// sides. A plate is not invisible edge-on.
#[allow(clippy::too_many_arguments)]
fn surface(
    mesh: &mut Mesh,
    counts: (usize, usize, usize),
    extent: [f64; 6],
    values: &[f64],
    lo: f64,
    hi: f64,
    signed: bool,
) -> usize {
    let (nx, ny, nz) = counts;
    // Stride so the face count stays inside the budget. Cubic, because faces go as the square of
    // the stride reduction on each of three axes and the surface is two-dimensional; the estimate
    // only has to be in the right decade to keep the file openable.
    let faces = 2 * (nx * ny + ny * nz + nz * nx);
    let stride = if faces > MAX_FACES {
        ((faces as f64 / MAX_FACES as f64).sqrt().ceil() as usize).max(1)
    } else {
        1
    };
    let taken = |n: usize| n.div_ceil(stride);
    let (tx, ty, tz) = (taken(nx), taken(ny), taken(nz));
    let at = |i: usize, j: usize, k: usize| -> Option<f64> {
        let v = *values.get(i + nx * (j + ny * k))?;
        v.is_finite().then_some(v)
    };
    // Sample coordinates: `capture` samples corner to corner, and a collapsed axis at its middle.
    let along = |i: usize, n: usize| {
        if n > 1 {
            i as f64 / (n - 1) as f64
        } else {
            0.5
        }
    };
    let size = [
        extent[3] - extent[0],
        extent[4] - extent[1],
        extent[5] - extent[2],
    ];
    // Half a cell on each axis, which is what makes a cell a box rather than a point. A collapsed
    // axis has no cell to be half of, so it gets a sliver: thin enough to read as a plane and
    // non-zero so it is not degenerate geometry.
    let half = |a: usize, n: usize| {
        if n > 1 {
            0.5 * size[a] / (n - 1) as f64 * stride as f64
        } else {
            let widest = size[0].abs().max(size[1].abs()).max(size[2].abs());
            (widest * 1e-3).max(f64::MIN_POSITIVE)
        }
    };
    let (hx, hy, hz) = (half(0, nx), half(1, ny), half(2, nz));
    let centre = |i: usize, j: usize, k: usize| {
        [
            extent[0] + size[0] * along(i, nx),
            extent[1] + size[1] * along(j, ny),
            extent[2] + size[2] * along(k, nz),
        ]
    };

    for kt in 0..tz {
        for jt in 0..ty {
            for it in 0..tx {
                let (i, j, k) = (it * stride, jt * stride, kt * stride);
                let Some(v) = at(i, j, k) else { continue };
                let c = colour(place(v, lo, hi, signed), signed);
                let m = centre(i, j, k);
                // **Clipped to the extent.** `capture` samples corner to corner, so the first and
                // last samples sit *on* the faces of the box rather than half a cell inside it.
                // Giving every sample a full cell either side made the exported object one whole
                // cell larger than the extent — 2x for a two-sample axis, and a 40 mm cube came
                // out 80 mm. An end node owns half a cell, which is the same arithmetic the
                // boundary defect in `Tube` and `Room` turned on, now for the third time and in
                // geometry rather than in a divergence.
                let clip = |lo: f64, hi: f64, a: usize| -> (f64, f64) {
                    if size[a].abs() > 0.0 {
                        let (elo, ehi) =
                            (extent[a].min(extent[a + 3]), extent[a].max(extent[a + 3]));
                        (lo.max(elo), hi.min(ehi))
                    } else {
                        (lo, hi)
                    }
                };
                let (cx0, cx1) = clip(m[0] - hx, m[0] + hx, 0);
                let (cy0, cy1) = clip(m[1] - hy, m[1] + hy, 1);
                let (cz0, cz1) = clip(m[2] - hz, m[2] + hz, 2);
                let (x0, x1) = (cx0 as f32, cx1 as f32);
                let (y0, y1) = (cy0 as f32, cy1 as f32);
                let (z0, z1) = (cz0 as f32, cz1 as f32);
                // A neighbour one stride away, in the strided lattice — the cells between are not
                // drawn, so a face is exposed when the *drawn* neighbour is missing.
                let exposed = |di: isize, dj: isize, dk: isize| -> bool {
                    let step = stride as isize;
                    let (ni, nj, nk) = (
                        i as isize + di * step,
                        j as isize + dj * step,
                        k as isize + dk * step,
                    );
                    if ni < 0 || nj < 0 || nk < 0 {
                        return true;
                    }
                    let (ni, nj, nk) = (ni as usize, nj as usize, nk as usize);
                    if ni >= nx || nj >= ny || nk >= nz {
                        return true;
                    }
                    at(ni, nj, nk).is_none()
                };
                // Winding is counter-clockwise seen from outside, which is what glTF's front face
                // is and what a renderer culling backfaces needs.
                if exposed(0, 0, -1) {
                    mesh.quad(
                        [[x0, y0, z0], [x0, y1, z0], [x1, y1, z0], [x1, y0, z0]],
                        [0.0, 0.0, -1.0],
                        c,
                    );
                }
                if exposed(0, 0, 1) {
                    mesh.quad(
                        [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]],
                        [0.0, 0.0, 1.0],
                        c,
                    );
                }
                if exposed(0, -1, 0) {
                    mesh.quad(
                        [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]],
                        [0.0, -1.0, 0.0],
                        c,
                    );
                }
                if exposed(0, 1, 0) {
                    mesh.quad(
                        [[x0, y1, z0], [x0, y1, z1], [x1, y1, z1], [x1, y1, z0]],
                        [0.0, 1.0, 0.0],
                        c,
                    );
                }
                if exposed(-1, 0, 0) {
                    mesh.quad(
                        [[x0, y0, z0], [x0, y0, z1], [x0, y1, z1], [x0, y1, z0]],
                        [-1.0, 0.0, 0.0],
                        c,
                    );
                }
                if exposed(1, 0, 0) {
                    mesh.quad(
                        [[x1, y0, z0], [x1, y1, z0], [x1, y1, z1], [x1, y0, z1]],
                        [1.0, 0.0, 0.0],
                        c,
                    );
                }
            }
        }
    }
    stride
}

/// The radius to draw a body at, in metres.
///
/// **Not in the data.** A body set records where each body is and one value to colour it by; it
/// carries no extent, because nothing in the physics needs one — a point mass has no radius and an
/// atom's is a convention. So this is a quarter of the median distance to the nearest neighbour,
/// which keeps spheres apart at any scale from an orbit to a lattice, and [`Exported::notes`] says
/// it is a drawing convention.
///
/// The median rather than the minimum: one close pair in a hundred bodies would otherwise shrink
/// every sphere to nothing.
fn body_radius(positions: &[[f64; 3]], bounds: &[f64; 6]) -> f64 {
    let span = (bounds[3] - bounds[0])
        .abs()
        .max((bounds[4] - bounds[1]).abs())
        .max((bounds[5] - bounds[2]).abs());
    let fallback = if span > 0.0 { span * 0.02 } else { 1.0 };
    if positions.len() < 2 {
        return fallback;
    }
    let mut nearest: Vec<f64> = Vec::with_capacity(positions.len());
    for (i, a) in positions.iter().enumerate() {
        let mut best = f64::MAX;
        for (j, b) in positions.iter().enumerate() {
            if i == j {
                continue;
            }
            let d = ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt();
            if d < best {
                best = d;
            }
        }
        if best.is_finite() {
            nearest.push(best);
        }
    }
    if nearest.is_empty() {
        return fallback;
    }
    nearest.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = nearest[nearest.len() / 2];
    if median > 0.0 {
        median * 0.25
    } else {
        fallback
    }
}

/// Rings and sectors of the sphere a body becomes. 8 by 12 is 96 quads, which reads as round at
/// the size a body is drawn and costs a fifth of what 16 by 24 would.
const SPHERE: (usize, usize) = (8, 12);

/// Append one smooth-shaded sphere.
///
/// Smooth, unlike a voxel face: a sphere *is* round, so the normal is the outward direction at the
/// vertex and the facets are a budget rather than a shape.
fn sphere(mesh: &mut Mesh, centre: [f64; 3], radius: f64, c: [f32; 4]) {
    let (rings, sectors) = SPHERE;
    let base = mesh.positions.len() as u32;
    for r in 0..=rings {
        let phi = std::f64::consts::PI * r as f64 / rings as f64;
        for s in 0..=sectors {
            let theta = std::f64::consts::TAU * s as f64 / sectors as f64;
            let n = [
                (phi.sin() * theta.cos()) as f32,
                phi.cos() as f32,
                (phi.sin() * theta.sin()) as f32,
            ];
            mesh.positions.push([
                (centre[0] + radius * n[0] as f64) as f32,
                (centre[1] + radius * n[1] as f64) as f32,
                (centre[2] + radius * n[2] as f64) as f32,
            ]);
            mesh.normals.push(n);
            mesh.colours.push(c);
        }
    }
    let row = (sectors + 1) as u32;
    for r in 0..rings as u32 {
        for s in 0..sectors as u32 {
            let a = base + r * row + s;
            let b = a + row;
            mesh.indices.extend([a, b, b + 1, a, b + 1, a + 1]);
        }
    }
}

/// The value range of a panel, and whether it straddles zero.
///
/// A panel holding one value gets a range one wide so nothing divides by zero; the colour that
/// results is the bottom of the scale, which is honest for a field that is not varying.
fn span(values: &[f64]) -> (f64, f64, bool) {
    let lo = values.iter().copied().fold(f64::MAX, f64::min);
    let hi = values.iter().copied().fold(f64::MIN, f64::max);
    if !(lo.is_finite() && hi.is_finite()) {
        return (0.0, 1.0, false);
    }
    if hi > lo {
        (lo, hi, ramp::is_signed(lo, hi))
    } else {
        (lo, lo + 1.0, false)
    }
}

/// Assemble the document.
///
/// The binary layout is positions, normals, colours, then indices, per mesh, each padded to a
/// four-byte boundary — glTF requires an accessor's offset to be a multiple of its component
/// size, and every component here is four bytes wide.
fn document(title: &str, meshes: &[Mesh]) -> String {
    let mut blob: Vec<u8> = Vec::new();
    let mut views = Vec::new();
    let mut accessors = Vec::new();
    let mut nodes = Vec::new();
    let mut mesh_json = Vec::new();

    for mesh in meshes {
        let pos_view = push_view(&mut blob, &mut views, floats3(&mesh.positions), 34962);
        let (min, max) = extent(&mesh.positions);
        let pos_accessor = accessors.len();
        accessors.push(format!(
            "{{\"bufferView\":{pos_view},\"componentType\":5126,\"count\":{},\"type\":\"VEC3\",\
             \"min\":[{:.9},{:.9},{:.9}],\"max\":[{:.9},{:.9},{:.9}]}}",
            mesh.positions.len(),
            min[0],
            min[1],
            min[2],
            max[0],
            max[1],
            max[2]
        ));

        let mut attributes = format!("\"POSITION\":{pos_accessor}");

        if !mesh.normals.is_empty() {
            let view = push_view(&mut blob, &mut views, floats3(&mesh.normals), 34962);
            let k = accessors.len();
            accessors.push(format!(
                "{{\"bufferView\":{view},\"componentType\":5126,\"count\":{},\"type\":\"VEC3\"}}",
                mesh.normals.len()
            ));
            attributes.push_str(&format!(",\"NORMAL\":{k}"));
        }

        let col_view = push_view(&mut blob, &mut views, floats4(&mesh.colours), 34962);
        let col_accessor = accessors.len();
        accessors.push(format!(
            "{{\"bufferView\":{col_view},\"componentType\":5126,\"count\":{},\"type\":\"VEC4\"}}",
            mesh.colours.len()
        ));
        attributes.push_str(&format!(",\"COLOR_0\":{col_accessor}"));

        let indices = if mesh.indices.is_empty() {
            String::new()
        } else {
            let view = push_view(&mut blob, &mut views, u32s(&mesh.indices), 34963);
            let k = accessors.len();
            accessors.push(format!(
                "{{\"bufferView\":{view},\"componentType\":5125,\"count\":{},\"type\":\"SCALAR\"}}",
                mesh.indices.len()
            ));
            format!(",\"indices\":{k}")
        };

        mesh_json.push(format!(
            "{{\"name\":{},\"primitives\":[{{\"attributes\":{{{attributes}}}{indices},\
             \"mode\":{},\"material\":0}}]}}",
            quote(&mesh.name),
            mesh.mode
        ));
        nodes.push(format!(
            "{{\"name\":{},\"mesh\":{}}}",
            quote(&mesh.name),
            mesh_json.len() - 1
        ));
    }

    let node_list: Vec<String> = (0..nodes.len()).map(|i| i.to_string()).collect();
    format!(
        "{{\n\
         \"asset\":{{\"version\":\"2.0\",\"generator\":{}}},\n\
         \"scene\":0,\n\
         \"scenes\":[{{\"name\":{},\"nodes\":[{}]}}],\n\
         \"nodes\":[{}],\n\
         \"meshes\":[{}],\n\
         \"materials\":[{{\"name\":\"vertex colour\",\"pbrMetallicRoughness\":\
         {{\"baseColorFactor\":[1,1,1,1],\"metallicFactor\":0,\"roughnessFactor\":0.55}},\
         \"doubleSided\":true}}],\n\
         \"buffers\":[{{\"byteLength\":{},\"uri\":\"data:application/octet-stream;base64,{}\"}}],\n\
         \"bufferViews\":[{}],\n\
         \"accessors\":[{}]\n\
         }}\n",
        quote(&format!("pantometry-view {}", env!("CARGO_PKG_VERSION"))),
        quote(title),
        node_list.join(","),
        nodes.join(","),
        mesh_json.join(","),
        blob.len(),
        base64(&blob),
        views.join(","),
        accessors.join(",")
    )
}

/// Append bytes to the blob, pad to four, and record a view over them.
fn push_view(blob: &mut Vec<u8>, views: &mut Vec<String>, bytes: Vec<u8>, target: u32) -> usize {
    while blob.len() % 4 != 0 {
        blob.push(0);
    }
    let offset = blob.len();
    let length = bytes.len();
    blob.extend(bytes);
    views.push(format!(
        "{{\"buffer\":0,\"byteOffset\":{offset},\"byteLength\":{length},\"target\":{target}}}"
    ));
    views.len() - 1
}

fn floats3(v: &[[f32; 3]]) -> Vec<u8> {
    v.iter()
        .flat_map(|p| p.iter().flat_map(|f| f.to_le_bytes()))
        .collect()
}

fn floats4(v: &[[f32; 4]]) -> Vec<u8> {
    v.iter()
        .flat_map(|p| p.iter().flat_map(|f| f.to_le_bytes()))
        .collect()
}

fn u32s(v: &[u32]) -> Vec<u8> {
    v.iter().flat_map(|i| i.to_le_bytes()).collect()
}

/// The tightest box round a set of positions. **Required** by the spec on a `POSITION` accessor,
/// and a viewer that cannot frame a scene is the usual symptom of leaving it out.
fn extent(positions: &[[f32; 3]]) -> ([f32; 3], [f32; 3]) {
    let mut min = [f32::MAX; 3];
    let mut max = [f32::MIN; 3];
    for p in positions {
        for a in 0..3 {
            min[a] = min[a].min(p[a]);
            max[a] = max[a].max(p[a]);
        }
    }
    if positions.is_empty() {
        return ([0.0; 3], [0.0; 3]);
    }
    (min, max)
}

/// Standard base64, written out rather than depended on. Twenty lines against a crate.
fn base64(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(bytes.len().div_ceil(3) * 4);
    for chunk in bytes.chunks(3) {
        let b = [
            chunk[0],
            chunk.get(1).copied().unwrap_or(0),
            chunk.get(2).copied().unwrap_or(0),
        ];
        let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
        out.push(ALPHABET[(n >> 18) as usize & 63] as char);
        out.push(ALPHABET[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 {
            ALPHABET[(n >> 6) as usize & 63] as char
        } else {
            '='
        });
        out.push(if chunk.len() > 2 {
            ALPHABET[n as usize & 63] as char
        } else {
            '='
        });
    }
    out
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            c if c.is_control() => out.push(' '),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// One panel's geometry, for a caller that wants a single domain rather than the frame.
///
/// The same document, with one node in it.
pub fn gltf_panel(title: &str, panel: &Panel) -> Exported {
    gltf(
        title,
        &Frame {
            time_s: 0.0,
            panels: vec![panel.clone()],
            readings: Vec::new(),
        },
    )
}
