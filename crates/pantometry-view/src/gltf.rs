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

use crate::mesh::{self, Surface};
use crate::ramp;
use pantometry_scene::{Frame, Panel, PanelData};

/// The most bodies that become spheres before the exporter falls back to a point cloud.
///
/// A sphere is 96 quads, so 256 bodies is about 49,000 triangles — a comfortable Blender scene.
/// Ten thousand bodies would be two million, which is a worse picture than points and a much
/// slower one.
pub const MAX_SPHERES: usize = 256;

/// The most cell faces one field may write. See [`crate::mesh::MAX_FACES`], which is where the
/// budget and the striding live now.
pub const MAX_FACES: usize = crate::mesh::MAX_FACES;

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

    // There was a `quad` here, and building faces is [`crate::mesh`]'s job now — which is the
    // whole point of that module: two exporters cannot disagree about where a cell's faces are if
    // only one of them decides. `from_surface` below is all that is left of the geometry here.
}

/// Export one frame.
///
/// Every panel that is geometry becomes a node under one scene, named after its domain, so a
/// reader opening the file sees the same names the run reported.
pub fn gltf(title: &str, frame: &Frame) -> Exported {
    gltf_with(title, frame, mesh::Surfaces::Boundary)
}

/// The same, choosing which surface a field becomes. See [`mesh::Surfaces`].
///
/// Separate rather than an extra argument on [`gltf`], because that one is published and the
/// files it writes are the ones every caller here already expects. A default nobody asked for is
/// how an export quietly changes shape between versions.
pub fn gltf_with(title: &str, frame: &Frame, surfaces: mesh::Surfaces) -> Exported {
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
                let (lo, hi, signed) = mesh::span(values);
                coloured = true;
                let mut strands = Mesh::new(&panel.name, 1);
                strands.positions.reserve(vertices.len());
                for (k, start) in starts.iter().enumerate() {
                    let from = *start;
                    let to = starts.get(k + 1).copied().unwrap_or(vertices.len());
                    let c = colour(
                        mesh::place(values.get(k).copied().unwrap_or(lo), lo, hi, signed),
                        signed,
                    );
                    let base = strands.positions.len() as u32;
                    for v in &vertices[from..to] {
                        strands
                            .positions
                            .push([v[0] as f32, v[1] as f32, v[2] as f32]);
                        strands.colours.push(c);
                    }
                    // One segment per adjacent pair. Indices rather than duplicated vertices, so a
                    // ray of n points costs n positions and not 2(n-1).
                    for i in 0..(to - from).saturating_sub(1) as u32 {
                        strands.indices.push(base + i);
                        strands.indices.push(base + i + 1);
                    }
                }
                if !strands.indices.is_empty() {
                    meshes.push(strands);
                }
            }
            PanelData::Points {
                positions,
                values,
                bounds,
                ..
            } => {
                let (lo, hi, signed) = mesh::span(values);
                coloured = true;
                if positions.len() <= MAX_SPHERES {
                    let r = mesh::body_radius(positions, bounds);
                    notes.push(format!(
                        "{}: bodies drawn as spheres of radius {r:.6} m — a **size this run does \
                         not carry**. A body set records positions and a value, not an extent, so \
                         the radius is a quarter of the median distance to the nearest neighbour, \
                         which keeps them apart and is a drawing convention rather than a \
                         measurement",
                        panel.name
                    ));
                    let surface = mesh::body_spheres(positions, r);
                    meshes.push(from_surface(&panel.name, &surface, values, lo, hi, signed));
                } else {
                    notes.push(format!(
                        "{}: {} bodies is over the {MAX_SPHERES} a sphere each is worth, so this \
                         is a point cloud — which takes no light and casts no shadow",
                        panel.name,
                        positions.len()
                    ));
                    let mut m = Mesh::new(&panel.name, 0);
                    for (i, p) in positions.iter().enumerate() {
                        m.positions.push([p[0] as f32, p[1] as f32, p[2] as f32]);
                        m.colours.push(colour(
                            mesh::place(values.get(i).copied().unwrap_or(lo), lo, hi, signed),
                            signed,
                        ));
                    }
                    meshes.push(m);
                }
            }
            PanelData::Field {
                nx,
                ny,
                nz,
                extent_m,
                values,
            } => {
                if [*nx, *ny, *nz].iter().filter(|&&n| n > 1).count() < 2 {
                    skipped.push(format!(
                        "{} is a {nx}x{ny}x{nz} field: a row of samples along a line is a graph, \
                         not geometry, so it is not in the scene",
                        panel.name
                    ));
                    continue;
                }
                let (lo, hi, signed) = mesh::span(values);
                let surface = surfaces.of((*nx, *ny, *nz), *extent_m, values);
                if surface.indices.is_empty() {
                    // **Two different silences, and they used to share one sentence.** An empty
                    // boundary means no cell holds a value. An empty *level* means the field is
                    // full of values and none of them is the one asked for, which is a number a
                    // person can correct -- so it says what the field actually spans.
                    skipped.push(match surfaces {
                        mesh::Surfaces::Boundary => format!(
                            "{} is a {nx}x{ny}x{nz} field with no cell that holds a value, so it has no surface",
                            panel.name
                        ),
                        mesh::Surfaces::At(level) => {
                            let (lo, hi, _) = mesh::span(values);
                            format!(
                                "{} never reaches {level}: the field spans {lo} to {hi}, so there is no surface at that level",
                                panel.name
                            )
                        }
                    });
                    continue;
                }
                if surface.stride > 1 {
                    notes.push(format!(
                        "{}: subsampled every {} cells — the full surface is over {MAX_FACES} \
                         faces",
                        panel.name, surface.stride
                    ));
                }
                coloured = true;
                meshes.push(from_surface(&panel.name, &surface, values, lo, hi, signed));
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

/// A [`Surface`] as a glTF mesh, coloured from one frame's values.
///
/// The geometry is [`crate::mesh`]'s, shared with the USD writer, so a solid does not come out one
/// size from one file and another size from the other.
fn from_surface(
    name: &str,
    surface: &Surface,
    values: &[f64],
    lo: f64,
    hi: f64,
    signed: bool,
) -> Mesh {
    Mesh {
        name: name.to_string(),
        positions: surface.positions.clone(),
        normals: surface.normals.clone(),
        colours: surface
            .source
            .iter()
            .map(|src| {
                let v = values.get(*src as usize).copied().unwrap_or(f64::NAN);
                colour(mesh::place(v, lo, hi, signed), signed)
            })
            .collect(),
        indices: surface.indices.clone(),
        mode: 4,
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
