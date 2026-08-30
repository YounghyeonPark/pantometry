//! Turning a panel into triangles, once, for every exporter that needs them — and carrying a
//! designed mesh across when the triangles were somebody else's to begin with.
//!
//! # Why this is its own module
//!
//! Two writers want the same geometry and want it in different files: `gltf` writes one frame as
//! glTF and `usd` writes a whole run as USD. Both need the surface of a field and a sphere for a
//! body, and neither should have its own opinion about where a cell's faces are — a solid exported
//! one size from one file and another size from the other is worse than either being wrong,
//! because the two disagree and nothing says which to believe.
//!
//! It is also the piece that makes a **run** exportable. glTF cannot animate a field, so `gltf`
//! takes one frame; USD can, and a run's colours change while its topology does not. Separating
//! the two is what lets USD write the topology once and the colours per frame, and it is why
//! [`Surface::source`] exists at all.
//!
//! # A field is the boundary of its material
//!
//! One quad — two triangles — per cell face whose neighbour is absent or outside the grid. A solid
//! 9×9×9 block is 486 quads and not 4,374, because 89% of the faces of a solid are between two
//! present cells and no camera can see them. A void inside a part produces a real interior
//! surface, which is the same rule read from the other side.
//!
//! Flat-shaded, four vertices to a quad, each carrying that face's normal. A voxel surface *is*
//! faceted; sharing corner vertices would average three perpendicular normals into a direction no
//! face points and shade a cube as a blob.
//!
//! # And one shape that is not derived from anything
//!
//! Three of the four producers here *compute* a surface — from the boundary of a field's cells,
//! from a level through its values, from a body's position and a drawing convention.
//! [`mesh_surface`] does not: a designed part arrived as triangles and the work is to carry them
//! across unchanged, flat-shaded, with the winding the file wrote. It is the only one whose input
//! is not a panel, which is why the module is described twice above.
//!
//! **Cells are clipped to the extent.** [`capture`](pantometry_scene::capture) samples corner to
//! corner, so the first and last samples sit *on* the faces of the box rather than half a cell
//! inside it. Giving every sample a full cell either side made a 40 mm cube come out 80 mm. An end
//! node owns half a cell — the same arithmetic whose absence made every room mode in this
//! workspace read low, when `Tube` and `Room` divided a wall sample's divergence by a whole `dx`.

use pantometry_scene::PanelData;

/// The most faces one field may write before it is subsampled.
///
/// A 100³ field has a 60,000-face surface, which is fine. A 100³ field *full of holes* can have
/// twenty times that, which is what this is really for: past the budget the field is walked with a
/// stride and the caller is told, because a silently decimated surface is a picture of a coarser
/// simulation than the one that ran.
pub const MAX_FACES: usize = 200_000;

/// Rings and sectors of the sphere a body becomes.
///
/// 8 by 12 reads as round at the size a body is drawn and costs a fifth of what 16 by 24 would.
/// The seam column and the pole rows are duplicated so each vertex can carry its own normal, so
/// the count is `(rings + 1) * (sectors + 1)` — 117, not 96.
pub const SPHERE: (usize, usize) = (8, 12);

/// Triangles, with a note beside each vertex of where its value came from.
#[derive(Clone, Debug, Default)]
pub struct Surface {
    /// Vertex positions, in metres.
    pub positions: Vec<[f32; 3]>,
    /// One normal per vertex. Empty only for geometry that takes no light.
    pub normals: Vec<[f32; 3]>,
    /// Triangle indices, three to a face.
    pub indices: Vec<u32>,
    /// For each **vertex**, the index in the panel's `values` that colours it.
    ///
    /// This is what makes the topology reusable across frames: a run's colours change while its
    /// cells do not, so a writer that can animate — USD — writes the positions once and looks the
    /// colours up here per frame. A writer that cannot — glTF — uses it once.
    pub source: Vec<u32>,
    /// How many cells were skipped per axis. One means every cell is drawn.
    pub stride: usize,
}

impl Surface {
    /// How many triangles.
    pub fn faces(&self) -> usize {
        self.indices.len() / 3
    }

    /// Add one flat-shaded quad: four corners in winding order, one normal, one source index.
    fn quad(&mut self, corners: [[f32; 3]; 4], normal: [f32; 3], source: u32) {
        let base = self.positions.len() as u32;
        for c in corners {
            self.positions.push(c);
            self.normals.push(normal);
            self.source.push(source);
        }
        self.indices
            .extend([base, base + 1, base + 2, base, base + 2, base + 3]);
    }

    /// The tightest box round the positions, as `[x0, y0, z0, x1, y1, z1]`.
    ///
    /// Both formats require it: glTF on a `POSITION` accessor, USD as a prim's `extent`. A viewer
    /// that cannot get a bounding box out of a file frames the scene at the origin, which looks
    /// exactly like the geometry being in the wrong place.
    pub fn bounds(&self) -> [f32; 6] {
        let mut b = [f32::MAX, f32::MAX, f32::MAX, f32::MIN, f32::MIN, f32::MIN];
        if self.positions.is_empty() {
            return [0.0; 6];
        }
        for p in &self.positions {
            for a in 0..3 {
                b[a] = b[a].min(p[a]);
                b[a + 3] = b[a + 3].max(p[a]);
            }
        }
        b
    }
}

/// Which surface an exporter draws for a field.
///
/// The two answer different questions and neither is a better version of the other.
/// [`field_surface`] draws the outside of the cells that hold a value: where the material is,
/// blocky and correct, and a picture of the grid as much as of the object. [`isosurface`] draws
/// where the values reach a number: where it is 100 degrees, where the pressure is zero, where
/// the melt front sits.
///
/// A **selection** rather than a second exporter, because everything after the mesh is
/// identical — the colours, the bounds, the accessors, the winding — and a second copy of that
/// is a second chance at the defects it has already had.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub enum Surfaces {
    /// The boundary of the cells that hold a value. What every export wrote before there was a
    /// choice, and so the default: a file that changed shape because a new variant appeared
    /// would be worse than one that never gained the option.
    #[default]
    Boundary,
    /// Where the field reaches this value.
    ///
    /// A level outside the field's range produces **nothing**, which the exporters report as a
    /// skipped panel rather than as an empty file. That is the same silence
    /// `pantometry-app`'s isosurface control had to learn to label: an absence a reader can see
    /// beats one they cannot.
    At(f64),
}

impl Surfaces {
    /// Build the surface this selects, from a field's grid and values.
    ///
    /// One place, so `gltf` and `usd` cannot come to different answers about the same panel —
    /// which is the argument this module's header already makes about the two writers.
    pub fn of(self, counts: (usize, usize, usize), extent: [f64; 6], values: &[f64]) -> Surface {
        match self {
            Surfaces::Boundary => field_surface(counts, extent, values),
            Surfaces::At(level) => isosurface(counts, extent, values, level),
        }
    }
}

/// The surface inside a field where it equals `level`.
///
/// The first thing in this workspace with a **surface** rather than a boundary. `field_surface`
/// draws the outside of the cells that hold a value — a staircase, correct and blocky, and a
/// picture of the grid as much as of the object. This draws where the *values* reach a number,
/// which is the question a person actually asks of a field: where is it 100 °C, where is the
/// pressure zero, where does the melt front sit.
///
/// # Marching tetrahedra, and why not cubes
///
/// Marching cubes needs a 256-entry case table, and several of those entries are ambiguous — the
/// same corner signs admit two different surfaces, and choosing wrong leaves a hole. This splits
/// every cell into **six tetrahedra** instead: four corners, sixteen sign patterns, and only two
/// shapes among them. There is no table to get wrong and no ambiguous case to resolve, and the
/// result is watertight by construction because adjacent cells split their shared face the same
/// way.
///
/// The cost is more triangles for the same surface. That is the right trade here: a large table is
/// a large surface to be quietly wrong on, and a hole in a closed shape is the kind of wrong that
/// looks like a rendering artefact rather than like a bug.
///
/// # A void is not a value
///
/// A cell with a non-finite sample is skipped entirely, along with every tetrahedron touching it.
/// `Solid3D::temperature_at` returns `NaN` for an emptied cell deliberately — "a zero or an
/// ambient is a value somebody would plot and believe" — and interpolating an edge that ends in
/// one would invent the very number that `NaN` exists to refuse.
///
/// # What comes back
///
/// Vertices are shared between triangles by the edge they sit on, so the mesh is watertight and
/// each vertex gets a normal averaged from the faces around it. [`Surface::source`] points at the
/// nearer of the two samples the vertex sits between, which is what colours it.
pub fn isosurface(
    counts: (usize, usize, usize),
    extent: [f64; 6],
    values: &[f64],
    level: f64,
) -> Surface {
    let (nx, ny, nz) = counts;
    let mut out = Surface {
        stride: 1,
        ..Surface::default()
    };
    // Two samples along every axis is the minimum that has an inside: a plane of samples has no
    // cell to march through.
    if nx < 2 || ny < 2 || nz < 2 || values.len() < nx * ny * nz || !level.is_finite() {
        return out;
    }

    let at = |i: usize, j: usize, k: usize| values[i + nx * (j + ny * k)];
    // Samples sit at cell centres of a grid spanning `extent`, corner to corner — the same
    // convention `field_surface` uses, and the one `mesh` had to be corrected to once when a
    // 40 mm cube exported 80 mm across.
    let step = |lo: f64, hi: f64, n: usize| {
        if n > 1 {
            (hi - lo) / (n - 1) as f64
        } else {
            0.0
        }
    };
    let (dx, dy, dz) = (
        step(extent[0], extent[3], nx),
        step(extent[1], extent[4], ny),
        step(extent[2], extent[5], nz),
    );
    let position = |i: usize, j: usize, k: usize| {
        [
            extent[0] + dx * i as f64,
            extent[1] + dy * j as f64,
            extent[2] + dz * k as f64,
        ]
    };

    // The six tetrahedra of a cube, as indices into its eight corners in the binary order bit 0 =
    // x, bit 1 = y, bit 2 = z. Every one of them contains the 0–7 diagonal, which is what makes
    // neighbouring cells agree about how their shared face is cut.
    const TETS: [[usize; 4]; 6] = [
        [0, 7, 1, 3],
        [0, 7, 3, 2],
        [0, 7, 2, 6],
        [0, 7, 6, 4],
        [0, 7, 4, 5],
        [0, 7, 5, 1],
    ];

    // One vertex per cut edge, keyed by the two samples it lies between, so triangles share it.
    let mut vertices: std::collections::HashMap<(usize, usize), u32> =
        std::collections::HashMap::new();

    for k in 0..nz - 1 {
        for j in 0..ny - 1 {
            for i in 0..nx - 1 {
                let corner = |c: usize| (i + (c & 1), j + ((c >> 1) & 1), k + ((c >> 2) & 1));
                let sample = |c: usize| {
                    let (a, b, d) = corner(c);
                    (a + nx * (b + ny * d), at(a, b, d))
                };
                let corners: [(usize, f64); 8] = std::array::from_fn(sample);
                // A cell touching a void is not marched. Skipping the whole cell rather than the
                // tetrahedra that touch the void keeps the surface closed: a partial cell would
                // leave an edge with nothing on the other side of it.
                if corners.iter().any(|(_, v)| !v.is_finite()) {
                    continue;
                }

                for tet in TETS {
                    let inside: [bool; 4] = std::array::from_fn(|n| corners[tet[n]].1 < level);
                    let count = inside.iter().filter(|x| **x).count();
                    if count == 0 || count == 4 {
                        continue;
                    }
                    // The cut edges are the ones whose ends disagree. One corner apart from the
                    // other three gives a triangle; two against two gives a quad, as two.
                    let mut cut: Vec<(usize, usize)> = Vec::with_capacity(4);
                    for a in 0..4 {
                        for b in a + 1..4 {
                            if inside[a] != inside[b] {
                                cut.push((tet[a], tet[b]));
                            }
                        }
                    }
                    // Behind the surface: the average of the corners that are below `level`.
                    // Every face of this tetrahedron must point away from it.
                    let mut behind = [0.0f64; 3];
                    let mut n_behind = 0.0f64;
                    for (n, c) in tet.iter().enumerate() {
                        if inside[n] {
                            let (a, b, d) = corner(*c);
                            let p = position(a, b, d);
                            for x in 0..3 {
                                behind[x] += p[x];
                            }
                            n_behind += 1.0;
                        }
                    }
                    for x in behind.iter_mut() {
                        *x /= n_behind.max(1.0);
                    }

                    let mut idx: Vec<u32> = Vec::with_capacity(cut.len());
                    for (a, b) in &cut {
                        let (ia, va) = corners[*a];
                        let (ib, vb) = corners[*b];
                        let key = (ia.min(ib), ia.max(ib));
                        let handle = *vertices.entry(key).or_insert_with(|| {
                            // Where the straight line between the two samples crosses `level`.
                            // Guarded: two equal samples that straddle nothing would divide by
                            // zero, and the midpoint is the only answer that is not a direction.
                            let t = if (vb - va).abs() > f64::MIN_POSITIVE {
                                ((level - va) / (vb - va)).clamp(0.0, 1.0)
                            } else {
                                0.5
                            };
                            let (pa, pb) = (
                                position(corner(*a).0, corner(*a).1, corner(*a).2),
                                position(corner(*b).0, corner(*b).1, corner(*b).2),
                            );
                            out.positions.push([
                                (pa[0] + (pb[0] - pa[0]) * t) as f32,
                                (pa[1] + (pb[1] - pa[1]) * t) as f32,
                                (pa[2] + (pb[2] - pa[2]) * t) as f32,
                            ]);
                            out.normals.push([0.0; 3]);
                            out.source.push(if t < 0.5 { ia as u32 } else { ib as u32 });
                            (out.positions.len() - 1) as u32
                        });
                        idx.push(handle);
                    }
                    // Three cut edges is one triangle; four is a quad, and the two triangles must
                    // share a diagonal of it rather than be picked independently.
                    match idx.len() {
                        3 => emit(&mut out, [idx[0], idx[1], idx[2]], behind),
                        4 => {
                            emit(&mut out, [idx[0], idx[1], idx[3]], behind);
                            emit(&mut out, [idx[0], idx[3], idx[2]], behind);
                        }
                        _ => {}
                    }
                }
            }
        }
    }

    // Vertex normals, accumulated from the faces that meet there. A face's own normal is its
    // cross product; summing and normalising is the standard average weighted by area, which is
    // what makes a curved surface look curved rather than faceted.
    for f in out.indices.chunks_exact(3) {
        let p = |n: u32| out.positions[n as usize];
        let (a, b, c) = (p(f[0]), p(f[1]), p(f[2]));
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        for h in f {
            let slot = &mut out.normals[*h as usize];
            for a in 0..3 {
                slot[a] += n[a];
            }
        }
    }
    for n in out.normals.iter_mut() {
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if len > f32::MIN_POSITIVE {
            for a in n.iter_mut() {
                *a /= len;
            }
        }
    }
    out
}

/// Add one triangle, wound so its face points **away from the low side** of the field.
///
/// Consistent winding is not decoration: it is what lets a caller compute an enclosed volume from
/// the mesh at all, and what stops a renderer with back-face culling showing a surface only from
/// inside. The reference is the tetrahedron's own geometry — the centroid of the corners below
/// `level` is behind the face, so a normal pointing away from it points out.
fn emit(out: &mut Surface, tri: [u32; 3], behind: [f64; 3]) {
    let p = |n: u32| out.positions[n as usize];
    let (a, b, c) = (p(tri[0]), p(tri[1]), p(tri[2]));
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    // Positive when the face already points away from what is behind it.
    let outward = (0..3)
        .map(|x| n[x] as f64 * (a[x] as f64 - behind[x]))
        .sum::<f64>();
    if outward >= 0.0 {
        out.indices.extend_from_slice(&tri);
    } else {
        out.indices.extend_from_slice(&[tri[0], tri[2], tri[1]]);
    }
}

/// The surface of the cells of a field that hold a value.
///
/// A field of one dimension returns an empty surface: a row of samples along a line is a graph and
/// not geometry, and the caller reports that rather than drawing dots in space.
pub fn field_surface(counts: (usize, usize, usize), extent: [f64; 6], values: &[f64]) -> Surface {
    let (nx, ny, nz) = counts;
    let mut out = Surface {
        stride: 1,
        ..Surface::default()
    };
    if nx == 0 || ny == 0 || nz == 0 || values.len() < nx * ny * nz {
        return out;
    }
    if [nx, ny, nz].iter().filter(|&&n| n > 1).count() < 2 {
        return out;
    }

    // Stride so the face count stays inside the budget. The surface goes as the square of the
    // linear resolution, so the square root of the overrun is the right reduction; it only has to
    // be in the right decade to keep the file openable.
    let faces = 2 * (nx * ny + ny * nz + nz * nx);
    out.stride = if faces > MAX_FACES {
        ((faces as f64 / MAX_FACES as f64).sqrt().ceil() as usize).max(1)
    } else {
        1
    };
    let stride = out.stride;

    let at = |i: usize, j: usize, k: usize| -> Option<f64> {
        let v = *values.get(i + nx * (j + ny * k))?;
        v.is_finite().then_some(v)
    };
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
    // Half a cell either side, which is what makes a cell a box rather than a point. A collapsed
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
    let clip = |lo: f64, hi: f64, a: usize| -> (f64, f64) {
        if size[a].abs() > 0.0 {
            let (elo, ehi) = (extent[a].min(extent[a + 3]), extent[a].max(extent[a + 3]));
            (lo.max(elo), hi.min(ehi))
        } else {
            (lo, hi)
        }
    };

    let taken = |n: usize| n.div_ceil(stride);
    for kt in 0..taken(nz) {
        for jt in 0..taken(ny) {
            for it in 0..taken(nx) {
                let (i, j, k) = (it * stride, jt * stride, kt * stride);
                if at(i, j, k).is_none() {
                    continue;
                }
                let src = (i + nx * (j + ny * k)) as u32;
                let m = [
                    extent[0] + size[0] * along(i, nx),
                    extent[1] + size[1] * along(j, ny),
                    extent[2] + size[2] * along(k, nz),
                ];
                let (a0, a1) = clip(m[0] - hx, m[0] + hx, 0);
                let (b0, b1) = clip(m[1] - hy, m[1] + hy, 1);
                let (c0, c1) = clip(m[2] - hz, m[2] + hz, 2);
                let (x0, x1) = (a0 as f32, a1 as f32);
                let (y0, y1) = (b0 as f32, b1 as f32);
                let (z0, z1) = (c0 as f32, c1 as f32);

                // A neighbour one stride away, in the strided lattice: the cells between are not
                // drawn, so a face is exposed when the **drawn** neighbour is missing.
                let exposed = |di: isize, dj: isize, dk: isize| -> bool {
                    let s = stride as isize;
                    let (ni, nj, nk) = (
                        i as isize + di * s,
                        j as isize + dj * s,
                        k as isize + dk * s,
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
                // Counter-clockwise seen from outside, which is what glTF's front face is and what
                // USD's default `rightHanded` orientation means.
                if exposed(0, 0, -1) {
                    out.quad(
                        [[x0, y0, z0], [x0, y1, z0], [x1, y1, z0], [x1, y0, z0]],
                        [0.0, 0.0, -1.0],
                        src,
                    );
                }
                if exposed(0, 0, 1) {
                    out.quad(
                        [[x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]],
                        [0.0, 0.0, 1.0],
                        src,
                    );
                }
                if exposed(0, -1, 0) {
                    out.quad(
                        [[x0, y0, z0], [x1, y0, z0], [x1, y0, z1], [x0, y0, z1]],
                        [0.0, -1.0, 0.0],
                        src,
                    );
                }
                if exposed(0, 1, 0) {
                    out.quad(
                        [[x0, y1, z0], [x0, y1, z1], [x1, y1, z1], [x1, y1, z0]],
                        [0.0, 1.0, 0.0],
                        src,
                    );
                }
                if exposed(-1, 0, 0) {
                    out.quad(
                        [[x0, y0, z0], [x0, y0, z1], [x0, y1, z1], [x0, y1, z0]],
                        [-1.0, 0.0, 0.0],
                        src,
                    );
                }
                if exposed(1, 0, 0) {
                    out.quad(
                        [[x1, y0, z0], [x1, y1, z0], [x1, y1, z1], [x1, y0, z1]],
                        [1.0, 0.0, 0.0],
                        src,
                    );
                }
            }
        }
    }
    out
}

/// The surface a designed mesh already is.
///
/// The other three producers in this module *derive* a surface: from the boundary of a field's
/// cells, from a level through its values, from a body's position and a drawing convention. A
/// designed mesh needs none of that — it arrived as triangles and the work is to carry them
/// across without changing them.
///
/// # Triangles, not a `Mesh`
///
/// `pantometry-shape` owns meshes and this crate does not depend on it, deliberately: layer 3 may
/// read layer 2's scene and nothing below it. So this takes the same plain arrays every other
/// producer here takes, and the caller — which links both — does the conversion. It costs one
/// `collect` at the boundary and keeps a rule that has held through eleven domains.
///
/// # Flat, and one `source`
///
/// Per-face normals, three vertices to a triangle, nothing shared. The argument is the one this
/// module's header makes about voxel faces and it survives the change of subject: a designed mesh
/// is faceted *by construction* — an STL has no curvature to preserve, only the tessellation of
/// one — and averaging normals across a machined edge would round a chamfer that is not there.
/// Smoothing here would draw a shape the file does not describe.
///
/// Every vertex carries the caller's `source`, one for the whole mesh, because a part is one
/// object with one material. Nothing in an STL varies along it, so a per-vertex index would be
/// three copies of the same number and an invitation to colour a solid by an accident.
///
/// # What is dropped, and what is not budgeted
///
/// A triangle whose normal is not finite and positive is skipped: zero area means no normal
/// exists, and a non-finite vertex poisons the cross product to `NaN`, so the one test catches
/// both. Exporters emit degenerate triangles routinely and a vertex normal of `NaN` shades a face
/// black. `a_designed_mesh.rs` pins that both go, and that a mesh of nothing but degenerates comes
/// out **empty** rather than malformed.
///
/// Unlike [`field_surface`] there is no [`MAX_FACES`] budget, because subsampling has no meaning
/// here. Dropping every other cell of a field draws a coarser field; dropping every other triangle
/// of a mesh punches holes in a solid. A million-triangle STL is the caller's to decimate, with a
/// tool that knows what the surface is.
pub fn mesh_surface(triangles: &[[[f64; 3]; 3]], source: u32) -> Surface {
    let mut out = Surface {
        stride: 1,
        ..Surface::default()
    };
    for t in triangles {
        let u = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
        let v = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
        // The winding is the file's. `a` to `b` to `c` counter-clockwise seen from outside is
        // what every STL writer means, and `(b - a) x (c - a)` is the direction that follows
        // from it -- so a mesh that was inside-out on disk is drawn inside-out here, which is
        // the honest outcome. The cube test asserts a positive enclosed volume to prove the
        // conversion does not flip it.
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        if !len.is_finite() || len <= 0.0 {
            continue;
        }
        let unit = [
            (n[0] / len) as f32,
            (n[1] / len) as f32,
            (n[2] / len) as f32,
        ];
        let base = out.positions.len() as u32;
        for p in t {
            out.positions.push([p[0] as f32, p[1] as f32, p[2] as f32]);
            out.normals.push(unit);
            out.source.push(source);
        }
        out.indices.extend([base, base + 1, base + 2]);
    }
    out
}

/// A sphere for each body.
///
/// Smooth-shaded, unlike a voxel face: a sphere *is* round, so the normal is the outward direction
/// at the vertex and the facets are a budget rather than a shape.
pub fn body_spheres(positions: &[[f64; 3]], radius: f64) -> Surface {
    let (rings, sectors) = SPHERE;
    let mut out = Surface {
        stride: 1,
        ..Surface::default()
    };
    for (b, centre) in positions.iter().enumerate() {
        let base = out.positions.len() as u32;
        for r in 0..=rings {
            let phi = std::f64::consts::PI * r as f64 / rings as f64;
            for s in 0..=sectors {
                let theta = std::f64::consts::TAU * s as f64 / sectors as f64;
                let n = [
                    (phi.sin() * theta.cos()) as f32,
                    phi.cos() as f32,
                    (phi.sin() * theta.sin()) as f32,
                ];
                out.positions.push([
                    (centre[0] + radius * n[0] as f64) as f32,
                    (centre[1] + radius * n[1] as f64) as f32,
                    (centre[2] + radius * n[2] as f64) as f32,
                ]);
                out.normals.push(n);
                out.source.push(b as u32);
            }
        }
        let row = (sectors + 1) as u32;
        for r in 0..rings as u32 {
            for s in 0..sectors as u32 {
                let a = base + r * row + s;
                let c = a + row;
                out.indices.extend([a, c, c + 1, a, c + 1, a + 1]);
            }
        }
    }
    out
}

/// The radius to draw a body at, in metres.
///
/// **Not in the data.** A body set records where each body is and one value to colour it by; it
/// carries no extent, because nothing in the physics needs one — a point mass has no radius and an
/// atom's is a convention. So this is a quarter of the median distance to the nearest neighbour,
/// which keeps spheres apart at any scale from an orbit to a lattice, and every caller states that
/// it is a drawing convention.
///
/// The median rather than the minimum: one close pair among a hundred bodies would otherwise
/// shrink every sphere to nothing.
pub fn body_radius(positions: &[[f64; 3]], bounds: &[f64; 6]) -> f64 {
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

/// The value range of a panel's samples, and whether it straddles zero.
///
/// A panel holding one value gets a range one wide so nothing divides by zero.
pub fn span(values: &[f64]) -> (f64, f64, bool) {
    let lo = values.iter().copied().fold(f64::MAX, f64::min);
    let hi = values.iter().copied().fold(f64::MIN, f64::max);
    if !(lo.is_finite() && hi.is_finite()) {
        return (0.0, 1.0, false);
    }
    if hi > lo {
        (lo, hi, crate::ramp::is_signed(lo, hi))
    } else {
        (lo, lo + 1.0, false)
    }
}

/// Where a value sits on its colour scale, in `0..=1`.
///
/// A signed range is placed about **zero** rather than about the middle of the range, so the
/// neutral colour lands on the value zero. For −100 to +300 those are a quarter of the scale apart.
pub fn place(v: f64, lo: f64, hi: f64, signed: bool) -> f64 {
    if signed {
        let reach = lo.abs().max(hi.abs()).max(f64::MIN_POSITIVE);
        0.5 + 0.5 * v / reach
    } else {
        (v - lo) / (hi - lo)
    }
}

/// Whether every frame agrees about which cells of a field hold a value.
///
/// The question a writer that animates has to ask before it writes topology once and colours many
/// times. Voids in this workspace are set when a scene is built and do not move, so the answer is
/// normally yes — but "normally" is not a thing to assume about a file somebody else will open, and
/// the alternative is geometry from one frame presented as the run's.
pub fn presence_is_constant(frames: &[&PanelData]) -> bool {
    let mut first: Option<Vec<bool>> = None;
    for data in frames {
        let PanelData::Field { values, .. } = data else {
            continue;
        };
        let here: Vec<bool> = values.iter().map(|v| v.is_finite()).collect();
        match &first {
            None => first = Some(here),
            Some(f) => {
                if *f != here {
                    return false;
                }
            }
        }
    }
    true
}
