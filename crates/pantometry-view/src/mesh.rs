//! Turning a panel into triangles, once, for every exporter that needs them.
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
