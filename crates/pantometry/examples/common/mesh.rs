//! Solids built from a path and from a point: a swept tube and a sphere.
//!
//! [`PanelData::Surface`](pantometry::scene::PanelData::Surface) takes triangles and says, in its
//! own documentation, that the winding is the caller's to get right and is not checked. This is
//! where the caller gets it right, and where the getting-it-right is measured: a mesh built here
//! is checked against a *closed form for the mesh itself* — the exact area and volume of a
//! regular prism, and of an icosahedron — rather than against a second implementation of the same
//! sweep. An inside-out mesh has negative volume and no picture shows it.
//!
//! What it does: a circular cross-section swept along a polyline with parallel-transport frames,
//! and a sphere by subdividing an icosahedron onto its circumsphere. Both come out **closed**:
//! every edge shared by exactly two triangles.
//! What it does not: booleans, unions, a solvent-excluded surface, or anything that needs a grid.
//! A union of spheres that reads as one skin is marching cubes, which is mesh generation of the
//! kind [the library deliberately does not
//! have](https://github.com/ypark-dev/pantometry/blob/main/CLAUDE.md) — and an example that grew
//! one would be the wrong place for it.

#![allow(dead_code)]

use std::collections::BTreeMap;

/// Triangles, the points they name, and a value per point to colour it by.
///
/// The three arrays are what `PanelData::surface` takes, in that order, so a mesh goes into a
/// panel without rearranging anything.
#[derive(Clone, Debug, Default)]
pub struct Mesh {
    /// Every vertex, in world coordinates.
    pub points: Vec<[f64; 3]>,
    /// Three indices into `points` per triangle, wound counter-clockwise seen from outside.
    pub faces: Vec<[u32; 3]>,
    /// One value per vertex.
    pub values: Vec<f64>,
}

impl Mesh {
    /// Put `other` alongside, renumbering its triangles.
    ///
    /// The result is one mesh holding two solids, which is what a protein and its ligand are.
    /// They stay separate *bodies* — nothing is merged, welded or intersected — so the
    /// closed-solid property holds per element and [`Mesh::unshared_edges`] is asked of each
    /// piece before they are joined, not after.
    pub fn append(&mut self, other: &Mesh) {
        let base = self.points.len() as u32;
        self.points.extend_from_slice(&other.points);
        self.values.extend_from_slice(&other.values);
        self.faces
            .extend(other.faces.iter().map(|t| t.map(|i| i + base)));
    }

    /// The total area of every triangle.
    pub fn area(&self) -> f64 {
        self.faces
            .iter()
            .map(|t| {
                let (a, b, c) = (
                    self.points[t[0] as usize],
                    self.points[t[1] as usize],
                    self.points[t[2] as usize],
                );
                0.5 * cross(sub(b, a), sub(c, a))
                    .iter()
                    .map(|x| x * x)
                    .sum::<f64>()
                    .sqrt()
            })
            .sum()
    }

    /// The volume the triangles enclose, by the divergence theorem.
    ///
    /// **Signed, on purpose.** `(1/6) * sum v0 . (v1 x v2)` is positive when the winding faces
    /// outward and negative when it does not, so this is the one number that can tell an
    /// inside-out solid from a right-side-out one. A renderer cannot: a reversed mesh lights as a
    /// surface either way.
    ///
    /// Taken about the mesh's own mean point rather than about the origin. For a closed surface
    /// the two are the same number in exact arithmetic — the volume of a closed thing does not
    /// know where zero is — and better conditioned about the mean in floating point, because the
    /// terms are cubes of the distance from whatever point they are taken about.
    ///
    /// **Measured, rather than asserted, and it buys less than the first version of this comment
    /// claimed.** Adenylate kinase sits 5.74e-9 m from the origin and is 4.04e-9 m across, not
    /// the 3e-9 and 1e-10 that stood here — 1e-10 m is one atom. The backbone tube's
    /// `sum|term| / |sum term|` is 37.2 about the origin and **10.4** about the mean: a factor of
    /// 3.6, not a rescue. Ten times cancellation is still there and no tolerance here carries it,
    /// which is safe only because the one volume asked of a curved tube is asked for its sign and
    /// a ten-per-cent band. For the solids checked against a closed form the ratio is 1.000
    /// either way, because they are built at the origin and the centring buys nothing at all.
    pub fn volume(&self) -> f64 {
        if self.points.is_empty() {
            return 0.0;
        }
        let mut mean = [0.0; 3];
        for p in &self.points {
            for k in 0..3 {
                mean[k] += p[k] / self.points.len() as f64;
            }
        }
        self.faces
            .iter()
            .map(|t| {
                let (a, b, c) = (
                    sub(self.points[t[0] as usize], mean),
                    sub(self.points[t[1] as usize], mean),
                    sub(self.points[t[2] as usize], mean),
                );
                dot(a, cross(b, c)) / 6.0
            })
            .sum()
    }

    /// How many edges are **not** shared by exactly two triangles. Zero means closed.
    ///
    /// Counted by *position* rather than by index, as `optical_bench` counts it: a crease is two
    /// faces that share a place and not an index, so counting indices reports every deliberate
    /// duplicate as a hole. Positions compare by bits, because the two copies are the same `f64`
    /// written twice and a tolerance here would be a choice about how close is the same.
    pub fn unshared_edges(&self) -> usize {
        let key = |p: [f64; 3]| p.map(f64::to_bits);
        let mut edges: BTreeMap<([u64; 3], [u64; 3]), usize> = BTreeMap::new();
        for t in &self.faces {
            for k in 0..3 {
                let (a, b) = (
                    key(self.points[t[k] as usize]),
                    key(self.points[t[(k + 1) % 3] as usize]),
                );
                *edges
                    .entry(if a <= b { (a, b) } else { (b, a) })
                    .or_default() += 1;
            }
        }
        edges.values().filter(|&&n| n != 2).count()
    }
}

/// A circular cross-section swept along `path`, capped at both ends.
///
/// `values` is one per path point and is carried to every vertex of that point's ring, so a tube
/// drawn along a backbone is coloured per residue with no interpolation invented here — the
/// renderer interpolates between rings, which is what makes a smooth gradient out of 214 numbers.
///
/// **The frame is parallel-transported, not rebuilt.** Choosing a normal afresh at each point by
/// crossing the tangent with a fixed axis makes the tube spin wherever the path turns through
/// that axis, and on a protein backbone it turns through every axis. Each ring's normal is the
/// previous one with its tangential part removed, which is the discrete parallel transport: the
/// minimum rotation that keeps it perpendicular. The seam therefore runs straight along the tube
/// and the triangles do not shear.
///
/// # Panics
///
/// If `path` is shorter than two points, `sides` is less than three, `values` is not one per
/// point, or two consecutive points are in the same place — a path with no direction at a point
/// has no ring there, and a silent `NaN` would travel all the way to a picture that looks like a
/// gap.
pub fn tube(path: &[[f64; 3]], values: &[f64], radius: f64, sides: usize) -> Mesh {
    assert!(
        path.len() >= 2,
        "a tube needs a path, and this has {} points",
        path.len()
    );
    assert!(sides >= 3, "a cross-section needs three sides, not {sides}");
    assert_eq!(
        values.len(),
        path.len(),
        "a tube is coloured per path point, and this has {} values for {} points",
        values.len(),
        path.len()
    );

    // Each segment's own direction, as a unit vector.
    let n = path.len();
    let along: Vec<[f64; 3]> = (0..n - 1)
        .map(|i| {
            let d = sub(path[i + 1], path[i]);
            let len = length(d);
            assert!(
                len > 0.0 && len.is_finite(),
                "points {i} and {} of this path are the same place, so the segment between them \
                 has no direction and no ring",
                i + 1
            );
            scale(d, 1.0 / len)
        })
        .collect();

    // The cut plane at an interior point is the one that **bisects** the two segment directions,
    // which is `unit(a) + unit(b)` and not the chord `p[i+1] - p[i-1]`. The chord is the bisector
    // only where the two steps are the same length: an alpha-carbon trace steps 3.8 A each time
    // but not exactly, and off the bisector the one ellipse below cannot be the right cut for
    // both segments at once. The ends have one segment and take it.
    let tangents: Vec<[f64; 3]> = (0..n)
        .map(|i| match i {
            0 => along[0],
            i if i == n - 1 => along[n - 2],
            i => {
                let s = add(along[i - 1], along[i]);
                assert!(
                    length(s) > 1e-9,
                    "this path doubles back exactly at point {i}, where the cut plane that \
                     bisects the two segments contains them both"
                );
                unit(s)
            }
        })
        .collect();

    // Seed the frame off whichever axis the first tangent is least aligned with, so the cross
    // product is never the one that is about to vanish.
    let t0 = tangents[0];
    let axis = {
        let m = t0.map(f64::abs);
        if m[0] <= m[1] && m[0] <= m[2] {
            [1.0, 0.0, 0.0]
        } else if m[1] <= m[2] {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        }
    };
    let mut normal = unit(cross(t0, axis));

    // The direction of the segment arriving at each point, which is the axis the tube has to keep
    // its radius about. The ends have only one segment, so theirs is the ring's own normal and
    // the miter below does nothing there.
    let arriving: Vec<[f64; 3]> = (0..n)
        .map(|i| if i == 0 { along[0] } else { along[i - 1] })
        .collect();

    let mut mesh = Mesh::default();
    let mut tightest = 1.0f64;
    for (i, &t) in tangents.iter().enumerate() {
        // Parallel transport: drop the part of the previous normal that has become tangential.
        // A turn of nearly 180 degrees would leave nothing to renormalise, which a backbone does
        // not do -- but a path that doubled back exactly would, so it is reseeded rather than
        // divided by zero.
        let projected = sub(normal, scale(t, dot(normal, t)));
        normal = if length(projected) > 1e-6 {
            unit(projected)
        } else {
            unit(cross(t, axis))
        };
        let binormal = cross(t, normal);

        // **The miter.** A circle laid flat on the cut plane is not the cut: a cylinder of radius
        // `radius` meets a plane tilted by `phi` in an *ellipse*, long by `1 / cos phi` along the
        // steepest direction of the tilt. Laying a circle there instead pinches the tube to
        // `radius * cos phi` at every turn, which is a solid 24% short on a protein backbone and
        // reads in a picture as a chain of beads.
        let cos_phi = dot(arriving[i], t);
        tightest = tightest.min(cos_phi);
        let steepest = sub(arriving[i], scale(t, cos_phi));
        let stretch = if length(steepest) > 1e-9 && cos_phi > 0.0 {
            Some((unit(steepest), 1.0 / cos_phi - 1.0))
        } else {
            None
        };

        for j in 0..sides {
            let theta = std::f64::consts::TAU * j as f64 / sides as f64;
            let mut r = scale(
                add(scale(normal, theta.cos()), scale(binormal, theta.sin())),
                radius,
            );
            if let Some((axis, extra)) = stretch {
                r = add(r, scale(axis, dot(r, axis) * extra));
            }
            mesh.points.push(add(path[i], r));
            mesh.values.push(values[i]);
        }
    }
    assert!(
        tightest > 0.2,
        "this path turns through {:.0} degrees at one point, and a miter there would stretch the \
         ring {:.1}x and fold the tube through itself — subdivide the path first",
        2.0 * tightest.acos().to_degrees(),
        1.0 / tightest
    );

    // The skin. With (normal, binormal, tangent) right-handed, going round the ring from normal
    // towards binormal is counter-clockwise seen from ahead, and this pair of triangles is the
    // one whose cross product comes out along the ring's own radius -- outward.
    let ring = sides as u32;
    for i in 0..(n as u32 - 1) {
        for j in 0..ring {
            let k = (j + 1) % ring;
            let (a, b) = (i * ring + j, i * ring + k);
            let (c, d) = ((i + 1) * ring + k, (i + 1) * ring + j);
            mesh.faces.push([a, b, c]);
            mesh.faces.push([a, c, d]);
        }
    }

    // Caps, each a fan from the path's own end point -- which is the polygon's centre, so the cap
    // is exactly the regular n-gon whose area the check below knows in closed form.
    let first = mesh.points.len() as u32;
    mesh.points.push(path[0]);
    mesh.values.push(values[0]);
    let last = mesh.points.len() as u32;
    mesh.points.push(path[n - 1]);
    mesh.values.push(values[n - 1]);
    let base = (n as u32 - 1) * ring;
    for j in 0..ring {
        let k = (j + 1) % ring;
        // Reversed at the start, because the outward direction there is against the path.
        mesh.faces.push([first, k, j]);
        mesh.faces.push([last, base + j, base + k]);
    }
    mesh
}

/// A sphere, as an icosahedron subdivided `level` times and pushed out onto its circumsphere.
///
/// Level 0 is the icosahedron itself. Each level quarters every triangle, so the count is
/// `20 * 4^level` and the area approaches `4 pi r^2` from below at the same rate.
///
/// Vertices **are** shared between the faces that meet at them, through a cache on the edge that
/// was split. Costs nothing — [`Mesh::unshared_edges`] counts by position and not by index, so
/// sharing is invisible to it — and saves a factor of five in points, which is the factor the
/// animation's file size is decided by: 57 atoms on 48 frames is the same sphere written 2736
/// times.
pub fn sphere(centre: [f64; 3], radius: f64, level: u32, value: f64) -> Mesh {
    // The icosahedron's twelve vertices are the corners of three golden rectangles.
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let v: Vec<[f64; 3]> = [
        [-1.0, phi, 0.0],
        [1.0, phi, 0.0],
        [-1.0, -phi, 0.0],
        [1.0, -phi, 0.0],
        [0.0, -1.0, phi],
        [0.0, 1.0, phi],
        [0.0, -1.0, -phi],
        [0.0, 1.0, -phi],
        [phi, 0.0, -1.0],
        [phi, 0.0, 1.0],
        [-phi, 0.0, -1.0],
        [-phi, 0.0, 1.0],
    ]
    .iter()
    .map(|&p| unit(p))
    .collect();
    #[rustfmt::skip]
    let faces: [[usize; 3]; 20] = [
        [0, 11, 5], [0, 5, 1], [0, 1, 7], [0, 7, 10], [0, 10, 11],
        [1, 5, 9], [5, 11, 4], [11, 10, 2], [10, 7, 6], [7, 1, 8],
        [3, 9, 4], [3, 4, 2], [3, 2, 6], [3, 6, 8], [3, 8, 9],
        [4, 9, 5], [2, 4, 11], [6, 2, 10], [8, 6, 7], [9, 8, 1],
    ];

    let mut points = v;
    let mut tris: Vec<[u32; 3]> = faces.iter().map(|&f| f.map(|i| i as u32)).collect();
    for _ in 0..level {
        // One new vertex per edge, found again by the pair of ends that made it -- so the two
        // triangles sharing that edge get the same index and the sphere stays one skin.
        let mut split: BTreeMap<(u32, u32), u32> = BTreeMap::new();
        let mut next = Vec::with_capacity(tris.len() * 4);
        for t in &tris {
            let mut middle = [0u32; 3];
            for k in 0..3 {
                let (a, b) = (t[k], t[(k + 1) % 3]);
                let edge = if a < b { (a, b) } else { (b, a) };
                middle[k] = *split.entry(edge).or_insert_with(|| {
                    // Back onto the unit sphere, which is what makes this converge to a sphere
                    // rather than to a finer icosahedron.
                    points.push(unit(add(points[a as usize], points[b as usize])));
                    points.len() as u32 - 1
                });
            }
            next.push([t[0], middle[0], middle[2]]);
            next.push([middle[0], t[1], middle[1]]);
            next.push([middle[2], middle[1], t[2]]);
            next.push(middle);
        }
        tris = next;
    }

    Mesh {
        values: vec![value; points.len()],
        points: points
            .iter()
            .map(|&p| add(centre, scale(p, radius)))
            .collect(),
        faces: tris,
    }
}

/// A path through the same points, with `per_segment` samples of a Catmull-Rom spline between
/// each consecutive pair, and the values carried along it.
///
/// **This invents no measurement.** The spline passes exactly through every control point — that
/// is the defining property of Catmull-Rom and the thing checked against it — so every alpha
/// carbon is still where the crystallographer put it, and what is interpolated is only the curve
/// between two atoms, which no atom occupies. A straight polyline through them invents a path
/// there too; this one just turns less sharply, which is what the miter needs: the ring at a
/// bend widens by `1 / cos phi`, so a trace that changes direction by 60 degrees in one 3.8 A
/// step makes a wide flat lozenge where a chain should be.
///
/// **Centripetal** — knots spaced by the square root of the step rather than evenly — which is
/// the parameterisation with no cusp and no self-intersection inside a segment, and costs three
/// lengths per segment to have.
///
/// It is *not* what made this backbone drawable, and the first version of this comment said it
/// was. The tightest curve on the smoothed trace measured 1.79 tube radii with uniform knots and
/// 1.79 with centripetal ones: the tight turns are in the data, and no parameterisation was going
/// to move them. What the raw polyline reported — 0.65, comfortably safe — was the same estimator
/// on a 3.8 Å sampling, and three points that far apart cannot resolve a radius of 1 Å. A coarser
/// grid reporting a gentler curve is not a gentler curve.
///
/// The ends are handled by reflecting the second point through the first, which is the standard
/// phantom-point rule and keeps the curve's end tangent along its own first segment.
pub fn smoothed(
    path: &[[f64; 3]],
    values: &[f64],
    per_segment: usize,
) -> (Vec<[f64; 3]>, Vec<f64>) {
    assert!(
        per_segment >= 1,
        "a segment cannot be cut into {per_segment} pieces"
    );
    assert_eq!(values.len(), path.len(), "one value per point");
    let n = path.len();
    let at = |i: isize| -> [f64; 3] {
        if i < 0 {
            sub(scale(path[0], 2.0), path[1])
        } else if i as usize >= n {
            sub(scale(path[n - 1], 2.0), path[n - 2])
        } else {
            path[i as usize]
        }
    };
    let mut out = Vec::with_capacity((n - 1) * per_segment + 1);
    let mut carried = Vec::with_capacity((n - 1) * per_segment + 1);
    for i in 0..n - 1 {
        let p = [at(i as isize - 1), path[i], path[i + 1], at(i as isize + 2)];
        // Knots spaced by the square root of the step: alpha = 1/2, the centripetal choice.
        let mut knot = [0.0; 4];
        for k in 1..4 {
            let step = length(sub(p[k], p[k - 1])).sqrt();
            // Not clamped to MIN_POSITIVE, which is a guard that produces NaN: a knot interval of
            // 2.2e-308 against a parameter of order 1e-5 divides to infinity, and the coordinates
            // come out NaN to be caught further downstream with a message about something else.
            assert!(
                step > 0.0 && step.is_finite(),
                "points {} and {} of this path are the same place, so the spline through them \
                 has no parameter",
                i + k - 1,
                i + k
            );
            knot[k] = knot[k - 1] + step;
        }
        for s in 0..per_segment {
            let u = s as f64 / per_segment as f64;
            let t = knot[1] + u * (knot[2] - knot[1]);
            // Barry and Goldman's pyramid: three lerps, then two, then one. At t = knot[1] the
            // last lerp returns p[1] untouched, which is the interpolating property this file
            // checks the spline against.
            let lerp = |a: [f64; 3], b: [f64; 3], t0: f64, t1: f64| {
                let w = (t - t0) / (t1 - t0);
                add(scale(a, 1.0 - w), scale(b, w))
            };
            let a1 = lerp(p[0], p[1], knot[0], knot[1]);
            let a2 = lerp(p[1], p[2], knot[1], knot[2]);
            let a3 = lerp(p[2], p[3], knot[2], knot[3]);
            let b1 = lerp(a1, a2, knot[0], knot[2]);
            let b2 = lerp(a2, a3, knot[1], knot[3]);
            out.push(lerp(b1, b2, knot[1], knot[2]));
            carried.push(values[i] * (1.0 - u) + values[i + 1] * u);
        }
    }
    out.push(path[n - 1]);
    carried.push(values[n - 1]);
    (out, carried)
}

/// The closest two samples of a path come, counting only pairs at least `apart` indices apart.
///
/// **[`crowding`] is local and this is not.** A tube folds at a bend, which is what crowding
/// measures; two *straight* stretches of the same path, far apart along it, can pass within a
/// tube diameter of each other while nothing bends at all. On a protein that is a helix packing
/// against a sheet, and it is the ordinary case rather than a pathology. A tube fat enough to
/// span the gap fuses the two into one body, which is closed, right-side-out, keeps its radius
/// everywhere and does not fold — every other check here passes on it.
///
/// `apart` exists because neighbouring samples are *supposed* to be within a diameter: that is
/// what makes the tube continuous. It has to exclude at least the stretch of path a tube's width
/// covers, and beyond that the exclusion is crowding's territory.
pub fn nearest_approach(path: &[[f64; 3]], apart: usize) -> f64 {
    let mut closest = f64::MAX;
    for i in 0..path.len() {
        for j in (i + apart)..path.len() {
            closest = closest.min(length(sub(path[j], path[i])));
        }
    }
    closest
}

/// The tube's radius against the tightest curve the path takes: below one it has not folded
/// through itself.
///
/// **The condition is curvature, not step length.** A swept tube of radius `r` following a curve
/// whose radius of curvature is `R` closes up on the inside of the bend and passes through its
/// own axis when `r` reaches `R`. Step length has nothing to do with it: consecutive rings on a
/// straight path are parallel discs and never meet however finely it is cut — a first version of
/// this compared the ring against the step, which is a criterion that *tightens* when a path is
/// subdivided and reported a smoothed trace as three times worse than the polyline it smoothed.
///
/// The radius of curvature at an interior point is the circumradius of the triangle its two
/// neighbours make with it, `abc / 4T`, which is exact for three points on a circle. Collinear
/// points give zero area and an infinite radius, which is the right answer there.
///
/// Every other check passes on a folded tube — it is closed, its rings keep their radius, and it
/// renders as a bead rather than a chain — so this is the only number that sees it, and it has to
/// be asked of **every conformation drawn** and not only the one at rest.
pub fn crowding(path: &[[f64; 3]], radius: f64) -> f64 {
    let tightest = (1..path.len().saturating_sub(1))
        .map(|i| {
            let (u, v) = (sub(path[i], path[i - 1]), sub(path[i + 1], path[i]));
            let twice_area = length(cross(u, v));
            if twice_area == 0.0 {
                f64::MAX
            } else {
                length(u) * length(v) * length(add(u, v)) / (2.0 * twice_area)
            }
        })
        .fold(f64::MAX, f64::min);
    radius / tightest
}

/// The largest distance by which any ring vertex misses the tube's radius about its own segment.
///
/// **This is what makes a tube a tube**, and the one property of a mitered sweep that a picture
/// cannot show: a tube that pinches at every turn still renders as a smooth solid, and the
/// closed-solid count still reads zero. For every vertex of the ring at path point `i`, the
/// distance to the *line* of the segment arriving there — and to the line of the segment leaving
/// — has to be exactly `radius`. Both, because the cut plane bisects, so one ellipse has to serve
/// both segments; against the chord instead of the bisector it serves neither, and against a
/// circle instead of an ellipse it serves neither.
///
/// Exact, not a limit: the only error in it is floating point.
pub fn worst_radius_error(path: &[[f64; 3]], mesh: &Mesh, radius: f64, sides: usize) -> f64 {
    let mut worst: f64 = 0.0;
    for i in 0..path.len() {
        for axis in [
            unit(sub(path[i.max(1)], path[i.max(1) - 1])),
            unit(sub(
                path[(i + 1).min(path.len() - 1)],
                path[i.min(path.len() - 2)],
            )),
        ] {
            for j in 0..sides {
                let v = sub(mesh.points[i * sides + j], path[i]);
                let perpendicular = sub(v, scale(axis, dot(v, axis)));
                worst = worst.max((length(perpendicular) - radius).abs());
            }
        }
    }
    worst
}

/// The exact area and volume of the prism [`tube`] draws along a straight path.
///
/// A regular `sides`-gon of circumradius `radius` has perimeter `2 n r sin(pi/n)` and area
/// `(n/2) r^2 sin(2 pi/n)`; swept `length` along its own normal, with a flat cap at each end,
/// that is the whole solid. Exact — not a limit — which is what makes it a check rather than a
/// second opinion, and what lets the same numbers be read again as `n` grows to watch them
/// approach `2 pi r L + 2 pi r^2` and `pi r^2 L`.
pub fn straight_tube_exactly(radius: f64, length: f64, sides: usize) -> (f64, f64) {
    let n = sides as f64;
    let perimeter = 2.0 * n * radius * (std::f64::consts::PI / n).sin();
    let cross_section = 0.5 * n * radius * radius * (std::f64::consts::TAU / n).sin();
    (
        perimeter * length + 2.0 * cross_section,
        cross_section * length,
    )
}

/// The exact area and volume of the icosahedron [`sphere`] starts from, at level 0.
///
/// With circumradius `r` the edge is `4r / sqrt(10 + 2 sqrt 5)`; the area is `5 sqrt(3) a^2` and
/// the volume `(5/12)(3 + sqrt 5) a^3`.
pub fn icosahedron_exactly(radius: f64) -> (f64, f64) {
    let a = 4.0 * radius / (10.0 + 2.0 * 5.0_f64.sqrt()).sqrt();
    (
        5.0 * 3.0_f64.sqrt() * a * a,
        (5.0 / 12.0) * (3.0 + 5.0_f64.sqrt()) * a * a * a,
    )
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

fn unit(a: [f64; 3]) -> [f64; 3] {
    scale(a, 1.0 / length(a))
}
