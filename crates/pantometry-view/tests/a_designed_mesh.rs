//! Triangles from a file, against a shape whose area and volume are known exactly.
//!
//! `mesh_surface` is the one producer in `mesh` that does not *derive* a surface — it carries
//! triangles across. That makes it the easiest to get subtly wrong and the easiest to check
//! honestly: a cube of side `s` has area `6s²` and encloses `s³`, both exact, and both computed
//! here from the emitted [`Surface`] rather than from the triangles that went in. Reading the
//! answer off the input would test nothing.
//!
//! The side is **40 mm** on purpose. A 40 mm cube came out 80 mm across once, when the cell
//! mapping in this module gave every sample a full cell either side, and that is the failure this
//! family of tests exists to catch.

use pantometry_view::mesh::{mesh_surface, Surface};

/// The cube's side, in metres. It spans `0..S` on every axis, with a corner **at** the origin —
/// which two of the checks below have to reason about, so it is stated here rather than inferred.
const S: f64 = 0.040;

/// The twelve triangles of a cube spanning `0..S` on every axis, wound counter-clockwise seen
/// from outside — the convention every STL writer means.
fn cube() -> Vec<[[f64; 3]; 3]> {
    let v = [
        [0.0, 0.0, 0.0],
        [S, 0.0, 0.0],
        [S, S, 0.0],
        [0.0, S, 0.0],
        [0.0, 0.0, S],
        [S, 0.0, S],
        [S, S, S],
        [0.0, S, S],
    ];
    // Each row is one face as a quad, already outward: the two triangles are (0,1,2) and (0,2,3)
    // of the row. Written this way so a reader can check six faces rather than twelve triangles.
    let quads = [
        [0, 3, 2, 1], // z = 0, outward -z
        [4, 5, 6, 7], // z = S, outward +z
        [0, 1, 5, 4], // y = 0, outward -y
        [3, 7, 6, 2], // y = S, outward +y
        [0, 4, 7, 3], // x = 0, outward -x
        [1, 2, 6, 5], // x = S, outward +x
    ];
    let mut out = Vec::new();
    for q in quads {
        out.push([v[q[0]], v[q[1]], v[q[2]]]);
        out.push([v[q[0]], v[q[2]], v[q[3]]]);
    }
    out
}

/// The faces of a surface, as triples of `f64` positions.
fn faces(s: &Surface) -> Vec<[[f64; 3]; 3]> {
    s.indices
        .chunks_exact(3)
        .map(|f| {
            let p = |i: u32| {
                let q = s.positions[i as usize];
                [q[0] as f64, q[1] as f64, q[2] as f64]
            };
            [p(f[0]), p(f[1]), p(f[2])]
        })
        .collect()
}

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Total area, from the emitted triangles.
fn area(s: &Surface) -> f64 {
    faces(s)
        .iter()
        .map(|t| {
            let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
            dot(n, n).sqrt() / 2.0
        })
        .sum()
}

/// Enclosed volume by the divergence theorem: a sixth of the sum of `a · (b × c)`.
///
/// **Signed**, which is what makes it a check of the winding as well as of the size. A face
/// pointing the wrong way subtracts its tetrahedron instead of adding it, so a mesh with mixed
/// winding collapses toward zero and one that is wholly inside-out comes out negative.
fn volume(s: &Surface) -> f64 {
    faces(s)
        .iter()
        .map(|t| dot(t[0], cross(t[1], t[2])) / 6.0)
        .sum::<f64>()
}

/// Positions cross to `f32`, whose mantissa is 24 bits — a relative error of `2⁻²⁴ ≈ 6e-8` per
/// coordinate. Area is quadratic in them and volume cubic, so a few multiples of that is the
/// floor; `1e-6` is comfortably above it and far below any error of construction.
const F32: f64 = 1e-6;

#[test]
fn a_cube_has_the_area_and_volume_of_a_cube() {
    let s = mesh_surface(&cube(), 0);
    assert_eq!(s.faces(), 12, "twelve triangles in, twelve out");

    let (want_area, want_volume) = (6.0 * S * S, S * S * S);
    let (got_area, got_volume) = (area(&s), volume(&s));

    assert!(
        (got_area - want_area).abs() / want_area < F32,
        "area {got_area} against {want_area}"
    );
    assert!(
        (got_volume - want_volume).abs() / want_volume < F32,
        "volume {got_volume} against {want_volume}"
    );
}

#[test]
fn the_winding_survives_the_crossing() {
    // The volume test above already fails on a flipped face, but it fails by the wrong *size*
    // and a reader could mistake that for a scaling bug. This says what the failure means: a
    // negative volume is an inside-out mesh, and one near zero is a mesh at war with itself.
    let good = volume(&mesh_surface(&cube(), 0));
    assert!(good > 0.0, "the mesh is inside-out");

    // Wholly inside-out: every term negates, so the answer is exactly the volume with a minus
    // sign. This is the unambiguous case and it does not depend on where the origin is.
    let inverted: Vec<_> = cube()
        .into_iter()
        .map(|mut t| {
            t.swap(1, 2);
            t
        })
        .collect();
    let back_to_front = volume(&mesh_surface(&inverted, 0));
    assert!(
        (back_to_front + good).abs() / good < F32,
        "inverted came to {back_to_front}, not {}",
        -good
    );

    // One face of twelve. **Triangle 2 and not triangle 0**, and the difference is a blind spot
    // worth naming: this cube has a corner *at* the origin, and `a · (b × c)` is zero whenever
    // `a` is the origin, so the two triangles touching it contribute nothing and flipping either
    // changes the sum by nothing. The check is exact for a closed mesh at any origin; a mesh
    // broken open by a single flip is origin-dependent, and this origin cannot see two of the
    // twelve. Triangle 2 is on the far face and carries `S³/6`, so flipping it costs a third.
    let mut one = cube();
    one[2].swap(1, 2);
    let partial = volume(&mesh_surface(&one, 0));
    let want = good - 2.0 * (S * S * S) / 6.0;
    assert!(
        (partial - want).abs() / good < F32,
        "flipping one face gave {partial}, not {want}"
    );
}

#[test]
fn every_face_points_away_from_the_inside() {
    // For a convex body about its own centre, outward is the direction the centroid lies in.
    // This is what back-face culling and every lighting model read, and it is not implied by the
    // volume being right -- a mesh can enclose the correct volume with its normals detached from
    // its winding, if the normals were computed some other way.
    let s = mesh_surface(&cube(), 0);
    let centre = [S / 2.0, S / 2.0, S / 2.0];
    for (f, t) in faces(&s).iter().enumerate() {
        let centroid = [
            (t[0][0] + t[1][0] + t[2][0]) / 3.0,
            (t[0][1] + t[1][1] + t[2][1]) / 3.0,
            (t[0][2] + t[1][2] + t[2][2]) / 3.0,
        ];
        let n = s.normals[s.indices[f * 3] as usize];
        let n = [n[0] as f64, n[1] as f64, n[2] as f64];
        let outward = dot(n, sub(centroid, centre));
        assert!(outward > 0.0, "face {f} points inward: {outward}");
        assert!(
            (dot(n, n) - 1.0).abs() < 1e-6,
            "face {f}'s normal is not a unit vector"
        );
    }
}

#[test]
fn faces_share_nothing_and_a_cube_has_six_directions() {
    let s = mesh_surface(&cube(), 0);
    assert_eq!(
        s.positions.len(),
        36,
        "three vertices to a triangle, unshared"
    );
    assert_eq!(s.normals.len(), 36);
    assert_eq!(s.source.len(), 36);

    // The three vertices of a face agree, and the six faces do not: a cube shaded with shared
    // corner normals has one averaged direction per corner and reads as a sphere.
    for t in 0..12 {
        let n = s.normals[t * 3];
        assert_eq!(s.normals[t * 3 + 1], n, "triangle {t} is not flat");
        assert_eq!(s.normals[t * 3 + 2], n);
    }
    let mut directions: Vec<[u32; 3]> = s
        .normals
        .iter()
        .map(|n| [n[0].to_bits(), n[1].to_bits(), n[2].to_bits()])
        .collect();
    directions.sort_unstable();
    directions.dedup();
    assert_eq!(directions.len(), 6, "a cube points six ways");
}

#[test]
fn a_triangle_with_no_normal_is_dropped() {
    // Three ways a face can have no direction. An exporter emits the first two routinely, and
    // any of them shades a face black if it reaches a renderer.
    let degenerate = [
        // Two vertices in the same place.
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        // Three points on a line.
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        // A vertex that is not a number. It poisons the cross product to NaN, which is why the
        // one guard on the normal catches this as well as the other two.
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, f64::NAN, 0.0]],
    ];
    for (i, t) in degenerate.iter().enumerate() {
        let s = mesh_surface(std::slice::from_ref(t), 0);
        assert_eq!(s.faces(), 0, "case {i} was kept");
        assert!(s.positions.is_empty(), "case {i} left vertices behind");
    }

    // A good triangle among bad ones survives, and the indices still address the right vertex --
    // dropping a face after allocating its base index is how a skip becomes a corruption.
    let mixed = [degenerate[0], cube()[0], degenerate[1]];
    let s = mesh_surface(&mixed, 0);
    assert_eq!(s.faces(), 1);
    assert_eq!(s.indices, vec![0, 1, 2]);
    assert_eq!(s.positions.len(), 3);
}

#[test]
fn a_mesh_of_nothing_comes_out_empty_rather_than_broken() {
    // The silent-failure shape: an empty surface is a legitimate answer here and a renderer will
    // draw nothing, so the contract is that it is *well formed* when empty. A `Surface` with
    // indices and no positions would panic in the exporter instead.
    let s = mesh_surface(&[], 0);
    assert_eq!(s.faces(), 0);
    assert!(s.positions.is_empty() && s.normals.is_empty() && s.indices.is_empty());
    assert_eq!(s.stride, 1, "stride one: nothing here is ever subsampled");
}

#[test]
fn every_vertex_carries_the_callers_source() {
    // One index for the whole mesh, not one per vertex. The exporters look each vertex's colour
    // up by this, and a part is one object with one material.
    let s = mesh_surface(&cube(), 7);
    assert!(s.source.iter().all(|&i| i == 7));
    assert_eq!(s.source.len(), s.positions.len(), "one per vertex");
}
