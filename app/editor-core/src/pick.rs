//! **Pointing at a model: which face of a domain a pixel is over, which cell of it, and the edit
//! that writes what was pointed at.**
//!
//! A cooled face, a held one, a patch where a part is bolted and a box of cells that dissipates
//! were all chosen by typing a word or six indices into the settings form, and reading them back
//! meant counting cells in one's head. A simulation tool is set up by pointing at the model, so
//! that is what this is: a ray from [`viewer_core::Camera::ray`], the domain's box in its own frame,
//! and the pose that places it. No event loop, so the browser editor can make the same calls.
//!
//! # In the domain's own frame
//!
//! A turned domain's faces are not the world's: the face a scene calls `x-max` is wherever the
//! pose has put it. So the ray is taken into the domain's frame — the pose undone, which is
//! [`pantometry::scene::Placed::apply`] read backwards — and the box is the one the grid is stated
//! in. A face found there is the face the file means, however the domain is turned.
//!
//! # What an area is
//!
//! A cooled face states its `area_cm2`, and the block divides it among the **solid** cells on
//! the face, or on the patch when there is one — a void cell takes no share. So the area a pick
//! writes is the solid cells it covers times a cell's face: the whole face of a plain block,
//! which is what [`pantometry_world::templates::entry_for`] writes for a new entry, to the bit; and
//! for a designed part, the cells the part actually puts there. Scene 29's bracket states
//! 16.5 cm² for a top whose grid face is 27.04, which is this rule, written down by hand.

use crate::edit;
use pantometry::scene::Placed;
use pantometry_world::templates::{self, Placing};
use pantometry_world::FaceSpec;

/// A domain a click can land on: a grid of cells, its box in its own frame, and where it is.
#[derive(Clone, Debug, PartialEq)]
pub struct Target {
    /// The domain's name.
    pub name: String,
    /// Its `kind`, as the file spells it: which lists can be pointed at is the format's answer for
    /// the kind, in [`templates::SITES`].
    pub kind: String,
    /// The box, `[x0, y0, z0, x1, y1, z1]`, in metres, in the domain's own frame.
    pub extent: [f64; 6],
    /// Cells along each axis.
    pub cells: [usize; 3],
    /// A cell's side, as the file states it, in millimetres.
    pub cell_mm: f64,
    /// Where the domain's frame sits in the world.
    pub place: Placed,
    /// Which cells hold nothing, `i + nx (j + ny k)`: from the build, or all solid for a domain
    /// nothing in the file can make void. `None` for a block with parts or a `void` region in a
    /// scene that does not build — its solid cells are not known, and [`set_face`] refuses to
    /// write an area for it rather than count every cell.
    pub void: Option<Vec<bool>>,
}

/// Where a ray met a box.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Hit {
    /// The face it came in through.
    pub face: FaceSpec,
    /// Where, in the domain's frame, on the face exactly.
    pub local: [f64; 3],
    /// How far along the ray, in metres.
    pub distance: f64,
}

/// A box of cells on a face, `from` inclusive and `to` exclusive, in the face's two axes.
pub type Patch = ([usize; 2], [usize; 2]);

/// A point in the world, in the frame `place` puts somewhere: [`Placed::apply`] undone.
pub fn to_local(place: Placed, p: [f64; 3]) -> [f64; 3] {
    turn_back(
        place,
        [
            p[0] - place.at_m[0],
            p[1] - place.at_m[1],
            p[2] - place.at_m[2],
        ],
    )
}

/// A direction in the world, in the frame `place` turns: the conjugate rotation.
fn turn_back(place: Placed, v: [f64; 3]) -> [f64; 3] {
    let [x, y, z, w] = place.turn;
    Placed {
        at_m: [0.0; 3],
        turn: [-x, -y, -z, w],
    }
    .apply(v)
}

/// Where `ray` first meets the box `extent`, which `place` puts in the world.
///
/// The slab method, in the box's own frame: along each axis the ray is inside between two
/// crossings, and inside the box where all three overlap. The face it enters through is the axis
/// whose crossing is last to begin. An eye inside the box sees the face it leaves through, which
/// is the face in front of it. `None` when the ray misses, when the box is behind the eye, and for a
/// ray that is not a number.
///
/// The point is put **on** the face: its coordinate across the face is the face's own, and the two
/// in it are held to the box, so a hit at an edge is on the face it says and a cell computed from
/// it is a cell of the grid.
pub fn hit(extent: [f64; 6], place: Placed, ray: &viewer_core::Ray) -> Option<Hit> {
    let o = to_local(place, ray.origin);
    let d = turn_back(place, ray.direction);
    let mut enter = (f64::NEG_INFINITY, None);
    let mut leave = (f64::INFINITY, None);
    for a in 0..3 {
        let (lo, hi) = (extent[a], extent[a + 3]);
        if d[a] == 0.0 {
            // Parallel to this pair of faces: inside them throughout, or never.
            if o[a] < lo || o[a] > hi {
                return None;
            }
            continue;
        }
        let (t_lo, t_hi) = ((lo - o[a]) / d[a], (hi - o[a]) / d[a]);
        let (near, far, near_face, far_face) = if d[a] > 0.0 {
            (t_lo, t_hi, (a, false), (a, true))
        } else {
            (t_hi, t_lo, (a, true), (a, false))
        };
        if near > enter.0 {
            enter = (near, Some(near_face));
        }
        if far < leave.0 {
            leave = (far, Some(far_face));
        }
    }
    // A ray that is not a number never sets a face above — every comparison with a NaN is false —
    // so it misses at the `?` below. This is the box being behind the eye or beside the ray.
    if enter.0 > leave.0 || leave.0 < 0.0 {
        return None;
    }
    let (t, (axis, high)) = if enter.0 >= 0.0 {
        (enter.0, enter.1?)
    } else {
        (leave.0, leave.1?)
    };
    let face = FaceSpec::across(axis, high);
    Some(Hit {
        face,
        local: on_face(
            extent,
            face,
            [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]],
        ),
        distance: t,
    })
}

/// Where `ray` meets the plane of `face`, held to the face: the far end of a drag that began on it.
///
/// A drag is a box on the face it started on, wherever the pointer wanders — off the edge, or over
/// the face beside it — so its end is the ray against that face's plane rather than whatever the
/// ray hits, and a pointer past the edge is the edge. `None` for a ray parallel to the plane or
/// meeting it behind the eye.
pub fn on_plane(
    extent: [f64; 6],
    place: Placed,
    face: FaceSpec,
    ray: &viewer_core::Ray,
) -> Option<[f64; 3]> {
    let o = to_local(place, ray.origin);
    let d = turn_back(place, ray.direction);
    let (axis, high) = face.normal();
    let plane = if high { extent[axis + 3] } else { extent[axis] };
    if d[axis] == 0.0 {
        return None;
    }
    let t = (plane - o[axis]) / d[axis];
    if t.is_nan() || t < 0.0 {
        return None;
    }
    Some(on_face(
        extent,
        face,
        [o[0] + t * d[0], o[1] + t * d[1], o[2] + t * d[2]],
    ))
}

/// `p` put on `face` of `extent`: across it, the face's coordinate; in it, held to the box.
fn on_face(extent: [f64; 6], face: FaceSpec, p: [f64; 3]) -> [f64; 3] {
    let (axis, high) = face.normal();
    let mut out = p;
    for a in 0..3 {
        out[a] = out[a].clamp(extent[a], extent[a + 3]);
    }
    out[axis] = if high { extent[axis + 3] } else { extent[axis] };
    out
}

/// The cell holding `x` along an axis from `lo` to `hi` cut into `n`, held to the grid: a point on
/// the far wall is in the last cell, not one past it.
pub fn cell_along(x: f64, lo: f64, hi: f64, n: usize) -> usize {
    if n == 0 || hi.is_nan() || lo.is_nan() || hi <= lo {
        return 0;
    }
    let f = ((x - lo) / (hi - lo) * n as f64).floor();
    if f >= 0.0 {
        (f as usize).min(n - 1)
    } else {
        0
    }
}

/// The patch two cells on a face span, in the face's axes: from the lower of each to one past the
/// higher, so the cells under both ends are in it.
pub fn patch_between(a: [usize; 2], b: [usize; 2]) -> Patch {
    (
        [a[0].min(b[0]), a[1].min(b[1])],
        [a[0].max(b[0]) + 1, a[1].max(b[1]) + 1],
    )
}

impl Target {
    /// Where `ray` first meets this domain's box.
    pub fn hit(&self, ray: &viewer_core::Ray) -> Option<Hit> {
        hit(self.extent, self.place, ray)
    }

    /// The cell holding a point of the domain's frame, held to the grid.
    pub fn cell_at(&self, local: [f64; 3]) -> [usize; 3] {
        [0, 1, 2].map(|a| cell_along(local[a], self.extent[a], self.extent[a + 3], self.cells[a]))
    }

    /// A cell's place on `face`, in the face's two axes — the reading `from` and `to` take.
    pub fn on(&self, face: FaceSpec, cell: [usize; 3]) -> [usize; 2] {
        let [a, b] = face.axes();
        [cell[a], cell[b]]
    }

    /// The whole of `face`, as a patch: every cell on it.
    pub fn whole(&self, face: FaceSpec) -> Patch {
        let [a, b] = face.axes();
        ([0, 0], [self.cells[a], self.cells[b]])
    }

    /// Whether a cell holds something.
    pub fn is_solid(&self, c: [usize; 3]) -> bool {
        let [nx, ny, _] = self.cells;
        self.void
            .as_ref()
            .and_then(|v| v.get(c[0] + nx * (c[1] + ny * c[2])))
            .is_none_or(|void| !*void)
    }

    /// The solid cells on `face` inside `patch`, which is the number a cooled face's area is
    /// divided among.
    pub fn solid_on(&self, face: FaceSpec, patch: Patch) -> usize {
        let (axis, high) = face.normal();
        let [a, b] = face.axes();
        let layer = if high {
            self.cells[axis].saturating_sub(1)
        } else {
            0
        };
        let mut n = 0;
        for u in patch.0[0]..patch.1[0].min(self.cells[a]) {
            for v in patch.0[1]..patch.1[1].min(self.cells[b]) {
                let mut c = [0; 3];
                c[axis] = layer;
                c[a] = u;
                c[b] = v;
                if self.is_solid(c) {
                    n += 1;
                }
            }
        }
        n
    }

    /// The solid cells on each whole face, in [`FaceSpec::ALL`]'s order: what
    /// [`crate::edit::add_entry_with`] picks a new face by.
    pub fn solid_faces(&self) -> [usize; 6] {
        FaceSpec::ALL.map(|f| self.solid_on(f, self.whole(f)))
    }

    /// The area a cooled `face` should state, in cm², and the solid cells it is made of: the
    /// patch's when there is one, the face's when not. See the module's account of what an area
    /// is.
    ///
    /// The same expression [`templates::entry_for`] writes for a new entry, so for a block with no
    /// void the two are one number to the bit.
    pub fn area_cm2(&self, face: FaceSpec, patch: Option<Patch>) -> (f64, usize) {
        let n = self.solid_on(face, patch.unwrap_or_else(|| self.whole(face)));
        (n as f64 * self.cell_mm * self.cell_mm / 100.0, n)
    }

    /// The four corners of a patch on `face`, in world metres, in order round it.
    pub fn patch_corners(&self, face: FaceSpec, patch: Patch) -> [[f64; 3]; 4] {
        let (axis, high) = face.normal();
        let [a, b] = face.axes();
        let at = |k: usize, i: usize| {
            self.extent[k]
                + (self.extent[k + 3] - self.extent[k]) * i as f64 / self.cells[k].max(1) as f64
        };
        let across = if high {
            self.extent[axis + 3]
        } else {
            self.extent[axis]
        };
        let corner = |u: usize, v: usize| {
            let mut p = [0.0; 3];
            p[axis] = across;
            p[a] = at(a, u);
            p[b] = at(b, v);
            self.place.apply(p)
        };
        let ((u0, v0), (u1, v1)) = ((patch.0[0], patch.0[1]), (patch.1[0], patch.1[1]));
        [
            corner(u0, v0),
            corner(u1, v0),
            corner(u1, v1),
            corner(u0, v1),
        ]
    }

    /// The eight corners of a box of cells, in world metres, in [`crate::PlacedBox`]'s order so
    /// [`crate::EDGES`] draws it.
    pub fn cells_corners(&self, from: [usize; 3], to: [usize; 3]) -> [[f64; 3]; 8] {
        let at = |k: usize, i: usize| {
            self.extent[k]
                + (self.extent[k + 3] - self.extent[k]) * i as f64 / self.cells[k].max(1) as f64
        };
        let mut out = [[0.0; 3]; 8];
        for (n, corner) in out.iter_mut().enumerate() {
            let p = [0, 1, 2].map(|k| at(k, if n & (1 << k) == 0 { from[k] } else { to[k] }));
            *corner = self.place.apply(p);
        }
        out
    }
}

/// What a selected entry states about where it is, read from the text.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Stated {
    /// A face, and a patch on it when the entry states one.
    Face {
        /// The face.
        face: FaceSpec,
        /// The patch, or `None` for the whole face.
        patch: Option<Patch>,
    },
    /// A box of cells.
    Cells {
        /// The first cell, inclusive.
        from: [usize; 3],
        /// One past the last.
        to: [usize; 3],
    },
}

/// The domain, list and index a tree path names, when it names an entry: `/domains/<d>/<key>/<j>`.
pub fn entry_of(path: &str) -> Option<(String, String, usize)> {
    let rest = path.strip_prefix("/domains/")?;
    let (head, index) = rest.rsplit_once('/')?;
    let (domain, key) = head.rsplit_once('/')?;
    let index = index.parse().ok()?;
    (!domain.is_empty() && !key.is_empty()).then(|| (domain.to_string(), key.to_string(), index))
}

/// How the entry at a tree path says where it is, if it says: the format's answer for its domain's
/// kind and its list, from [`templates::SITES`].
pub fn placing_of(text: &str, path: &str) -> Option<Placing> {
    let (domain, key, _) = entry_of(path)?;
    let root: serde_json::Value = serde_json::from_str(text).ok()?;
    let (_, d) = edit::domain_at(&root, &domain).ok()?;
    templates::placing(d.get("kind")?.as_str()?, &key)
}

/// A JSON array of `N` whole numbers.
fn indices<const N: usize>(v: Option<&serde_json::Value>) -> Option<[usize; N]> {
    let list: Vec<usize> = v?
        .as_array()?
        .iter()
        .map(|x| x.as_u64().map(|n| n as usize))
        .collect::<Option<_>>()?;
    list.try_into().ok()
}

/// What the entry at a tree path states about where it is. `None` for a path that is not an entry,
/// an entry of a list that does not say where, and one whose words are not the format's.
pub fn stated(text: &str, path: &str) -> Option<Stated> {
    let placing = placing_of(text, path)?;
    let (domain, key, index) = entry_of(path)?;
    let root: serde_json::Value = serde_json::from_str(text).ok()?;
    let (_, d) = edit::domain_at(&root, &domain).ok()?;
    let entry = d.get(&key)?.as_array()?.get(index)?;
    match placing {
        Placing::Face | Placing::FacePatch => {
            let word = entry.get("face")?.as_str()?;
            let face = FaceSpec::ALL.into_iter().find(|f| f.spelled() == word)?;
            let patch = match (
                indices::<2>(entry.get("from")),
                indices::<2>(entry.get("to")),
            ) {
                (Some(a), Some(b)) => Some((a, b)),
                _ => None,
            };
            Some(Stated::Face { face, patch })
        }
        Placing::Cells => Some(Stated::Cells {
            from: indices::<3>(entry.get("from"))?,
            to: indices::<3>(entry.get("to"))?,
        }),
    }
}

/// `member` put in the object at `object`: on its line when the object is written on one, which is
/// how a cooled face is written, and on a line of its own otherwise.
fn insert_member(
    text: &str,
    object: std::ops::Range<usize>,
    member: &str,
) -> Result<String, String> {
    if text[object.clone()].contains('\n') {
        return edit::append_member(text, object, member);
    }
    let inner = &text[object.start + 1..object.end - 1];
    if inner.trim().is_empty() {
        return Ok(edit::splice(
            text,
            object.start + 1..object.end - 1,
            &format!(" {member} "),
        ));
    }
    let last = object.start + 1 + inner.trim_end().len();
    Ok(edit::splice(text, last..last, &format!(", {member}")))
}

/// The object at `object` without its member `key`, and without the comma that held it: the bytes
/// [`insert_member`] and `append_member` write, so adding a member and taking it out gives back
/// the file it started from.
fn cut_member(text: &str, object: &str, key: &str) -> Result<String, String> {
    let pointer = format!("{object}/{}", edit::escape(key));
    let Some(value) = edit::value_span(text, &pointer) else {
        return Ok(text.to_string());
    };
    let b = text.as_bytes();
    // Back from the value: the colon, then the key's closing and opening quotes. The keys cut here
    // are the format's own `from` and `to`, which hold no escaped quote.
    let mut i = value.start;
    while i > 0 && b[i - 1].is_ascii_whitespace() {
        i -= 1;
    }
    if i == 0 || b[i - 1] != b':' {
        return Err(format!("{pointer}: no colon before the value"));
    }
    i -= 1;
    while i > 0 && b[i - 1].is_ascii_whitespace() {
        i -= 1;
    }
    if i == 0 || b[i - 1] != b'"' {
        return Err(format!("{pointer}: no key before the value"));
    }
    i -= 1;
    while i > 0 && b[i - 1] != b'"' {
        i -= 1;
    }
    if i == 0 {
        return Err(format!("{pointer}: the key has no opening quote"));
    }
    let key_start = i - 1;
    let mut before = key_start;
    while before > 0 && b[before - 1].is_ascii_whitespace() {
        before -= 1;
    }
    if before > 0 && b[before - 1] == b',' {
        // A member before it: the comma that ends that one is the one to take.
        return Ok(edit::splice(text, before - 1..value.end, ""));
    }
    let mut after = value.end;
    while after < b.len() && b[after].is_ascii_whitespace() {
        after += 1;
    }
    if after < b.len() && b[after] == b',' {
        let mut next = after + 1;
        while next < b.len() && b[next].is_ascii_whitespace() {
            next += 1;
        }
        return Ok(edit::splice(text, key_start..next, ""));
    }
    Ok(edit::splice(text, key_start..value.end, ""))
}

/// `key` of the object at `object` set to the literal `rendered`, added when it is not there.
fn put_member(text: &str, object: &str, key: &str, rendered: &str) -> Result<String, String> {
    let pointer = format!("{object}/{}", edit::escape(key));
    if let Some(span) = edit::value_span(text, &pointer) {
        return Ok(edit::splice(text, span, rendered));
    }
    let span = edit::value_span(text, object).ok_or_else(|| format!("{object}: no such entry"))?;
    let quoted = serde_json::to_string(key).map_err(|e| e.to_string())?;
    insert_member(text, span, &format!("{quoted}: {rendered}"))
}

/// The entry `index` of `target`'s list `key`: where it is in the file, and the object.
fn entry_in(
    root: &serde_json::Value,
    target: &Target,
    key: &str,
    index: usize,
) -> Result<(String, serde_json::Value), String> {
    let (i, d) = edit::domain_at(root, &target.name)?;
    let entry = d
        .get(key)
        .and_then(|l| l.as_array())
        .and_then(|l| l.get(index))
        .ok_or_else(|| format!("{}.{key}[{index}]: no such entry", target.name))?;
    Ok((
        format!("/domains/{i}/{}/{index}", edit::escape(key)),
        entry.clone(),
    ))
}

/// A patch as a person reads it: the cells it runs over in each of the face's axes.
fn spans(face: FaceSpec, patch: Patch) -> String {
    let names = ["i", "j", "k"];
    let [a, b] = face.axes();
    format!(
        "{} {}..{}, {} {}..{}",
        names[a], patch.0[0], patch.1[0], names[b], patch.0[1], patch.1[1]
    )
}

/// Point entry `index` of `target`'s list `key` at `face`, and at `patch` on it when given.
///
/// What this writes is what the format reads:
///
/// - `face`, always;
/// - `from` and `to` for a patch, and neither for the whole face — a patch a list does not take is
///   dropped and said, since a held face is the whole face;
/// - `area_cm2`, when the entry states one, as the solid cells it now covers times a cell's face.
///   An area was stated for the face it named, and is the wrong number for any other.
///
/// `Ok(None)` when the entry already says exactly this, so pointing again at what is set changes
/// nothing — not even an area a person stated by hand. `Ok(Some((text, said)))` otherwise, `said`
/// being the sentence for the status line. Refused, and the file left alone, when the face or
/// patch holds no solid cell: that cools nothing, and a scene that says otherwise is the silence
/// every refusal in the format's cooling exists to stop.
pub fn set_face(
    text: &str,
    target: &Target,
    key: &str,
    index: usize,
    face: FaceSpec,
    patch: Option<Patch>,
) -> Result<Option<(String, String)>, String> {
    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    let placing = templates::placing(&target.kind, key);
    if !matches!(placing, Some(Placing::Face | Placing::FacePatch)) {
        return Err(format!(
            "{}.{key} does not say where by a face",
            target.name
        ));
    }
    let (pointer, entry) = entry_in(&root, target, key, index)?;
    let site = format!("{} {key} {}", target.name, index + 1);
    let mut dropped = false;
    let patch = match (placing, patch) {
        (Some(Placing::FacePatch), Some(p)) if p == target.whole(face) => None,
        (Some(Placing::FacePatch), p) => p,
        (_, p) => {
            dropped = p.is_some();
            None
        }
    };
    let now = match stated(text, &format!("/domains/{}/{key}/{index}", target.name)) {
        Some(Stated::Face { face, patch }) => Some((face, patch)),
        _ => None,
    };
    if now == Some((face, patch)) {
        return Ok(None);
    }

    let mut out = edit::set_text(text, &format!("{pointer}/face"), face.spelled())?;
    let mut said = format!("{site}: {}", face.spelled());
    match patch {
        Some(p) => {
            said.push_str(&format!(", cells {}", spans(face, p)));
            let (from, to) = (
                format!("[{}, {}]", p.0[0], p.0[1]),
                format!("[{}, {}]", p.1[0], p.1[1]),
            );
            let has = |key: &str| edit::value_span(&out, &format!("{pointer}/{key}")).is_some();
            if !has("from") && !has("to") {
                // **Together, as one member's worth of text**, which is how a patch is written —
                // scene 29's pad states both on one line — so taking the patch off again, which
                // cuts each with its comma, gives back the bytes it was given.
                let span = edit::value_span(&out, &pointer)
                    .ok_or_else(|| format!("{pointer}: no span"))?;
                out = insert_member(&out, span, &format!(r#""from": {from}, "to": {to}"#))?;
            } else {
                out = put_member(&out, &pointer, "from", &from)?;
                out = put_member(&out, &pointer, "to", &to)?;
            }
        }
        None => {
            out = cut_member(&out, &pointer, "to")?;
            out = cut_member(&out, &pointer, "from")?;
        }
    }
    if let Some(was) = entry.get("area_cm2").and_then(|a| a.as_f64()) {
        if target.void.is_none() {
            return Err(format!(
                "{site}: {} does not build yet, so which cells of {} are solid is not known and \
                 its area cannot be written — fix the scene, then point at the face again",
                target.name,
                face.spelled()
            ));
        }
        let (area, solid) = target.area_cm2(face, patch);
        if solid == 0 {
            return Err(format!(
                "{site}: {} has no solid cell of {} to cool{}",
                face.spelled(),
                target.name,
                if patch.is_some() {
                    " under that box"
                } else {
                    " — the part does not reach it"
                }
            ));
        }
        let span = edit::value_span(&out, &format!("{pointer}/area_cm2"))
            .ok_or_else(|| format!("{pointer}/area_cm2: no span"))?;
        out = edit::splice(&out, span, &edit::as_float(area));
        let of = if patch.is_some() {
            "solid cells under it"
        } else {
            "solid cells on it"
        };
        said.push_str(&format!(
            ", {} cm² over {solid} {of} (was {})",
            edit::as_float(area),
            edit::as_float(was)
        ));
    }
    if dropped {
        said.push_str(&format!(
            " — a {key} face is a whole face, so the box was not kept"
        ));
    }
    Ok(Some((out, said)))
}

/// What a cooled entry's stated area is beside the area its solid cells give, when the two differ:
/// the sentence to add when pointing again at what is set writes nothing, and the one the entry's
/// form shows.
///
/// `None` when they agree, when the entry states no area or no face, and when the solid cells are
/// not known. A stated area that differs may be meant — scene 29's bracket states the 16.5 cm² of
/// its STL's top against 16.4 rasterised — so this says so and changes nothing. A face with no
/// solid cell on it never gets here: the format refuses it, so the scene does not build and the
/// void map is not known.
pub fn area_note(text: &str, target: &Target, key: &str, index: usize) -> Option<String> {
    let root: serde_json::Value = serde_json::from_str(text).ok()?;
    let (_, entry) = entry_in(&root, target, key, index).ok()?;
    let stated_area = entry.get("area_cm2")?.as_f64()?;
    target.void.as_ref()?;
    let Some(Stated::Face { face, patch }) =
        stated(text, &format!("/domains/{}/{key}/{index}", target.name))
    else {
        return None;
    };
    let (area, solid) = target.area_cm2(face, patch);
    ((area - stated_area).abs() > 1e-9 * area.abs().max(stated_area.abs())).then(|| {
        format!(
            "; its area_cm2 is {} where its {solid} solid cells make {}",
            edit::as_float(stated_area),
            edit::as_float(area)
        )
    })
}

/// Set the box of entry `index` of `target`'s list `key` to `patch` on `face`, in the face's two
/// axes, and leave its third axis as it was.
///
/// A drag on a face says where in two axes and nothing about the third, so the depth the box
/// already had is the depth it keeps — a heat source drawn on the top of a die stays the layers it
/// was. `Ok(None)` when the box already spans exactly this.
pub fn set_cells(
    text: &str,
    target: &Target,
    key: &str,
    index: usize,
    face: FaceSpec,
    patch: Patch,
) -> Result<Option<(String, String)>, String> {
    let root: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("this scene does not parse: {e}"))?;
    if templates::placing(&target.kind, key) != Some(Placing::Cells) {
        return Err(format!("{}.{key} is not a box of cells", target.name));
    }
    let (pointer, _) = entry_in(&root, target, key, index)?;
    let Some(Stated::Cells { from, to }) =
        stated(text, &format!("/domains/{}/{key}/{index}", target.name))
    else {
        return Err(format!(
            "{} {key} {}: from and to are not three whole numbers each",
            target.name,
            index + 1
        ));
    };
    let [a, b] = face.axes();
    let (mut next_from, mut next_to) = (from, to);
    next_from[a] = patch.0[0];
    next_from[b] = patch.0[1];
    next_to[a] = patch.1[0];
    next_to[b] = patch.1[1];
    if (next_from, next_to) == (from, to) {
        return Ok(None);
    }
    let mut out = text.to_string();
    for (k, v) in [(a, patch.0[0]), (b, patch.0[1])] {
        out = edit::set_number(&out, &format!("{pointer}/from/{k}"), v as f64)?;
    }
    for (k, v) in [(a, patch.1[0]), (b, patch.1[1])] {
        out = edit::set_number(&out, &format!("{pointer}/to/{k}"), v as f64)?;
    }
    let (axis, _) = face.normal();
    let names = ["i", "j", "k"];
    Ok(Some((
        out,
        format!(
            "{} {key} {}: cells {} on {}, {} {}..{} as before",
            target.name,
            index + 1,
            spans(face, patch),
            face.spelled(),
            names[axis],
            from[axis],
            to[axis]
        ),
    )))
}
