//! The GUI-free half of the scene editor, for the same reason `viewer-core` is the GUI-free
//! half of the viewer: everything that could be got wrong the same way twice lives where a
//! test can reach it, and the shell is left holding a text box, a canvas and an event loop.
//!
//! What lives here:
//!
//! - **Checking** — the same two steps `pantometry check` runs, parse then build, with the
//!   parse error carried as `line:column` because that is what an editor puts a squiggle
//!   under.
//! - **Placement geometry** — every placed extent as eight posed corners, ready to wireframe,
//!   with the union bounds a camera fits to, and every designed part's triangles posed the same
//!   way. Both go through [`pantometry::core::Pose::point_to_world`], which was written here
//!   before any scene stated a pose and needed no change when one did — and they have to go
//!   through the *same* one, because a part drawn away from its own box is two pictures of one
//!   object with nothing to say which moved. See [`PlacedMesh`].
//! - **Editing** — the inspector's half of the loop: which values a selected row lets a person
//!   change, and the splices that change them in the text without touching any other byte —
//!   a number, a string, a domain added or removed, and a placement, which is the one that
//!   *creates* a key rather than replacing one. See [`edit`] for why they are splices and not a
//!   `Value` round trip.
//! - **The arithmetic behind a handle** — [`drag_along_axis`], which turns a drag on the screen
//!   into a distance along an axis in the world. Here rather than in the shell because it is a
//!   function of a camera and two vectors, and because its error has an *order* that a test can
//!   measure — which is how it was caught moving things 77% of the way.
//! - **Running and verifying** — thin passes over [`World::run`] and
//!   [`pantometry_world::verify::verify`], returning the run's JSON (which `viewer-core` reads)
//!   and the battery's rendered report.
//!
//! # The two halves, and which one this is
//!
//! ARCHITECTURE.md's platform rules split an editor into an authoring half, which is the
//! composition root and may name domains, and an inspection half, which must dispatch on the
//! shape of the data so that an eleventh physics costs no editor edit. This crate is the
//! authoring half's machinery: it consumes `Scene` and `DomainSpec` through `pantometry-world`'s
//! public API — `DomainSpec::placement()` is where domain knowledge already legitimately
//! lives — and hands the shell *shapes*: boxes, points, paths, readings. The shell's painting
//! code never sees a domain name except to label things, which is what keeps the viewport
//! open to domains that do not exist yet.

#![deny(missing_docs)]

pub mod edit;

pub use edit::{
    add_domain, domain_named, drag_along_axis, editable, pose_of, remove_domain, set_number,
    set_pose, set_text, set_turn, turn_of, Editable, Value,
};

/// A starting example of every domain the format defines, re-exported so a shell talks to this
/// crate and not past it. See [`pantometry_world::templates`] for how the list is kept in step.
pub use pantometry_world::templates::TEMPLATES;

use pantometry::units::LengthVec;
use pantometry_world::{DomainSpec, Parts, Scene, World};

/// Where a scene's `parts` come from, re-exported so a shell talks to this crate and not past it.
///
/// The native editor has a filesystem and uses [`OnDisk`]; the browser has uploads and uses
/// [`Uploaded`]. Both are the same scene format and the same builder — see [`Parts`].
pub use pantometry_world::{Beside, OnDisk, Uploaded};

/// The material names a scene can use without declaring them, re-exported for a shell's menus.
///
/// Asked for rather than copied: a list of substances typed into a dropdown is a list that drifts,
/// and the failure is a menu offering something the builder then refuses.
pub use pantometry_world::MATERIALS;

/// One placed extent, as the eight corners of its box in world coordinates, in metres.
///
/// Corner order is the binary one — bit 0 is x, bit 1 is y, bit 2 is z, low corner for a
/// clear bit — which is the same order `Camera::fit` walks and the order [`EDGES`] indexes.
#[derive(Clone, Debug, PartialEq)]
pub struct PlacedBox {
    /// The domain this box belongs to, for a label beside the wireframe.
    pub name: String,
    /// Eight corners, `[x, y, z]` each, in metres.
    pub corners: [[f64; 3]; 8],
}

/// One designed part's triangles, posed into world coordinates, in metres.
///
/// The wireframe beside it is the *grid* the part was rasterised onto; this is the part. A
/// staircase and the shape it approximates are different pictures and a person authoring an
/// assembly wants the second one, before there is a run and therefore before there is a field
/// to draw a surface of.
///
/// # Placed the way the rasteriser places it, or not at all
///
/// `Voxels::onto` reads the mesh's coordinates as written — an STL carries absolute positions —
/// against a grid whose origin is the block's local zero, and the *domain* is what a `Pose`
/// moves. So the only transform here is that same pose, and it has to be that one: a mesh drawn
/// a millimetre from its own voxels is two pictures of one part that disagree, which `mesh`'s
/// header calls worse than either being wrong, because nothing says which to believe.
#[derive(Clone, Debug, PartialEq)]
pub struct PlacedMesh {
    /// The domain the part belongs to. The same string [`PlacedBox::name`] carries, so hiding a
    /// domain hides its parts with it rather than leaving them floating over a hidden box.
    pub name: String,
    /// `name/parts[n]` — the site the build's notes and `verify`'s findings already use, so a
    /// reader can carry a rasterisation loss straight to the thing on screen.
    pub site: String,
    /// The file, for a label.
    pub stl: String,
    /// Where the domain's own coordinates sit in the world.
    ///
    /// **Stated, not baked**, which is the convention `Panel::place` established for a run and
    /// which this now shares. A rigid motion loses nothing when it is multiplied into a vertex
    /// list — unlike a box, where it does — so baking would not be *wrong* here. It would be a
    /// second convention in a workspace that has just spent four commits removing one.
    pub place: pantometry::scene::Placed,
    /// Triangles in the **domain's own** metres, three vertices each, wound as the file wound
    /// them. [`PlacedMesh::place`] is where they go.
    pub triangles: Vec<[[f64; 3]; 3]>,
}

/// The twelve edges of a box, as index pairs into [`PlacedBox::corners`].
///
/// Written out rather than generated because twelve constant pairs cannot be wrong quietly,
/// and a loop that generates them can — an edge between corners differing in two bits draws a
/// diagonal, and a wireframe with a diagonal in it still looks like geometry.
pub const EDGES: [(usize, usize); 12] = [
    (0, 1),
    (2, 3),
    (4, 5),
    (6, 7), // along x
    (0, 2),
    (1, 3),
    (4, 6),
    (5, 7), // along y
    (0, 4),
    (1, 5),
    (2, 6),
    (3, 7), // along z
];

/// What checking the editor's text produced. Every field is present on both success and
/// failure, because an editor redraws from whatever state it has: an error should not blank
/// the viewport that was drawn from the last good parse.
#[derive(Clone, Debug, Default)]
pub struct Checked {
    /// The parse or build error, `None` when the scene is runnable. A parse error leads with
    /// `line:column`, which is what an editor puts a squiggle under.
    pub error: Option<String>,
    /// One line of what the scene is — title, domains, duration, frames — for the header.
    pub summary: Option<String>,
    /// What the build dismissed and why — a stated condition a domain correctly ignored,
    /// with the measurement that earns it. Shown under the summary, because a dismissal
    /// nobody can see is the silence it exists to replace.
    pub notes: Vec<String>,
    /// Every placed extent, posed into world coordinates.
    pub boxes: Vec<PlacedBox>,
    /// Every designed part's triangles, posed the same way. Empty for a scene with no `parts`,
    /// and empty for one whose STLs could not be read — in which case `error` says so, because
    /// this is laid out from the same `files` the build used and fails for the same reasons.
    pub meshes: Vec<PlacedMesh>,
    /// The union of every box, `[x0, y0, z0, x1, y1, z1]`, for the camera to fit. `None` when
    /// nothing in the scene has geometry — which the shell must say rather than draw an empty
    /// viewport that reads as "framed on nothing".
    pub bounds: Option<[f64; 6]>,
}

/// Every designed part a scene names, as triangles in its domain's own frame.
///
/// # Why this is public, and why the STL is read twice
///
/// [`check`] draws these in the editor's viewport and the CLI hands them to `.gltf` and `.usda`,
/// which are two callers with one question. `World::build_with` already read **and** rasterised
/// every one of them, so a second parse is a small fraction of what a build costs, and the
/// alternative — carrying every mesh out of the build — would put a million triangles inside every
/// `World` for the sake of the callers that draw them.
///
/// # What the exporters do with it
///
/// They draw it beside the cells it became. That is the only picture that shows the rasterisation
/// loss instead of printing it: [`pantometry_world::Rasterised`] has reported a volume error per
/// part since parts existed, and a number in a terminal is not a thing a designer looks at.
///
/// It is drawn **uncoloured**, and deliberately so. The solver ran on the cells; a smooth surface
/// tinted from the field would claim a resolution the computation never had, which is the same
/// mistake as renormalising a colour scale per frame in a different dress. The editor already made
/// this decision — see the grey in its viewport, chosen to be unlike any scale in this workspace —
/// and this is the exporters agreeing with it.
///
/// A part whose file is missing or unreadable is **skipped silently here**, because this runs on a
/// scene that may not build; `World::build_with` is where a bad part is refused, and `check` calls
/// it first.
pub fn designed(scene: &Scene, files: &dyn Parts) -> Vec<PlacedMesh> {
    let placed = scene.placements();
    let mut out = Vec::new();
    for spec in &scene.domains {
        let DomainSpec::Block { parts, .. } = spec else {
            continue;
        };
        let placement = placed
            .get(spec.name())
            .copied()
            .unwrap_or_else(|| spec.placement());
        for (n, part) in parts.iter().enumerate() {
            let Ok(bytes) = files.bytes(&part.stl) else {
                continue;
            };
            let Ok(mesh) = pantometry::shape::Mesh::from_stl(&bytes) else {
                continue;
            };
            let triangles = mesh
                .triangles()
                .iter()
                .map(|t| {
                    // By component rather than by vertex: a `Triangle` holds `glam::DVec3`, which
                    // `pantometry` does not re-export, and adding a dependency to name a type in a
                    // closure signature is the wrong reason to add one.
                    [
                        [t.a.x, t.a.y, t.a.z],
                        [t.b.x, t.b.y, t.b.z],
                        [t.c.x, t.c.y, t.c.z],
                    ]
                })
                .collect();
            out.push(PlacedMesh {
                name: spec.name().to_string(),
                site: format!("{}/parts[{n}]", spec.name()),
                stl: part.stl.clone(),
                place: pantometry::scene::Placed::of(placement.pose),
                triangles,
            });
        }
    }
    out
}

/// A point in a domain's own frame, in the world's.
///
/// A thin pass to [`pantometry::scene::Placed::apply`], kept as a free function because the
/// editor's call sites read better without naming the type. It used to be the quaternion sandwich
/// written out again; `viewer-core` holds the only other copy, and that one is paid for -- the wire
/// format is the boundary between the two crates, which `the_wire_format_is_enough` explains.
pub fn place(at: pantometry::scene::Placed, p: [f64; 3]) -> [f64; 3] {
    at.apply(p)
}

/// Parse and build the text, and lay out its geometry.
///
/// The same two steps `pantometry check` runs, in the same order, so the editor and the
/// CLI cannot disagree about what a valid scene is. Geometry is laid out from the *parsed*
/// scene even when the build fails, because "the beam's face count disagrees with the bar's"
/// is exactly when a person wants to be looking at the boxes.
pub fn check(text: &str, files: &dyn Parts) -> Checked {
    let scene: Scene = match serde_json::from_str(text) {
        Ok(s) => s,
        Err(e) => {
            return Checked {
                error: Some(format!("{}:{}: {e}", e.line(), e.column())),
                ..Checked::default()
            }
        }
    };

    let (error, notes) = match World::build_with(scene.clone(), files) {
        Ok(world) => (None, world.notes().to_vec()),
        Err(e) => (Some(e), Vec::new()),
    };
    let mut out = Checked {
        error,
        notes,
        summary: Some(format!(
            "{}: {} domain(s), {:.3} s in {} frames",
            scene.title,
            scene.domains.len(),
            scene.duration_s,
            scene.frames
        )),
        ..Checked::default()
    };

    // The scene's own placements, so a pose the file states moves the box on screen — and so
    // this shell and `World` cannot disagree about where anything is.
    let placed = scene.placements();
    // Parts before the extents, because a part is geometry whether or not its domain states a box
    // to wireframe -- and because an assembly is the picture somebody authoring one is looking at.
    out.meshes = designed(&scene, files);
    for spec in &scene.domains {
        let placement = placed
            .get(spec.name())
            .copied()
            .unwrap_or_else(|| spec.placement());
        let Some(extent) = placement.extent else {
            continue;
        };
        let (lo, hi) = (extent.min.to_si(), extent.max.to_si());
        let mut corners = [[0.0; 3]; 8];
        for (i, corner) in corners.iter_mut().enumerate() {
            let local = LengthVec::m(
                if i & 1 == 0 { lo.x } else { hi.x },
                if i & 2 == 0 { lo.y } else { hi.y },
                if i & 4 == 0 { lo.z } else { hi.z },
            );
            let world = placement.pose.point_to_world(local).to_si();
            *corner = [world.x, world.y, world.z];
        }
        out.boxes.push(PlacedBox {
            name: spec.name().to_string(),
            corners,
        });
    }

    let mut bounds = [
        f64::INFINITY,
        f64::INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    ];
    for b in &out.boxes {
        for c in &b.corners {
            for a in 0..3 {
                bounds[a] = bounds[a].min(c[a]);
                bounds[a + 3] = bounds[a + 3].max(c[a]);
            }
        }
    }
    // Meshes too, so a domain that states parts and no extent is still framed. A part is inside
    // its grid by construction -- `Voxels::onto` refuses one that is not -- so this widens
    // nothing when both are present, and is the whole of the answer when only the parts are.
    for m in &out.meshes {
        for t in &m.triangles {
            for c in t {
                // In the domain's frame, so the placement has to be applied before a corner is a
                // corner of anything. It used to be baked in and this loop read it straight.
                let w = place(m.place, *c);
                for a in 0..3 {
                    bounds[a] = bounds[a].min(w[a]);
                    bounds[a + 3] = bounds[a + 3].max(w[a]);
                }
            }
        }
    }
    if bounds[0].is_finite() {
        out.bounds = Some(bounds);
    }
    out
}

/// Measure a ladder of grids for the uploaded parts, and say which row to use.
///
/// `pantometry_world::fit` with the shell's own framing: the names come from `files` rather than from a
/// scene, because this is the step *before* there is a scene — somebody has dropped CAD on the
/// window and does not yet know what `cells` and `cell_mm` to write.
pub fn fit(files: &Uploaded, budget_cells: usize, material: &str) -> Result<String, String> {
    let mut parts = Vec::new();
    for name in files.names() {
        let bytes = files.bytes(name)?;
        let mesh = pantometry::shape::Mesh::from_stl(&bytes).map_err(|e| format!("{name}: {e}"))?;
        parts.push((name.to_string(), mesh));
    }
    let fit = pantometry_world::fit::propose(&parts, budget_cells)?;
    let picked = fit.recommended(0.5);
    Ok(serde_json::json!({
        "table": fit.render(),
        "fragment": picked.map(|c| fit.scene_fragment(c, material)),
        "cell_mm": picked.map(|c| c.cell_m * 1e3),
        "cells": picked.map(|c| [c.counts.0, c.counts.1, c.counts.2]),
    })
    .to_string())
}

/// Run the scene and return the run as JSON — the same bytes `pantometry run scene.json out.json`
/// writes, which is the format `viewer-core` reads. A violation is the error, worded by the
/// kernel.
pub fn run(text: &str, files: &dyn Parts) -> Result<String, String> {
    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let title = scene.title.clone();
    let mut world = World::build_with(scene, files)?;
    let frames = world.run().map_err(|v| {
        format!(
            "the audit stopped the run at t = {:.4} s: {v}",
            world.time().to_si()
        )
    })?;
    Ok(pantometry::view::to_json(&title, &frames))
}

/// How a streamed run ended, when it did not fail.
///
/// Its own type rather than a bare `Ok(())`, because a stopped run and a finished one must
/// not be confusable: a partial run that reads as complete is a picture of something that did
/// not happen, which is this workspace's oldest failure shape wearing a scrub slider.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunEnd {
    /// Every frame the scene asked for was captured.
    Finished,
    /// The stop flag was raised between frames; what was emitted is a prefix, and the caller
    /// must say so wherever the frames are shown.
    Stopped,
}

/// Run the scene, emitting the run-so-far as JSON after every captured frame.
///
/// This is what makes a run **watchable while it happens**: each `emit` payload is a complete,
/// readable run — `viewer-core` parses every one — containing the frames captured so far, and
/// the last payload is byte-identical to what [`run`] returns for the same text, which the
/// tests pin. Intermediate payloads are the run *unsettled*: `settle_framing` runs once at the
/// end, exactly as [`World::run`] does, so a shell scrubbing mid-run sees each frame's own
/// framing and the picture settles when the run does — the honest rendering of a run that is
/// not finished yet.
///
/// `stop` is read between frames. A violation still emits nothing extra: the frames already
/// emitted stand, and the error carries the kernel's own words, so a shell can leave the
/// partial run on screen *beside* the reason it ended — which is precisely the view somebody
/// debugging a violation wants.
pub fn run_streaming(
    text: &str,
    files: &dyn Parts,
    stop: &std::sync::atomic::AtomicBool,
    mut emit: impl FnMut(String),
) -> Result<RunEnd, String> {
    use std::sync::atomic::Ordering;

    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let title = scene.title.clone();
    let mut world = World::build_with(scene.clone(), files)?;
    let dt = pantometry::units::Time::from_si(scene.duration_s / scene.frames as f64);
    let placed = world.placements();

    let mut frames = vec![pantometry::scene::capture(world.simulation(), &placed)];
    emit(pantometry::view::to_json(&title, &frames));
    for _ in 0..scene.frames {
        if stop.load(Ordering::Relaxed) {
            return Ok(RunEnd::Stopped);
        }
        world.advance(dt).map_err(|v| {
            format!(
                "the audit stopped the run at t = {:.4} s: {v}",
                world.time().to_si()
            )
        })?;
        frames.push(pantometry::scene::capture(world.simulation(), &placed));
        emit(pantometry::view::to_json(&title, &frames));
    }
    pantometry::scene::settle_framing(&mut frames);
    emit(pantometry::view::to_json(&title, &frames));
    Ok(RunEnd::Finished)
}

/// Run the verification battery and return its rendered report, with the findings count so
/// the shell can be as loud as the CLI's exit code.
pub fn verify(text: &str, deep: bool, files: &dyn Parts) -> Result<(String, usize), String> {
    let scene: Scene =
        serde_json::from_str(text).map_err(|e| format!("{}:{}: {e}", e.line(), e.column()))?;
    let battery = pantometry_world::verify::verify_with(&scene, deep, files)?;
    Ok((battery.render(), battery.findings.len()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const ROOM: &str = r#"{
  "title": "a room",
  "duration_s": 0.005,
  "frames": 2,
  "domains": [
    { "kind": "room", "name": "room", "width_m": 4.0, "height_m": 2.0,
      "cells_across": 21 }
  ]
}"#;

    /// The check mirrors `--check`: a valid scene has no error, a summary, and geometry.
    #[test]
    fn a_valid_scene_checks_clean_and_lays_out() {
        let c = check(ROOM, &OnDisk);
        assert!(c.error.is_none(), "{:?}", c.error);
        assert_eq!(c.boxes.len(), 1);
        let b = c.bounds.expect("a room has geometry");
        assert_eq!(b[3] - b[0], 4.0, "the box spans the stated width");
        assert_eq!(b[4] - b[1], 2.0);
    }

    /// A parse error carries `line:column`, which is what an editor puts a squiggle under —
    /// and the failure must not blank the rest of the state, it IS the state.
    #[test]
    fn a_parse_error_names_its_line_and_column() {
        let c = check("{ \"title\": ", &OnDisk);
        let e = c.error.expect("truncated JSON is an error");
        assert!(
            e.starts_with("1:"),
            "no position an editor could squiggle: {e}"
        );
        assert!(c.boxes.is_empty() && c.bounds.is_none());
    }

    /// A build error keeps the parsed geometry, because a cross-domain disagreement is
    /// exactly when a person wants to be looking at the boxes.
    #[test]
    fn a_build_error_keeps_the_geometry() {
        let two_names = ROOM.replace(
            "\"domains\": [",
            r#""domains": [
    { "kind": "room", "name": "room", "width_m": 1.0, "height_m": 1.0, "cells_across": 5 },"#,
        );
        let c = check(&two_names, &OnDisk);
        assert!(
            c.error
                .as_deref()
                .is_some_and(|e| e.contains("both called")),
            "{:?}",
            c.error
        );
        assert_eq!(c.boxes.len(), 2, "the boxes survive the refusal");
    }

    /// Every edge joins corners differing in exactly one bit — an edge, not a diagonal.
    #[test]
    fn the_edges_are_edges_and_cover_the_box() {
        for (a, b) in EDGES {
            assert_eq!((a ^ b).count_ones(), 1, "{a}-{b} is a diagonal");
        }
        // Each corner is touched by exactly three edges, as a cube's corners are.
        for corner in 0..8 {
            let touching = EDGES
                .iter()
                .filter(|(a, b)| *a == corner || *b == corner)
                .count();
            assert_eq!(touching, 3);
        }
    }

    /// The run's JSON is the wire format: `viewer-core` — which never links `pantometry` — reads
    /// it back whole. This is the editor standing on the same contract the viewer proved.
    #[test]
    fn a_run_round_trips_through_the_wire_format() {
        let json = run(ROOM, &OnDisk).expect("the room runs");
        let run = viewer_core::Run::from_json(&json).expect("the viewer's reader accepts it");
        assert_eq!(run.frames.len(), 3, "two frames plus the initial capture");
        assert!(run.frames[0].panels.iter().any(|p| p.name() == "room"));
    }

    /// The verify pass returns the same report the CLI prints, findings counted.
    #[test]
    fn verify_reports_and_counts_findings() {
        let (report, findings) = verify(ROOM, false, &OnDisk).expect("the battery runs");
        assert_eq!(findings, 0, "{report}");
        assert!(report.contains("determinism     two runs, identical bytes"));
    }

    /// Streaming emits one readable run per capture, and its last payload is byte-identical
    /// to [`run`]'s — the stream lands exactly where the batch does, so watching a run happen
    /// and reading it afterwards are the same run.
    #[test]
    fn a_streamed_run_is_watchable_and_lands_on_the_batch_answer() {
        let stop = std::sync::atomic::AtomicBool::new(false);
        let mut payloads = Vec::new();
        let end = run_streaming(ROOM, &OnDisk, &stop, |j| payloads.push(j)).expect("the room runs");
        assert_eq!(end, RunEnd::Finished);
        // The initial capture, one per frame, and the settled final.
        assert_eq!(payloads.len(), 4);
        for p in &payloads {
            viewer_core::Run::from_json(p).expect("every payload is a whole, readable run");
        }
        assert_eq!(
            payloads.last().unwrap(),
            &run(ROOM, &OnDisk).unwrap(),
            "the stream must land on the batch answer to the byte"
        );
    }

    /// The stop flag ends a run between frames, and the ending says Stopped — a prefix must
    /// never be confusable with the run.
    #[test]
    fn a_stopped_run_says_so() {
        use std::sync::atomic::{AtomicBool, Ordering};
        let stop = AtomicBool::new(false);
        let mut emitted = 0;
        let end = run_streaming(ROOM, &OnDisk, &stop, |_| {
            emitted += 1;
            stop.store(true, Ordering::Relaxed);
        })
        .expect("stopping is not a violation");
        assert_eq!(end, RunEnd::Stopped);
        assert_eq!(emitted, 1, "stopped after the first capture");
    }
}

/// One drawn cell of a field: where it is in the world, and what colour it is.
///
/// World space and no camera, deliberately. Two shells draw these now — the native window and
/// the browser — and the projection is `viewer-core`'s in both, so what lives here is the half
/// that could be got wrong the same way twice: **where the samples sit inside the placed box,
/// and what colour a value is.** The camera was already shared for that reason; this is the
/// same argument one layer along.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Splat {
    /// Where, in world metres.
    pub at: [f64; 3],
    /// Colour and coverage, straight sRGB with an alpha.
    pub rgba: [u8; 4],
}

/// What a field's splats are, and the sentence a reader needs beside them.
#[derive(Clone, Debug)]
pub struct Splatted {
    /// The cells to draw, in grid order. A shell projects, sorts by depth and paints.
    pub splats: Vec<Splat>,
    /// Cells skipped per axis. One means every cell is drawn.
    pub stride: usize,
    /// Whether the colours are Planck's or the conventional ramp — see [`field_splats`].
    pub physical: bool,
    /// Whether the conventional ramp used was the diverging one, so a shell's colour bar can say
    /// so and label its ends with the deflection rather than with the data's range.
    pub signed: bool,
    /// The note the canvas must carry, saying which of the two a reader is looking at.
    pub note: &'static str,
}

/// How a field's values become colours: the decision, made once and shared.
///
/// Two callers need it and must not disagree — [`field_splats`], which the flat painter used, and
/// [`field_shell`], which the shaded viewport uses. The decision is not a preference: whether a
/// field glows is a fact about its temperature, and a viewport that answered it differently from
/// the picture beside it would be showing two physics.
#[derive(Clone, Debug)]
pub struct Colouring {
    /// Whether the colours are Planck's rather than a ramp.
    physical: bool,
    /// The glow fraction of the hottest cell, which brightness is relative to.
    peak_glow: f64,
    /// How to read a value as kelvin, or `None` if the unit is not a temperature.
    kelvin_offset: Option<f64>,
    scale: Option<(f64, f64)>,
    signed: bool,
}

impl Colouring {
    /// Decide, from the unit and what the field actually holds.
    ///
    /// The dispatch is on the **unit**, which is data rather than a domain name, so a physics added
    /// next year is coloured correctly with no edit here.
    pub fn of(unit: &str, values: &[f64], scale: Option<(f64, f64)>) -> Colouring {
        let kelvin_offset = match unit {
            "K" => Some(0.0),
            "C" => Some(273.15),
            _ => None,
        };
        let hottest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let peak_glow = kelvin_offset
            .filter(|_| hottest.is_finite())
            .map_or(0.0, |o| pantometry::view::glow_fraction(hottest + o));
        Colouring {
            physical: peak_glow > 1e-6,
            peak_glow,
            kelvin_offset,
            scale,
            signed: scale_is_signed(scale),
        }
    }

    /// Whether these colours are Planck's. Public because the note a canvas carries depends on it.
    pub fn physical(&self) -> bool {
        self.physical
    }

    /// Whether the ramp in use is the diverging one, so a colour bar can label its ends right.
    pub fn signed(&self) -> bool {
        self.signed
    }

    /// A value's place in the run-wide range, `0..=1`. `0.5` when there is no range.
    pub fn place(&self, v: f64) -> f64 {
        match self.scale {
            Some((lo, hi)) if hi > lo => ((v - lo) / (hi - lo)).clamp(0.0, 1.0),
            _ => 0.5,
        }
    }

    /// The colour of a value, straight sRGB.
    pub fn srgb(&self, v: f64) -> [u8; 3] {
        if self.physical {
            let kelvin = v + self.kelvin_offset.unwrap_or(0.0);
            pantometry::view::blackbody_srgb(kelvin)
        } else {
            value_colour(v, self.scale)
        }
    }

    /// How bright a glowing cell is relative to this field's hottest, `0..=1`.
    ///
    /// A cool corner of a glowing block should be dark rather than merely bluer, which is what a
    /// photograph of it looks like. One for a field that is not glowing: a ramp carries its own
    /// lightness and dimming it would undo the property the ramp was built for.
    pub fn brightness(&self, v: f64) -> f64 {
        if !self.physical {
            return 1.0;
        }
        let kelvin = v + self.kelvin_offset.unwrap_or(0.0);
        (pantometry::view::glow_fraction(kelvin) / self.peak_glow.max(f64::MIN_POSITIVE))
            .clamp(0.0, 1.0)
            .sqrt()
    }

    /// The colour of a value in **linear** RGB, dimmed by [`Colouring::brightness`].
    ///
    /// Linear because a shader interpolates across a triangle and multiplies by a light, and both
    /// of those are wrong in sRGB — the glTF exporter shipped that defect and every export came out
    /// about 2.3× too bright in the midtones.
    pub fn linear(&self, v: f64) -> [f32; 3] {
        let srgb = self.srgb(v);
        let b = self.brightness(v) as f32;
        std::array::from_fn(|a| {
            let c = srgb[a] as f32 / 255.0;
            let linear = if c <= 0.04045 {
                c / 12.92
            } else {
                ((c + 0.055) / 1.055).powf(2.4)
            };
            linear * b
        })
    }

    /// The sentence a canvas must carry, saying which of the two colourings a reader is looking at.
    pub fn note(&self, stride: usize) -> &'static str {
        match (self.physical, stride) {
            (true, 1) => "field: colour is Planck's, not a palette",
            (true, _) => "field: Planck colour, subsampled — see the report for every cell",
            (false, 1) => "field: false colour — nothing here is hot enough to glow",
            (false, _) => "field: false colour, subsampled — see the report for every cell",
        }
    }
}

/// The most cells one field may draw in a frame.
///
/// A 100³ field is a million splats and a painter that stops painting. Past this the field is
/// subsampled with a stride and the note says so — a silently decimated picture is a picture of
/// a coarser simulation than the one that ran.
pub const MAX_SPLATS: usize = 8000;

/// Turn a field panel into cells to draw, coloured by physics where physics gives a colour.
///
/// # Two colourings, and the physics decides which
///
/// A temperature field is drawn in the colour a body at that temperature **actually is** —
/// Planck's law through the CIE matching functions, from [`pantometry::view::colour`] — whenever
/// anything in it is hot enough to emit visible light. That is not a palette: a melting block
/// glows the orange a melting block glows, and nothing here picked it.
///
/// Below that, physics gives no colour. A body at 300 K emits nothing visible and this
/// workspace holds no visible reflectance for it to have instead, so the field falls back to a
/// conventional ramp — which says *more* and *less* and does not pretend to say *looks like* —
/// and [`Splatted::note`] states which, because a false colour mistaken for a real one is a
/// wrong answer that looks right.
///
/// The dispatch is on the panel's **unit**, which is data rather than a domain name, so a
/// physics added next year is coloured correctly with no edit here.
pub fn field_splats(
    corners: &[[f64; 3]; 8],
    counts: (usize, usize, usize),
    values: &[f64],
    unit: &str,
    scale: Option<(f64, f64)>,
) -> Splatted {
    let (nx, ny, nz) = counts;
    if nx == 0 || ny == 0 || nz == 0 || values.len() < nx * ny * nz {
        return Splatted {
            splats: Vec::new(),
            stride: 1,
            physical: false,
            signed: false,
            note: "field: the panel's values do not fill its grid — not drawn",
        };
    }

    let colouring = Colouring::of(unit, values, scale);
    let physical = colouring.physical();

    let total = nx * ny * nz;
    let stride = if total > MAX_SPLATS {
        ((total as f64 / MAX_SPLATS as f64).cbrt().ceil() as usize).max(1)
    } else {
        1
    };

    // **Where each sample is, from the one place that knows.** This function carried the walk —
    // corner 0 as origin, `i / (n - 1)` along each axis, non-finite skipped — and then the viewer
    // needed the same thing to draw a field at all. Two copies of *where is sample `i`* are two
    // pictures of one run that can disagree about it, so there is one, in the crate that owns the
    // panel: [`viewer_core::field_points`].
    let mut splats = Vec::new();
    for (at, v) in viewer_core::field_points(corners, counts, values, stride) {
        // One scale across the whole run, never per frame, for the reason
        // `viewer-core` states.
        let s = colouring.place(v);
        let [r, g, b] = colouring.srgb(v);
        let rgba = if physical {
            [r, g, b, ((colouring.brightness(v) * 235.0) as u8).max(6)]
        } else {
            // Opacity climbs with the value's place in the range, so the quiet bulk
            // of a field clears out of the way. The colour is not computed from
            // that place: a signed field's colour is placed about **zero** and its
            // opacity about the middle of the deflection, and those are not the same
            // point unless the range happens to be symmetric.
            let a = if colouring.signed() {
                (2.0 * s - 1.0).abs()
            } else {
                s
            };
            [r, g, b, (30.0 + 200.0 * a * a) as u8]
        };
        splats.push(Splat { at, rgba });
    }

    Splatted {
        splats,
        stride,
        physical,
        signed: colouring.signed(),
        note: colouring.note(stride),
    }
}

/// A field's boundary as a shaded surface, in world metres, ready for a GPU.
///
/// # Why the boundary and not the cells
///
/// The flat painter drew every cell as a translucent circle sorted far to near. It composites
/// correctly and it reads as a point cloud, which is the criticism the HTML report's first draft
/// earned and the difference between this viewport and the one in usdview. A solid has a surface,
/// and one quad per cell face whose neighbour is absent is that surface: on a solid 9 cubed that is
/// 486 quads against 4374, so 89% of the work was hidden anyway.
///
/// # The geometry is the library's
///
/// [`pantometry::view::mesh::field_surface`] builds it -- the same function that writes glTF and
/// USD, so this viewport cannot disagree with an export about how big a solid is. That arithmetic
/// has been wrong once: a 40 mm cube exported 80 mm across, because a field sampled corner to
/// corner means the end node owns half a cell. A second copy here would be a second chance at it.
///
/// The box may be **placed**, so it is not axis-aligned in general. The mesh is built in the unit
/// cube and mapped through the box's own three axes; the normals are mapped as the cross products
/// of the axis pairs, which is exact for a rotated and non-uniformly scaled box where transforming
/// the normal by the same matrix would not be.
/// # Two surfaces, one mapping
///
/// `level` chooses **which** surface, and nothing else changes. `None` is the boundary of the
/// cells that hold a value — the outside of the object. `Some(v)` is the isosurface where the
/// field reaches `v`, which is a different question: not *where is the block* but *where is it
/// 100 °C*.
///
/// One function rather than two because everything after the mesh is the same, and the part that
/// is the same is the part that has been wrong before: mapping the unit cube through a placed
/// box's own axes, with normals as cross products of the axis pairs. A second copy of that would
/// be a second chance at the 40 mm cube that exported 80 mm across.
pub fn field_shell(
    corners: &[[f64; 3]; 8],
    counts: (usize, usize, usize),
    values: &[f64],
    unit: &str,
    scale: Option<(f64, f64)>,
    level: Option<f64>,
) -> Shelled {
    let colouring = Colouring::of(unit, values, scale);
    let unit_cube = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0];
    let mesh = match level {
        Some(v) => pantometry::view::mesh::isosurface(counts, unit_cube, values, v),
        None => pantometry::view::mesh::field_surface(counts, unit_cube, values),
    };
    if mesh.indices.is_empty() {
        return Shelled {
            positions: Vec::new(),
            normals: Vec::new(),
            colours: Vec::new(),
            source: Vec::new(),
            indices: Vec::new(),
            stride: mesh.stride,
            // A field of one dimension is a graph, not geometry -- `field_surface` says so by
            // returning nothing, and saying nothing here would leave a reader with an empty
            // viewport and no reason for it.
            // **The same condition `field_surface` uses**, not a paraphrase of it. The first
            // spelling here asked whether *any* axis had more than one cell, which a 32-by-1-by-1
            // row of samples satisfies -- so a graph was reported as an empty solid.
            note: if [counts.0, counts.1, counts.2]
                .iter()
                .filter(|&&n| n > 1)
                .count()
                < 2
            {
                "field: one dimension -- a graph, not geometry"
            } else {
                "field: nothing present to draw a surface of"
            },
            physical: colouring.physical(),
            signed: colouring.signed(),
        };
    }

    // The box's origin and three edge vectors. Corner 0 is the low one and bits 0, 1, 2 step one
    // axis each, which is the order `EDGES` is written against.
    let o = corners[0];
    let axis = |c: [f64; 3]| [c[0] - o[0], c[1] - o[1], c[2] - o[2]];
    let e = [axis(corners[1]), axis(corners[2]), axis(corners[4])];

    // The world normal of each local face direction. `+x` faces along `e1 x e2`, and so round,
    // with the sign following the local direction.
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let unit_of = |v: [f64; 3]| {
        let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        if n > 0.0 {
            [(v[0] / n) as f32, (v[1] / n) as f32, (v[2] / n) as f32]
        } else {
            [0.0, 0.0, 1.0]
        }
    };
    let faces = [
        unit_of(cross(e[1], e[2])),
        unit_of(cross(e[2], e[0])),
        unit_of(cross(e[0], e[1])),
    ];

    let mut out = Shelled {
        positions: Vec::with_capacity(mesh.positions.len()),
        normals: Vec::with_capacity(mesh.positions.len()),
        colours: Vec::with_capacity(mesh.positions.len()),
        source: mesh.source.clone(),
        indices: mesh.indices.clone(),
        stride: mesh.stride,
        note: colouring.note(mesh.stride),
        physical: colouring.physical(),
        signed: colouring.signed(),
    };
    for (i, p) in mesh.positions.iter().enumerate() {
        let (u, v, w) = (p[0] as f64, p[1] as f64, p[2] as f64);
        out.positions.push([
            o[0] + e[0][0] * u + e[1][0] * v + e[2][0] * w,
            o[1] + e[0][1] * u + e[1][1] * v + e[2][1] * w,
            o[2] + e[0][2] * u + e[1][2] * v + e[2][2] * w,
        ]);

        // The local normal is one of the six axis directions, so the axis it names is the one
        // whose component is largest, and the sign carries straight through.
        let n = mesh.normals.get(i).copied().unwrap_or([0.0, 0.0, 1.0]);
        let a = (0..3)
            .max_by(|x, y| n[*x].abs().total_cmp(&n[*y].abs()))
            .unwrap_or(2);
        let s = if n[a] < 0.0 { -1.0 } else { 1.0 };
        out.normals
            .push([faces[a][0] * s, faces[a][1] * s, faces[a][2] * s]);

        let value = mesh
            .source
            .get(i)
            .and_then(|c| values.get(*c as usize))
            .copied()
            .unwrap_or(f64::NAN);
        out.colours.push(if value.is_finite() {
            colouring.linear(value)
        } else {
            // A vertex whose cell has no value is not a vertex `field_surface` emits -- a cell has
            // to be present to have a face. Mid-ramp rather than black, so if that ever changes
            // the picture says something instead of drawing a hole.
            colouring.linear(0.0)
        });
    }
    out
}

/// The most surface vertices one field offers the cursor.
///
/// A readout has to walk its candidates, and a 128 cubed boundary is 200 000 faces and 800 000
/// vertices — a walk per pointer position at sixty frames a second. The old splat path had the same
/// bound for the same reason and called it `MAX_SPLATS`; this is the same idea about a surface, and
/// the readout names the nearest *sampled* vertex rather than pretending to the nearest cell.
pub const MAX_PROBES: usize = 3000;

impl Shelled {
    /// Vertices to offer the cursor, as `(world position, index into the panel's values)`.
    ///
    /// Strided, so a fine field costs the same as a coarse one. Deduplicating by cell would be the
    /// nicer answer and is not free either: a flat-shaded quad has four vertices of one cell, so the
    /// stride is four times coarser than it looks and that is stated rather than hidden.
    pub fn probes(&self) -> Vec<([f64; 3], u32)> {
        if self.positions.is_empty() {
            return Vec::new();
        }
        let stride = self.positions.len().div_ceil(MAX_PROBES).max(1);
        self.positions
            .iter()
            .zip(&self.source)
            .step_by(stride)
            .map(|(p, c)| (*p, *c))
            .collect()
    }
}

/// A field's boundary, mapped into the world.
#[derive(Clone, Debug)]
pub struct Shelled {
    /// Vertex positions in world metres.
    pub positions: Vec<[f64; 3]>,
    /// Unit world normals, one per vertex.
    pub normals: Vec<[f32; 3]>,
    /// Linear RGB, one per vertex, already dimmed by the glow where the colour is Planck's.
    pub colours: Vec<[f32; 3]>,
    /// For each vertex, the index into the panel's values that coloured it -- for a readout.
    pub source: Vec<u32>,
    /// Triangles, three indices to a face.
    pub indices: Vec<u32>,
    /// Cells skipped per axis. One means every cell's faces are drawn.
    pub stride: usize,
    /// The note the canvas must carry.
    pub note: &'static str,
    /// Whether the colours are Planck's.
    pub physical: bool,
    /// Whether the ramp is the diverging one.
    pub signed: bool,
}

/// The colour of a value on a run-wide scale, for a quantity physics gives no colour to.
///
/// **`pantometry::view::ramp`, not a formula written here.** What was here was a straight line in
/// sRGB from blue to red, and both shells kept their own copy of it, so the workspace held four
/// spellings of "cool to warm". Measured over 256 steps, that line covered **16.6 L\*** — 44.1 to
/// 60.7 — against the 74 the library's scale covers. Lightness is where a scale survives a
/// greyscale print, a projector, and the eight per cent of men with a colour-vision deficiency,
/// and a scale carrying its whole signal on the blue-red axis is the one pair those readers
/// cannot separate. Seventy-five of its 255 steps also ran backwards.
///
/// The library's scale is built in CIE LCh with lightness linear in the value, and its properties
/// are pinned by tests there rather than asserted here.
///
/// **Which scale is chosen by the numbers.** A range that straddles zero is a signed quantity —
/// a pressure, a component of **E** — and gets the diverging scale with its neutral at the value
/// **zero**, not at the middle of the range: for -100 to +300 those are a quarter of the scale
/// apart. Anything else gets the sequential one. Nothing here reads the domain's name.
pub fn value_colour(value: f64, scale: Option<(f64, f64)>) -> [u8; 3] {
    let Some((lo, hi)) = scale.filter(|(lo, hi)| hi > lo) else {
        // A panel holding one constant has no scale to place a value on. Mid-ramp rather than
        // invisible: the cell is there and the picture should say so.
        return pantometry::view::ramp::sequential(0.5);
    };
    if pantometry::view::ramp::is_signed(lo, hi) {
        let reach = lo.abs().max(hi.abs()).max(f64::MIN_POSITIVE);
        pantometry::view::ramp::diverging(0.5 + 0.5 * value / reach)
    } else {
        pantometry::view::ramp::sequential((value - lo) / (hi - lo))
    }
}

/// Whether a run-wide scale is a swing about zero, so a shell can say which scale its bar shows.
pub fn scale_is_signed(scale: Option<(f64, f64)>) -> bool {
    scale.is_some_and(|(lo, hi)| hi > lo && pantometry::view::ramp::is_signed(lo, hi))
}

/// The value at a position `u` in `0..=1` along the colour bar, inverting [`value_colour`].
///
/// A bar a reader cannot read a number off is a decoration. For a diverging scale this spans
/// `-reach..=+reach` rather than `lo..=hi`, because that is the range the scale actually covers:
/// it is symmetric by construction and labelling it with the data's range would mislabel it.
pub fn bar_value(u: f64, scale: Option<(f64, f64)>) -> f64 {
    let Some((lo, hi)) = scale.filter(|(lo, hi)| hi > lo) else {
        return 0.0;
    };
    if pantometry::view::ramp::is_signed(lo, hi) {
        let reach = lo.abs().max(hi.abs());
        (2.0 * u - 1.0) * reach
    } else {
        lo + u * (hi - lo)
    }
}

/// Indices of `depths`, furthest from the eye first.
///
/// **Painter's algorithm needs the depth the projection computed, and one shell had one that did
/// not.** `Camera::project` returns a `depth` — "distance from the eye, larger is further" — and
/// the native shell discarded it, then recovered a stand-in by projecting each point a second
/// time one millimetre along **world z** and taking the reciprocal of how far the two landed
/// apart on screen.
///
/// That separation is proportional to how far off the view axis a point is, not to how far away
/// it is. On a 6x6x6 field of splats it ordered **6,026 of 23,220 pairs backwards — 26%** — so a
/// quarter of a translucent volume composited in the wrong order, and worst at the centre of the
/// screen, which is where the object is. The browser shell sorted by the real depth and was
/// right; the native one kept its own copy and was not. This is that copy, written once.
///
/// A `NaN` sorts to the far end rather than panicking: `sort_by` over a partial order containing
/// one is a panic inside a paint loop, and a splat with no depth is a splat to draw first.
pub fn far_to_near(depths: &[f64]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..depths.len()).collect();
    order.sort_by(|&a, &b| {
        let (x, y) = (depths[a], depths[b]);
        y.partial_cmp(&x)
            .unwrap_or_else(|| y.is_nan().cmp(&x.is_nan()))
    });
    order
}

/// A number with its magnitude intact, for a panel of readings.
///
/// **Significant figures, not decimal places.** The native shell printed readings `{:.4}`, so a
/// cavity holding 3.19e-10 J showed `0.0000` beside a field the same run reported at 921 V/m —
/// the same erasure `readings_csv` was fixed for, and then `to_json`, and now here in the third
/// writer. A scene format spanning nanoseconds to hours and femtojoules to kilowatts has no one
/// scale at which a fixed-point format is right.
pub fn magnitude(v: f64) -> String {
    if !v.is_finite() {
        return "-".to_string();
    }
    if v == 0.0 {
        return "0".to_string();
    }
    let a = v.abs();
    if !(1e-3..1e6).contains(&a) {
        format!("{v:.4e}")
    } else if a >= 100.0 {
        format!("{v:.2}")
    } else if a >= 1.0 {
        format!("{v:.4}")
    } else {
        format!("{v:.6}")
    }
}

/// What kind of thing a node in the outliner is.
///
/// Dispatched on the **shape of the data**, never on a domain's name — the same rule the viewport
/// paints by, and the reason an eleventh physics costs no edit here. `Placed` is the exception
/// that proves it: a box a scene declared but no run has filled yet is a shape too.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeKind {
    /// The scene itself.
    Scene,
    /// A grouping row with no geometry of its own.
    Group,
    /// A placed extent from the scene text, before or without a run.
    Placed,
    /// A field sampled on a grid.
    Field,
    /// Bodies at positions.
    Points,
    /// Runs of connected points.
    Paths,
    /// A named scalar a domain reports.
    Reading,
}

impl NodeKind {
    /// A short word for the type column, the way an outliner labels one.
    pub fn label(self) -> &'static str {
        match self {
            NodeKind::Scene => "Scene",
            NodeKind::Group => "Group",
            NodeKind::Placed => "Extent",
            NodeKind::Field => "Field",
            NodeKind::Points => "Bodies",
            NodeKind::Paths => "Paths",
            NodeKind::Reading => "Scalar",
        }
    }
}

/// One row of the outliner.
///
/// Flat with parent indices rather than nested, because that is what an immediate-mode shell can
/// walk without allocating a closure per level, and because a stable index is the only selection
/// handle that survives a rebuild — the tree is rebuilt on every check and every streamed frame,
/// and a selection keyed by pointer or by row number would move under the reader's hand.
///
/// The handle is [`Node::path`]: a slash-joined name like `/room/pressure`, stable across
/// rebuilds for as long as the scene calls the thing that.
#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    /// The stable handle, `/`-joined. What a selection is stored as.
    pub path: String,
    /// The name shown in the row.
    pub name: String,
    /// What it is.
    pub kind: NodeKind,
    /// Index of the parent row, or `None` for the root.
    pub parent: Option<usize>,
    /// How deep to indent.
    pub depth: usize,
    /// The box it occupies in world metres, if it has one.
    pub bounds: Option<[f64; 6]>,
    /// What its values are in, for the inspector's header.
    pub unit: String,
    /// The lines the inspector shows for this row, in order: label and value.
    ///
    /// Built here rather than in the shell because they are the answer to "what is this", and
    /// two shells asking that question should not get two answers.
    pub detail: Vec<(String, String)>,
}

/// The outliner, and what the inspector says about each row.
///
/// Built from the scene's placed extents and, when there is one, the run's panels and readings.
/// A scene that has been checked but not run still has a tree: the boxes are there, and being
/// able to select and inspect one before anything is computed is most of what an outliner is for.
#[derive(Clone, Debug, Default)]
pub struct Tree {
    /// The rows, parents before children, in the order an outliner draws them.
    pub nodes: Vec<Node>,
}

impl Tree {
    /// The row with this path, if it is still in the tree.
    pub fn find(&self, path: &str) -> Option<&Node> {
        self.nodes.iter().find(|n| n.path == path)
    }

    /// Whether `path` is `ancestor` or sits under it, by path prefix.
    ///
    /// A prefix test on `/`-joined names, with the separator required — so `/room` contains
    /// `/room/pressure` and does **not** contain `/roomier`, which a bare `starts_with` would
    /// have said it did.
    pub fn contains(ancestor: &str, path: &str) -> bool {
        path == ancestor
            || (path.starts_with(ancestor) && path.as_bytes().get(ancestor.len()) == Some(&b'/'))
    }
}

fn fmt_box(b: [f64; 6]) -> String {
    format!(
        "({:.4}, {:.4}, {:.4}) .. ({:.4}, {:.4}, {:.4})",
        b[0], b[1], b[2], b[3], b[4], b[5]
    )
}

fn fmt_span(b: [f64; 6]) -> String {
    format!(
        "{} x {} x {}",
        metres(b[3] - b[0]),
        metres(b[4] - b[1]),
        metres(b[5] - b[2])
    )
}

/// A length in metres, on the unit an engineer would say it in.
///
/// Trailing zeros count as trailing only **after a decimal point**: stripping them without that
/// check turns 500 into 5, which captioned a 500 um channel as a 5 um one in the HTML report
/// before it was found.
pub fn metres(m: f64) -> String {
    let a = m.abs();
    let body = |x: f64| {
        let t = format!("{x:.3}");
        if t.contains('.') {
            t.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            t
        }
    };
    if !a.is_finite() {
        "-".to_string()
    } else if a >= 1e4 {
        format!("{m:.3e} m")
    } else if a >= 1.0 || a == 0.0 {
        format!("{} m", body(m))
    } else if a >= 1e-3 {
        format!("{} mm", body(m * 1e3))
    } else {
        format!("{} um", body(m * 1e6))
    }
}

/// Build the outliner from a checked scene and, if there is one, a run.
///
/// A scene that has been checked but never run still has a tree — the placed extents are there,
/// and selecting one to read its size before anything is computed is most of what an outliner is
/// for. A run adds what it actually produced, which is not the same list: a domain can be placed
/// and draw nothing, and a domain with no placement at all can still report scalars.
///
/// `frame` picks which frame's numbers the inspector shows. The **range** it shows is the run's,
/// not the frame's, for the reason the colour scale is: a per-frame range makes a decaying mode
/// look constant.
pub fn tree(checked: &Checked, run: Option<&viewer_core::Run>, frame: usize) -> Tree {
    let mut nodes: Vec<Node> = Vec::new();
    let mut push = |n: Node| -> usize {
        nodes.push(n);
        nodes.len() - 1
    };

    let root = push(Node {
        path: "/".into(),
        name: checked
            .summary
            .as_deref()
            .map(|s| s.split(" \u{2014} ").next().unwrap_or(s).to_string())
            .unwrap_or_else(|| "scene".into()),
        kind: NodeKind::Scene,
        parent: None,
        depth: 0,
        bounds: checked.bounds,
        unit: String::new(),
        detail: {
            let mut d = vec![
                ("domains".into(), checked.boxes.len().to_string()),
                (
                    "frames".into(),
                    run.map_or("-".into(), |r| r.frames.len().to_string()),
                ),
            ];
            if let Some(b) = checked.bounds {
                d.push(("bounds".into(), fmt_box(b)));
                d.push(("size".into(), fmt_span(b)));
            }
            if let Some(e) = &checked.error {
                d.push(("error".into(), e.clone()));
            }
            for (i, note) in checked.notes.iter().enumerate() {
                d.push((format!("note {}", i + 1), note.clone()));
            }
            d
        },
    });

    // Every domain the scene placed, whether or not a run has filled it. Grouped under one row so
    // a scene with thirty extents does not bury the run's output.
    if !checked.boxes.is_empty() {
        let group = push(Node {
            path: "/extents".into(),
            name: "Placed extents".into(),
            kind: NodeKind::Group,
            parent: Some(root),
            depth: 1,
            bounds: checked.bounds,
            unit: String::new(),
            detail: vec![("count".into(), checked.boxes.len().to_string())],
        });
        for b in &checked.boxes {
            let bb = box_of(&b.corners);
            push(Node {
                path: format!("/extents/{}", b.name),
                name: b.name.clone(),
                kind: NodeKind::Placed,
                parent: Some(group),
                depth: 2,
                bounds: Some(bb),
                unit: String::new(),
                detail: vec![
                    ("bounds".into(), fmt_box(bb)),
                    ("size".into(), fmt_span(bb)),
                    ("origin".into(), fmt_point(b.corners[0])),
                ],
            });
        }
    }

    let Some(run) = run else {
        return Tree { nodes };
    };
    let Some(f) = run
        .frames
        .get(frame.min(run.frames.len().saturating_sub(1)))
    else {
        return Tree { nodes };
    };

    if !f.panels.is_empty() {
        let group = push(Node {
            path: "/run".into(),
            name: "Run".into(),
            kind: NodeKind::Group,
            parent: Some(root),
            depth: 1,
            bounds: None,
            unit: String::new(),
            detail: vec![
                ("frames".into(), run.frames.len().to_string()),
                (
                    "frame".into(),
                    format!("{} of {}", frame + 1, run.frames.len()),
                ),
                ("t".into(), format!("{} s", magnitude(f.t))),
            ],
        });
        for panel in &f.panels {
            let scale = run.scale_of(panel.name());
            let bounds = panel.bounds();
            let mut detail: Vec<(String, String)> = vec![
                ("unit".into(), panel.unit().to_string()),
                ("samples".into(), panel.values().len().to_string()),
            ];
            let kind = match panel {
                viewer_core::Panel::Field {
                    nx,
                    ny,
                    nz,
                    extent_m,
                    ..
                } => {
                    detail.push(("grid".into(), format!("{nx} x {ny} x {nz}")));
                    match extent_m {
                        Some(e) => {
                            detail.push(("extent".into(), fmt_box(*e)));
                            detail.push(("size".into(), fmt_span(*e)));
                        }
                        None => detail.push((
                            "extent".into(),
                            "not in this run file — framed in cell units".into(),
                        )),
                    }
                    NodeKind::Field
                }
                viewer_core::Panel::Points { positions, .. } => {
                    detail.push(("bodies".into(), (positions.len() / 3).to_string()));
                    detail.push(("bounds".into(), fmt_box(bounds)));
                    NodeKind::Points
                }
                viewer_core::Panel::Paths {
                    starts, vertices, ..
                } => {
                    detail.push(("paths".into(), starts.len().to_string()));
                    detail.push(("vertices".into(), (vertices.len() / 3).to_string()));
                    detail.push(("bounds".into(), fmt_box(bounds)));
                    NodeKind::Paths
                }
            };
            // The run's range, and this frame's — both, labelled, because they answer different
            // questions and a picture drawn on one is often read as the other.
            match scale {
                Some((lo, hi)) => detail.push((
                    "range (run)".into(),
                    format!("{} .. {}", magnitude(lo), magnitude(hi)),
                )),
                None => detail.push(("range (run)".into(), "one value throughout".into())),
            }
            let finite: Vec<f64> = panel
                .values()
                .iter()
                .copied()
                .filter(|v| v.is_finite())
                .collect();
            if finite.is_empty() {
                detail.push(("this frame".into(), "no finite samples".into()));
            } else {
                let lo = finite.iter().copied().fold(f64::MAX, f64::min);
                let hi = finite.iter().copied().fold(f64::MIN, f64::max);
                let mean = finite.iter().sum::<f64>() / finite.len() as f64;
                detail.push((
                    "this frame".into(),
                    format!(
                        "{} .. {}, mean {}",
                        magnitude(lo),
                        magnitude(hi),
                        magnitude(mean)
                    ),
                ));
            }
            let absent = panel.values().len() - finite.len();
            if absent > 0 {
                // A hole is not a value, and a count of them is the difference between "the
                // picture has gaps" and "the picture is wrong".
                detail.push(("absent".into(), format!("{absent} sample(s) with no value")));
            }
            detail.push((
                "colour".into(),
                if scale_is_signed(scale) {
                    "diverging, neutral at zero".into()
                } else {
                    "sequential".into()
                },
            ));
            push(Node {
                path: format!("/run/{}", panel.name()),
                name: panel.name().to_string(),
                kind,
                parent: Some(group),
                depth: 2,
                bounds: Some(bounds),
                unit: panel.unit().to_string(),
                detail,
            });
        }
    }

    // Readings, grouped by the domain that reports them — including domains with no geometry at
    // all, which is the half of a scene an outliner built from boxes alone would lose.
    if !f.readings.is_empty() {
        let group = push(Node {
            path: "/readings".into(),
            name: "Readings".into(),
            kind: NodeKind::Group,
            parent: Some(root),
            depth: 1,
            bounds: None,
            unit: String::new(),
            detail: vec![("count".into(), f.readings.len().to_string())],
        });
        let mut domains: Vec<&str> = Vec::new();
        for r in &f.readings {
            if !domains.contains(&r.domain.as_str()) {
                domains.push(&r.domain);
            }
        }
        for domain in domains {
            let mine: Vec<_> = f.readings.iter().filter(|r| r.domain == domain).collect();
            let d = push(Node {
                path: format!("/readings/{domain}"),
                name: domain.to_string(),
                kind: NodeKind::Group,
                parent: Some(group),
                depth: 2,
                bounds: checked
                    .boxes
                    .iter()
                    .find(|b| b.name == domain)
                    .map(|b| box_of(&b.corners)),
                unit: String::new(),
                detail: vec![("scalars".into(), mine.len().to_string())],
            });
            for r in mine {
                // The whole history, not just now: the range over the run and where this frame
                // sits in it, which is what a scalar's row is actually asked.
                let series: Vec<f64> = run
                    .frames
                    .iter()
                    .filter_map(|fr| {
                        fr.readings
                            .iter()
                            .find(|q| q.domain == r.domain && q.label == r.label)
                            .map(|q| q.value)
                    })
                    .filter(|v| v.is_finite())
                    .collect();
                let mut detail = vec![
                    ("value".into(), format!("{} {}", magnitude(r.value), r.unit)),
                    ("unit".into(), r.unit.clone()),
                ];
                if !series.is_empty() {
                    let lo = series.iter().copied().fold(f64::MAX, f64::min);
                    let hi = series.iter().copied().fold(f64::MIN, f64::max);
                    detail.push((
                        "range (run)".into(),
                        format!("{} .. {}", magnitude(lo), magnitude(hi)),
                    ));
                    detail.push((
                        "first, last".into(),
                        format!(
                            "{}, {}",
                            magnitude(series[0]),
                            magnitude(series[series.len() - 1])
                        ),
                    ));
                }
                if series.len() < run.frames.len() {
                    detail.push((
                        "absent".into(),
                        format!("{} frame(s)", run.frames.len() - series.len()),
                    ));
                }
                push(Node {
                    path: format!("/readings/{}/{}", r.domain, r.label),
                    name: r.label.clone(),
                    kind: NodeKind::Reading,
                    parent: Some(d),
                    depth: 3,
                    bounds: None,
                    unit: r.unit.clone(),
                    detail,
                });
            }
        }
    }

    Tree { nodes }
}

fn box_of(corners: &[[f64; 3]; 8]) -> [f64; 6] {
    let mut b = [f64::MAX, f64::MAX, f64::MAX, f64::MIN, f64::MIN, f64::MIN];
    for c in corners {
        for a in 0..3 {
            b[a] = b[a].min(c[a]);
            b[a + 3] = b[a + 3].max(c[a]);
        }
    }
    b
}

fn fmt_point(p: [f64; 3]) -> String {
    format!("({:.4}, {:.4}, {:.4})", p[0], p[1], p[2])
}

/// The eight corners of an axis-aligned box, in the order [`PlacedBox::corners`] uses.
///
/// Bit 0 is x, bit 1 is y, bit 2 is z, and a clear bit takes the low face — the order [`EDGES`]
/// is written against, so a box built here wireframes correctly with no further agreement.
///
/// Here because a **run** can now say where its field was sampled: `extent_m` travels in the wire
/// format, so a shell drawing a run needs a box from six numbers without a scene beside it. That
/// is the case the standalone viewer has always been in and the editor is in whenever a run is
/// opened without the file that produced it.
pub fn corners_of(extent: [f64; 6]) -> [[f64; 3]; 8] {
    let pick = |i: usize, a: usize| {
        if i & (1 << a) == 0 {
            extent[a]
        } else {
            extent[a + 3]
        }
    };
    std::array::from_fn(|i| [pick(i, 0), pick(i, 1), pick(i, 2)])
}
