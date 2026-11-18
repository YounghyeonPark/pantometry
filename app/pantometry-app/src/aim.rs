//! **Pointing at the model in the view.**
//!
//! What a click in the viewport sets depends on what is selected, the way a simulation tool's
//! graphics window assigns whatever boundary condition the model tree has open:
//!
//! - **an entry that says where** — a cooled face, a held or pressed one, a heat source, a region
//!   — is pointed at its domain's faces. A click on a face sets an entry's `face`; a shift-drag
//!   across one draws a patch on it for a cooled face, and the box's footprint for a heat source or
//!   a region, keeping its depth;
//! - **the results panel**, on the probe, the line, the slice or the region, is pointed at the
//!   domain it is reading: a click puts the probe there, takes the plane through it, or sets a
//!   line's start and then its end; a shift-drag draws the line, or the region's box.
//!
//! The arithmetic is [`editor_core::pick`]'s, so the browser editor can make the same calls; this is
//! the part that knows about a pointer and a painter. The ray is [`viewer_core::Camera::ray`], the
//! inverse of the projection everything is drawn with, so what is drawn is what is hit.
//!
//! **Only the aimed-at domain can be hit.** Two domains often share a box — scene 25's structure is
//! its block's own grid — and a click meant for the selected entry's domain should not be taken by
//! the box in front of it, so the others are looked through.

use crate::results::{Aimed, Tab};
use editor_core::pick::{self, Hit, Patch, Stated, Target};
use eframe::egui;
use pantometry::scene::Placed;
use pantometry_world::templates::Placing;
use pantometry_world::FaceSpec;

/// What a click or a shift-drag in the view would set.
pub enum Aim {
    /// An entry of a list that says where on its domain it is.
    Entry {
        /// The domain, with its grid.
        target: Target,
        /// The list.
        key: String,
        /// Which entry of it.
        index: usize,
        /// How the format says it is placed.
        placing: Placing,
        /// What it states now, for drawing.
        stated: Option<Stated>,
    },
    /// The results panel's probe, line, box or plane.
    Results(Aimed),
}

/// A shift-drag that began on a face and is drawing on it rather than turning the view.
#[derive(Clone, Copy, Debug)]
pub struct Drag {
    /// The face it began on.
    pub face: FaceSpec,
    /// Where it began, in the domain's frame.
    pub start: [f64; 3],
    /// Where it is now.
    pub end: [f64; 3],
}

/// `viewer_core`'s placement as the library's: the same two fields, in a crate the viewer does not
/// link.
fn placed(p: viewer_core::Placed) -> Placed {
    Placed {
        at_m: p.at_m,
        turn: p.turn,
    }
}

impl Aim {
    /// Its box, in its own frame.
    pub fn extent(&self) -> [f64; 6] {
        match self {
            Aim::Entry { target, .. } => target.extent,
            Aim::Results(a) => a.extent,
        }
    }

    /// Where that frame sits.
    pub fn place(&self) -> Placed {
        match self {
            Aim::Entry { target, .. } => target.place,
            Aim::Results(a) => placed(a.place),
        }
    }

    /// Where `ray` meets the domain aimed at.
    pub fn hit(&self, ray: &viewer_core::Ray) -> Option<Hit> {
        pick::hit(self.extent(), self.place(), ray)
    }

    /// Whether a shift-drag draws something here.
    pub fn drags(&self) -> bool {
        match self {
            Aim::Entry { .. } => true,
            Aim::Results(a) => matches!(a.tab, Tab::Line | Tab::Region),
        }
    }

    /// Where a drag begun on `face` has reached, for the pointer's `ray`.
    ///
    /// A box stays on the face it was begun on, so its end is the ray against that face's plane,
    /// held to the face. A line may cross the domain to any face, so its end is wherever the ray
    /// meets it, and the face's plane when it meets nothing.
    pub fn end(&self, face: FaceSpec, ray: &viewer_core::Ray) -> Option<[f64; 3]> {
        if let Aim::Results(a) = self {
            if a.tab == Tab::Line {
                if let Some(h) = self.hit(ray) {
                    return Some(h.local);
                }
            }
        }
        pick::on_plane(self.extent(), self.place(), face, ray)
    }

    /// The cells a drag on a face spans, for an entry: from the cell under its start to the cell
    /// under its end, inclusive.
    pub fn patch(&self, drag: &Drag) -> Option<Patch> {
        let Aim::Entry { target, .. } = self else {
            return None;
        };
        let a = target.on(drag.face, target.cell_at(drag.start));
        let b = target.on(drag.face, target.cell_at(drag.end));
        Some(pick::patch_between(a, b))
    }

    /// What a click and a shift-drag would do here, for the readout's second line: the place a
    /// person looks when they want to know what the pointer is for.
    pub fn hint(&self) -> &'static str {
        match self {
            Aim::Entry {
                placing: Placing::FacePatch,
                ..
            } => "click: this face · shift-drag: a patch on it",
            Aim::Entry {
                placing: Placing::Face,
                ..
            } => "click: this face",
            Aim::Entry {
                placing: Placing::Cells,
                ..
            } => "shift-drag: the box's footprint here",
            Aim::Results(a) => match a.tab {
                Tab::Point => "click: the probe here",
                Tab::Slice => "click: the plane through here",
                Tab::Line => "click: a line end · shift-drag: the line",
                Tab::Region => "shift-drag: the region's box",
                Tab::Readings => "",
            },
        }
    }

    /// What is under the pointer, as the readout beside it says.
    pub fn readout(&self, hit: &Hit) -> String {
        match self {
            Aim::Entry {
                target, placing, ..
            } => {
                let c = target.cell_at(hit.local);
                let [u, v] = target.on(hit.face, c);
                // A void cell is the box's and not the part's: a face pointed at there is a face
                // the part may not reach, which the area a click writes then says.
                let void = if target.is_solid(c) { "" } else { " · void" };
                match placing {
                    Placing::Face | Placing::FacePatch => {
                        format!("{} · cell {u}, {v}{void}", hit.face.spelled())
                    }
                    Placing::Cells => {
                        format!(
                            "{} · cell {}, {}, {}{void}",
                            hit.face.spelled(),
                            c[0],
                            c[1],
                            c[2]
                        )
                    }
                }
            }
            Aim::Results(a) => {
                let (unit, per) = a.unit;
                // The panel's own spelling, so what the pointer is over reads the way the panel
                // and the status line will say it once it is clicked.
                match a.tab {
                    Tab::Slice => format!(
                        "{} = {} {unit}",
                        post_axis(a.axis),
                        crate::results::short(hit.local[a.axis] / per)
                    ),
                    _ => crate::results::spelled(hit.local, unit, per),
                }
            }
        }
    }

    /// The faces of the aimed-at domain that face the eye, each with its centre: where its name
    /// is written while an entry is being pointed.
    pub fn facing(&self, eye: [f64; 3]) -> Vec<(FaceSpec, [f64; 3])> {
        facing(self.extent(), self.place(), eye)
    }

    /// The four corners, in world metres, of what a drag would set: a patch of cells for an
    /// entry, the rectangle it spans on the face for the region. `None` for a line.
    pub fn outline(&self, drag: &Drag) -> Option<[[f64; 3]; 4]> {
        match self {
            Aim::Entry { target, .. } => Some(target.patch_corners(drag.face, self.patch(drag)?)),
            Aim::Results(a) if a.tab == Tab::Region => {
                let [p, q] = drag.face.axes();
                let place = self.place();
                let corner = |u: f64, v: f64| {
                    let mut c = drag.start;
                    c[p] = u;
                    c[q] = v;
                    place.apply(c)
                };
                let (s, t) = (drag.start, drag.end);
                Some([
                    corner(s[p], s[q]),
                    corner(t[p], s[q]),
                    corner(t[p], t[q]),
                    corner(s[p], t[q]),
                ])
            }
            Aim::Results(_) => None,
        }
    }
}

/// The faces of a box that face the eye, each with its centre in world metres: those whose outward
/// normal, turned into the world, points to the eye's side of the face.
pub fn facing(e: [f64; 6], place: Placed, eye: [f64; 3]) -> Vec<(FaceSpec, [f64; 3])> {
    let turn = Placed {
        at_m: [0.0; 3],
        turn: place.turn,
    };
    FaceSpec::ALL
        .into_iter()
        .filter_map(|face| {
            let (axis, high) = face.normal();
            let mut centre = [0.0; 3];
            let mut normal = [0.0; 3];
            for k in 0..3 {
                centre[k] = 0.5 * (e[k] + e[k + 3]);
            }
            centre[axis] = if high { e[axis + 3] } else { e[axis] };
            normal[axis] = if high { 1.0 } else { -1.0 };
            let (c, n) = (place.apply(centre), turn.apply(normal));
            let toward = (0..3).map(|k| n[k] * (eye[k] - c[k])).sum::<f64>();
            (toward > 0.0).then_some((face, c))
        })
        .collect()
}

/// An axis's name, the results panel's spelling.
fn post_axis(axis: usize) -> &'static str {
    editor_core::post::AXIS_NAMES[axis.min(2)]
}

/// What the selected entry states, drawn on the model: its face filled, with its patch inside it
/// when there is one, or the box of cells as twelve edges. A face turned away from the eye is
/// drawn as an outline only, since a fill on it would read as the face in front.
pub fn draw_stated(
    painter: &egui::Painter,
    target: &Target,
    stated: &Stated,
    eye: [f64; 3],
    to_screen: &dyn Fn([f64; 3]) -> Option<egui::Pos2>,
) {
    let ink = egui::Color32::from_rgb(255, 190, 60);
    match stated {
        Stated::Face { face, patch } => {
            let whole = target.whole(*face);
            let front = facing(target.extent, target.place, eye)
                .iter()
                .any(|(f, _)| f == face);
            let fill = if front {
                egui::Color32::from_rgba_unmultiplied(255, 190, 60, 40)
            } else {
                egui::Color32::TRANSPARENT
            };
            quad(
                painter,
                &target.patch_corners(*face, whole),
                fill,
                egui::Stroke::new(1.5_f32, ink),
                to_screen,
            );
            if let Some(p) = patch {
                let strong = if front {
                    egui::Color32::from_rgba_unmultiplied(255, 190, 60, 110)
                } else {
                    egui::Color32::TRANSPARENT
                };
                quad(
                    painter,
                    &target.patch_corners(*face, *p),
                    strong,
                    egui::Stroke::new(2.0_f32, ink),
                    to_screen,
                );
            }
        }
        Stated::Cells { from, to } => {
            let corners = target.cells_corners(*from, *to);
            for (a, b) in editor_core::EDGES {
                if let (Some(p), Some(q)) = (to_screen(corners[a]), to_screen(corners[b])) {
                    painter.line_segment([p, q], egui::Stroke::new(2.0_f32, ink));
                }
            }
        }
    }
}

/// A quadrilateral in the world, filled and outlined, when all four corners are on screen.
pub fn quad(
    painter: &egui::Painter,
    corners: &[[f64; 3]; 4],
    fill: egui::Color32,
    stroke: egui::Stroke,
    to_screen: &dyn Fn([f64; 3]) -> Option<egui::Pos2>,
) {
    let points: Option<Vec<egui::Pos2>> = corners.iter().map(|c| to_screen(*c)).collect();
    let Some(points) = points else {
        return;
    };
    painter.add(egui::Shape::convex_polygon(points, fill, stroke));
}

/// The readout beside the pointer: the probe's box, so the two read as one kind of thing, held
/// inside the view.
pub fn draw_readout(
    painter: &egui::Painter,
    rect: egui::Rect,
    at: egui::Pos2,
    text: &str,
    visuals: &egui::Visuals,
) {
    let font = egui::FontId::monospace(11.0);
    let galley = painter.layout_no_wrap(text.to_string(), font, visuals.strong_text_color());
    let pad = egui::vec2(6.0, 3.0);
    let size = galley.size() + pad * 2.0;
    let mut origin = at + egui::vec2(12.0, -size.y - 12.0);
    origin.x = origin.x.min(rect.right() - size.x).max(rect.left());
    origin.y = origin.y.max(rect.top()).min(rect.bottom() - size.y);
    let box_rect = egui::Rect::from_min_size(origin, size);
    painter.rect_filled(box_rect, 3.0, visuals.extreme_bg_color);
    painter.rect_stroke(
        box_rect,
        3.0,
        egui::Stroke::new(1.0_f32, visuals.weak_text_color()),
        egui::StrokeKind::Inside,
    );
    painter.galley(origin + pad, galley, visuals.strong_text_color());
}
