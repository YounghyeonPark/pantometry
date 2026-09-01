//! **Every panel is in its own frame, and says where that frame is.**
//!
//! It used to be two frames. `capture` sampled a field in the domain's coordinates and dropped
//! the pose — `let _ = pose` — while bodies had it multiplied into their positions. One file,
//! local and world coordinates, and nothing in it saying which a panel used. `Panel::bounds`' own
//! doc said so and declined to promise anything about overlaying them, which is a hazard written
//! down rather than closed.
//!
//! It was measurable and it was measured: two blocks half a metre apart exported to glTF **on top
//! of each other**, both at their local origin. Nothing caught it because **no shipped scene
//! stated a pose** — under the identity, local and world are the same thing, and all twenty-nine
//! agreed. Scene 30 is the one that does not, and it was written for that reason.
//!
//! So every test here uses a real pose, which is the whole point. The domain is invented in this
//! file, the way `knows_no_physics.rs` invents one: what is being checked is the scene layer, and
//! a real domain would bring its own reasons for numbers to move.

use pantometry_core::{Bodies, Domain, Exchange, Reading, ScalarField, Schedule, Simulation};
use pantometry_core::{Pose, Violation};
use pantometry_scene::{capture, Extent, PanelData, Placed, Placement};
use pantometry_units::{LengthVec, Time};
use std::collections::BTreeMap;

/// Two motes a metre apart on the local x axis, and a field that reads back the x it is asked
/// about — so every number below is a coordinate and nothing is a physical result.
struct Two;

impl Domain for Two {
    fn name(&self) -> &str {
        "two"
    }
    fn step(&mut self, _dt: Time, _t: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        Ok(())
    }
    fn as_field(&self) -> Option<&dyn ScalarField> {
        Some(self)
    }
    fn as_bodies(&self) -> Option<&dyn Bodies> {
        Some(self)
    }
    fn readings(&self) -> Vec<Reading> {
        Vec::new()
    }
}

impl ScalarField for Two {
    fn at(&self, p: LengthVec, _t: Time) -> f64 {
        p.to_si().x
    }
    fn unit(&self) -> &'static str {
        "m"
    }
}

impl Bodies for Two {
    fn count(&self) -> usize {
        2
    }
    fn position(&self, i: usize) -> LengthVec {
        LengthVec::m(i as f64, 0.0, 0.0)
    }
    fn value(&self, _i: usize) -> f64 {
        1.0
    }
    fn value_unit(&self) -> &'static str {
        "1"
    }
    fn cell(&self) -> Option<(LengthVec, LengthVec)> {
        Some((LengthVec::m(0.0, 0.0, 0.0), LengthVec::m(1.0, 1.0, 1.0)))
    }
}

/// One frame, with the domain placed as given.
fn framed(placement: Placement) -> pantometry_scene::Frame {
    let mut sim = Simulation::new(Schedule::Staggered).with(Two);
    sim.advance(Time::s(1.0)).expect("it claims nothing");
    let mut placed = BTreeMap::new();
    placed.insert("two".to_string(), placement);
    capture(&sim, &placed)
}

/// The extent the field is sampled over: the unit interval along x, two samples.
fn along_x() -> Extent {
    Extent::new(
        LengthVec::m(0.0, 0.0, 0.0),
        LengthVec::m(1.0, 0.0, 0.0),
        2,
        1,
        1,
    )
}

#[test]
fn the_identity_is_exact() {
    // Twenty-nine of the thirty shipped scenes are this case, so it has to be exact rather than
    // near, and not approximately so: a writer decides
    // by `is_here` whether to emit a transform at all, and "close enough" would make that
    // decision depend on rounding.
    assert!(Placed::HERE.is_here());
    assert_eq!(Placed::HERE.at_m, [0.0; 3]);
    assert_eq!(
        Placed::HERE.turn,
        [0.0, 0.0, 0.0, 1.0],
        "a unit quaternion in x, y, z, w"
    );
    assert_eq!(Placed::default(), Placed::HERE);

    let frame = framed(Placement::default());
    for panel in &frame.panels {
        assert!(
            panel.place.is_here(),
            "{} moved without being placed",
            panel.name
        );
    }
}

#[test]
fn a_body_stays_in_its_own_frame_and_the_panel_carries_the_pose() {
    // These used to come out at 0.5 and 1.5 with the panel saying nothing. They come out at 0
    // and 1 with the panel saying "half a metre along x" — the same placement stated instead of
    // baked, which is what lets a reader put a field beside them.
    let frame = framed(Placement::at(Pose::at(LengthVec::m(0.5, 0.0, 0.0))));
    let panel = frame
        .panels
        .iter()
        .find(|p| matches!(p.data, PanelData::Points { .. }))
        .expect("bodies");

    assert_eq!(panel.place.at_m, [0.5, 0.0, 0.0]);
    assert!(!panel.place.is_here());

    let PanelData::Points {
        positions, bounds, ..
    } = &panel.data
    else {
        unreachable!()
    };
    assert_eq!(positions[0], [0.0, 0.0, 0.0], "the pose was baked in");
    assert_eq!(positions[1], [1.0, 0.0, 0.0]);
    assert_eq!(bounds[0], 0.0, "the cell was baked in");
    assert_eq!(bounds[3], 1.0);
}

#[test]
fn a_rotated_cell_keeps_its_shape_instead_of_becoming_a_bounding_box() {
    // **The information the old shape lost.** Bodies used to pose the cell's two corners and take
    // the min and max, which for any rotation is a box *around* the cell rather than the cell. A
    // quarter turn of a unit cube gave a bounding box √2 across on two axes, and nothing written
    // down could recover the cube.
    let quarter = std::f64::consts::FRAC_PI_2;
    let frame = framed(Placement::at(Pose::turned(glam::DQuat::from_rotation_z(
        quarter,
    ))));
    let panel = frame
        .panels
        .iter()
        .find(|p| matches!(p.data, PanelData::Points { .. }))
        .expect("bodies");

    let PanelData::Points { bounds, .. } = &panel.data else {
        unreachable!()
    };
    for a in 0..2 {
        let across = bounds[a + 3] - bounds[a];
        assert!(
            (across - 1.0).abs() < 1e-12,
            "axis {a} is {across} across, not 1 — the cell became a bounding box"
        );
    }

    // And the turn is on the panel, exactly: a quarter turn about z is sin and cos of an eighth.
    let (s, c) = ((quarter / 2.0).sin(), (quarter / 2.0).cos());
    for (a, want) in [0.0, 0.0, s, c].into_iter().enumerate() {
        assert!(
            (panel.place.turn[a] - want).abs() < 1e-12,
            "quaternion component {a} is {}, not {want}",
            panel.place.turn[a]
        );
    }
}

#[test]
fn a_field_is_still_sampled_in_its_own_coordinates() {
    // Unchanged, and deliberately: how long the bar is should not change when the bar is moved.
    // What changed is that the panel says where those coordinates are instead of dropping the
    // pose on the floor.
    let mut placement = Placement::at(Pose::at(LengthVec::m(9.0, 0.0, 0.0)));
    placement.extent = Some(along_x());
    let frame = framed(placement);
    let panel = frame
        .panels
        .iter()
        .find(|p| matches!(p.data, PanelData::Field { .. }))
        .expect("a field");

    let PanelData::Field {
        extent_m, values, ..
    } = &panel.data
    else {
        unreachable!()
    };
    assert_eq!(
        extent_m[0], 0.0,
        "the extent moved with the pose; it is the domain's own box"
    );
    assert_eq!(extent_m[3], 1.0);
    assert_eq!(
        values[0], 0.0,
        "the field was sampled at the posed position"
    );
    assert_eq!(panel.place.at_m, [9.0, 0.0, 0.0]);
}

#[test]
fn two_domains_in_two_places_come_back_in_two_places() {
    // **The measured defect, in one frame.** Two blocks half a metre apart exported to glTF on
    // top of each other, both at their local origin, because the pose never reached the file.
    // Each panel must now report its own placement -- and the *values* must be identical,
    // because the two domains are identical and only their placement differs.
    //
    // A domain produces one panel, not two: `capture` takes a field if the domain offers one and
    // falls through to bodies otherwise. So this needs two domains, which is what the defect
    // needed too.
    let mut sim = Simulation::new(Schedule::Staggered)
        .with(Named("near"))
        .with(Named("far"));
    sim.advance(Time::s(1.0)).expect("it claims nothing");

    let mut placed = BTreeMap::new();
    placed.insert("near".to_string(), Placement::default());
    placed.insert(
        "far".to_string(),
        Placement::at(Pose::at(LengthVec::m(0.5, 0.0, 0.0))),
    );
    let frame = capture(&sim, &placed);

    let near = frame
        .panels
        .iter()
        .find(|p| p.name == "near")
        .expect("near");
    let far = frame.panels.iter().find(|p| p.name == "far").expect("far");

    assert!(near.place.is_here(), "near was moved: {:?}", near.place);
    assert_eq!(
        far.place.at_m,
        [0.5, 0.0, 0.0],
        "far does not say where it is"
    );

    // The same object, so the same numbers. Before, `far`'s were shifted and its placement was
    // nowhere; the shift has moved from the values to the placement and nothing else changed.
    assert_eq!(
        near.bounds(),
        far.bounds(),
        "two identical domains disagree about their own extent"
    );
}

/// The same domain under a chosen name, so one frame can hold two of them.
struct Named(&'static str);

impl Domain for Named {
    fn name(&self) -> &str {
        self.0
    }
    fn step(&mut self, _dt: Time, _t: Time, _bus: &mut Exchange) -> Result<(), Violation> {
        Ok(())
    }
    fn as_field(&self) -> Option<&dyn ScalarField> {
        None
    }
    fn as_bodies(&self) -> Option<&dyn Bodies> {
        Some(&Two)
    }
    fn readings(&self) -> Vec<Reading> {
        Vec::new()
    }
}

#[test]
fn a_turn_comes_back_as_the_axis_and_angle_it_was_written_as() {
    // **A round trip through the format's two vocabularies.** A scene states a rotation as an axis
    // and an angle, because a quaternion's four numbers carry a constraint between them that no
    // file can express and no person can check by eye. A run carries the quaternion, because that
    // is what a glTF node and a USD prim take. Anything showing a placement to a person has to
    // come back, and this is the closed form for that: what goes in comes out.
    let cases: [([f64; 3], f64); 5] = [
        ([0.0, 0.0, 1.0], 90.0),
        ([0.0, 0.0, 1.0], 60.0),
        ([1.0, 0.0, 0.0], 45.0),
        ([0.0, 1.0, 0.0], 30.0),
        ([1.0, 1.0, 1.0], 120.0),
    ];
    for (axis, degrees) in cases {
        let n = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        let unit = [axis[0] / n, axis[1] / n, axis[2] / n];
        let q = glam::DQuat::from_axis_angle(
            glam::DVec3::new(unit[0], unit[1], unit[2]),
            degrees.to_radians(),
        );
        let (got_axis, got_deg) = Placed::of(Pose::turned(q)).axis_angle();
        assert!(
            (got_deg - degrees).abs() < 1e-9,
            "{degrees}° about {axis:?} came back as {got_deg}°"
        );
        for a in 0..3 {
            assert!(
                (got_axis[a] - unit[a]).abs() < 1e-9,
                "axis came back as {got_axis:?}, not {unit:?}"
            );
        }
    }
}

#[test]
fn a_half_turn_is_where_acos_would_have_lost_its_precision() {
    // **Why `atan2` and not `acos(w)`, measured rather than remembered.** The comment this
    // replaced said `acos` loses precision near a *half turn*. It is the other end: at a half turn
    // `w` is zero and `|d acos/dw|` is 1, the best-conditioned point there is. At a thousandth of
    // a degree `w` -> 1 and the derivative is 1.1e5.
    //
    //     turn        acos        atan2
    //     0.0001°     4.6e-5      0
    //     0.001°      4.0e-7      2.2e-16
    //     1°          2.8e-13     0
    //     180°        0           0
    //
    // A thousandth of a degree is a real thing to ask of a rotation ring, and 180° is not a case
    // that needed guarding at all.
    let tiny = 1e-3_f64;
    let q = glam::DQuat::from_axis_angle(glam::DVec3::Z, tiny.to_radians());
    let (axis, degrees) = Placed::of(Pose::turned(q)).axis_angle();
    assert!(
        (degrees - tiny).abs() / tiny < 1e-9,
        "{tiny}° came back as {degrees}°, off by {:.2e} relative",
        (degrees - tiny).abs() / tiny
    );
    assert!((axis[2] - 1.0).abs() < 1e-9);

    // And the half turn itself, where `w` is zero and the axis is all that is left.
    let half = glam::DQuat::from_axis_angle(glam::DVec3::X, std::f64::consts::PI);
    let (axis, degrees) = Placed::of(Pose::turned(half)).axis_angle();
    assert!((degrees - 180.0).abs() < 1e-9, "{degrees}");
    assert!((axis[0].abs() - 1.0).abs() < 1e-9, "{axis:?}");
}

#[test]
fn the_identity_has_no_axis_and_is_given_one_anyway() {
    // Normalising a zero vector is a `NaN` per component, and a caption reading
    // `turned NaN° about (NaN, NaN, NaN)` is worse than one reading nothing. So the identity
    // answers `+z` and zero degrees -- a direction that means nothing because the angle is
    // nothing, which is the only honest thing to return.
    let (axis, degrees) = Placed::HERE.axis_angle();
    assert_eq!(degrees, 0.0);
    assert_eq!(axis, [0.0, 0.0, 1.0]);
    assert!(axis.iter().all(|v| v.is_finite()));
}
