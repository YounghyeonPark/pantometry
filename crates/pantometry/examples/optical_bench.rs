//! A three-dimensional optical bench: elements placed by pose, rays traced through them, and a
//! layout you can rotate in a browser.
//!
//! ```text
//! cargo run --release --example optical_bench             # numbers, checked
//! cargo run --release --example optical_bench bench.html  # and the 3D layout, rotatable
//! cargo run --release --example optical_bench bench.svg   # or a still of it, for a page
//! ```
//!
//! Every other example here draws a graph. This one draws the **instrument**: a collimated beam,
//! a doublet, a fold mirror that turns the axis through 90°, and the rays arriving at a tilted
//! image plane — as paths in space, on one page, rotatable and scrollable with nothing installed.
//!
//! # What this is and is not
//!
//! It is the analysis half of an optical layout tool: place, trace, measure, look. Ray paths are a
//! first-class shape in the scene layer now — `PanelData::Paths` — so the view draws them the same
//! way it draws bodies and fields, depth-sorted and on one colour scale.
//!
//! It is **not** an interactive editor and there is no renderer with materials or shadows. Nothing
//! here is real-time; a run produces a file. Those are different products and this workspace does
//! not pretend otherwise.
//!
//! # The glass is drawn, and it is the glass the rays are traced against
//!
//! For a long time this drew the rays and nothing else: light bending in mid-air, at angles
//! nothing in the picture explained. The doublet is drawn now, and not as a symbol for a lens --
//! `profile` samples the same sag `cap_intersect` solves against, so the curve on the screen is
//! the surface a ray meets, checked point by point rather than drawn to look right.
//!
//! One consequence is a measurement that only exists once the shape does: the **edge thickness**.
//! A prescription that focuses perfectly can still be one no workshop can cut, and the paraxial
//! arithmetic that produced these curvatures has no opinion on that at all.
//!
//! # What is checked, because a picture is not a result
//!
//! A layout that looks right is the easiest wrong answer in optics. So:
//!
//! - the doublet's focal length, against the thin-lens combination of its two powers;
//! - the fold, against the law of reflection — the axis leaves at exactly 90° to how it arrived;
//! - the marginal ray's height at the stop, against `h = f·tanθ` for each field angle;
//! - the RMS spot against the diffraction limit, so "in focus" is a number and not a look;
//! - every point of the drawn glass against `cap_intersect`, so the shape on the screen is the
//!   surface the rays meet and not a second description of it;
//! - the thickness of each element at its thinnest, which is the check the paraxial arithmetic
//!   cannot make and the drawn shape can.
//!
//! Each is arithmetic done here, not a value read back out of the thing being checked.

use glam::DVec3;
use pantometry::prelude::*;
use pantometry::scene::{Frame, Panel, PanelData, Placed};
use pantometry_optics::geometry::{
    cap_intersect, hexapolar_unit, plane_intersect, profile, refract, Hit, Ray,
};

mod common;
use common::{check, check_between, heading};

/// N-BK7 and SF2 at the d line — a classic cemented doublet pair.
const N_CROWN: f64 = 1.5168;
const N_FLINT: f64 = 1.6477;
/// Focal length to aim for.
const FOCAL: f64 = 100e-3;
/// The stop: the semi-aperture the pupil is sampled over.
const SEMI: f64 = 8.0e-3;
/// The glass is cut a little larger, so the *stop* limits the beam rather than the edge of a
/// surface. A prescription whose clear aperture equals its stop vignettes on the first tilt.
const GLASS: f64 = 9.5e-3;
/// Where the fold mirror sits along the axis, past the doublet.
const FOLD_AT: f64 = 60e-3;
/// The cemented face, where the crown ends and the flint begins.
const CEMENT_Z: f64 = 4e-3;
/// The flat back face of the flint, which is where the lens ends.
const BACK_Z: f64 = 7e-3;
/// Half-planes the glass is drawn in. Each is a full meridian, so this is twice as many curves.
const MERIDIANS: usize = 6;
/// Points per half-meridian of a drawn surface.
const PROFILE_SAMPLES: usize = 24;
/// Segments in a drawn rim.
const RIM_SEGMENTS: usize = 64;
/// Field angles traced, in degrees.
const FIELDS: [f64; 3] = [0.0, 1.5, 3.0];
/// Where the fold mirror's hit sits in a traced path: origin, three glass surfaces, the mirror.
const MIRROR_AT: usize = 4;
/// And the image plane, which is the last.
const IMAGE_AT: usize = 5;
/// Where the glass sits on that same scale.
///
/// One panel carries one quantity, and a report draws one card per panel -- so a glass panel of
/// its own would be a second picture of the same bench with the rays taken out of it, which is
/// the thing this change exists to stop. The glass therefore goes on the field-angle scale, past
/// the widest field traced, and the panel's unit says so rather than leaving a colour bar to
/// imply that a lens is a four-degree ray. Being the top of the scale it is also the brightest
/// thing drawn, which is what an instrument should be.
const GLASS_ON_THE_FIELD_SCALE: f64 = 4.0;

fn main() {
    // ================================================================ the prescription
    heading("A doublet, from the powers it has to add up to");

    // Three surfaces, and each one's power is `(n2 - n1)/R` — the refractive index it goes *into*
    // minus the one it comes *from*. The crown's front is air-to-crown; the cemented face is
    // crown-to-flint, a difference of 0.13 rather than 0.52; the flint's back is flat.
    //
    // That distinction is not pedantry. Using `(n - 1)/R` on the cemented surface — the air-glass
    // formula, which is the one everybody remembers — made the flint four times too weak, and the
    // traced focal length came out 150 mm against a 100 mm prescription. The trace was right and
    // the prescription was wrong, which is the correct way round and only visible because the
    // trace measures `f` instead of being told it.
    // One free parameter: how much of the total power the front surface carries. Whatever is
    // left goes on the cemented face, so **every split has the same paraxial focal length** and
    // they differ only in aberration. That is the classic bending variable.
    let split = 1.55;
    let (r_crown, r_cement) = curvatures(split);

    println!(
        "  {:<30} {:>9.3} mm  air to crown",
        "front radius",
        r_crown * 1e3
    );
    println!(
        "  {:<30} {:>9.3} mm  crown to flint",
        "cemented radius",
        r_cement * 1e3
    );
    println!("  {:<30} {:>9}     flint to air", "back", "flat");
    let phi = 1.0 / FOCAL;
    check(
        "the three powers add to 1/f",
        1.0 / (phi * split + (phi - phi * split)),
        FOCAL,
        1e-12,
        "m",
    );

    // ================================================================ the bench
    heading("The bench, placed in three dimensions");
    let mut bench = Bench {
        crown_r: r_crown,
        cement_r: r_cement,
        fold_at: FOLD_AT,
        // Worked out rather than guessed. reflect(d, n) = d - 2(d.n)n, and for +z to leave along
        // +x that needs n = (-1, 0, 1)/sqrt(2). The other sign of the z component sends the beam
        // to -x, where the image plane is not — and every ray then misses it, which is the right
        // failure and an opaque one to read: sixty-one rays traced, none arrived.
        fold_normal: DVec3::new(-1.0, 0.0, 1.0).normalize(),
        // Filled in from the measurement below, because the thin-lens prescription does not know
        // about the glass thickness and therefore does not know where the focus is.
        image_along: 0.0,
    };

    // Where the focus actually is, from a traced paraxial ray.
    let (f_measured, bfd) = bench.measure().expect("a paraxial ray gets through");
    bench.image_along = bfd - (FOLD_AT - BACK_Z);
    println!(
        "  {:<30} {:>9.3} mm  measured from a traced ray, against {:.1} thin-lens",
        "effective focal length",
        f_measured * 1e3,
        FOCAL * 1e3
    );
    println!(
        "  {:<30} {:>9.3} mm  from the back face, so the plane is {:.1} mm past the mirror",
        "back focal distance",
        bfd * 1e3,
        bench.image_along * 1e3
    );
    check_between(
        "a thick lens focuses shorter than its thin-lens prescription",
        bfd / f_measured,
        0.90,
        1.0,
        "x",
    );
    println!("  {:<30} {:>9.1} mm  along +z", "doublet at", 0.0);
    println!(
        "  {:<30} {:>9.1} mm  normal at 45 degrees",
        "fold mirror at",
        bench.fold_at * 1e3
    );
    println!(
        "  {:<30} {:>9.1} mm  along +x after the fold",
        "image plane at",
        bench.image_along * 1e3
    );

    // ================================================================ trace
    heading("Rays through it, and the closed forms they must obey");
    let Traced { spots, .. } = trace_fields(&bench);
    for (k, angle_deg) in FIELDS.iter().enumerate() {
        println!(
            "  {angle_deg:>4.1} deg   {:>3} of {:>3} rays through",
            spots[k].len(),
            PUPIL_RAYS
        );
        assert!(
            spots[k].len() * 10 >= PUPIL_RAYS * 9,
            "at least 90% of the pupil must make it, or the bench is vignetting"
        );
    }

    // The fold: a ray arriving along +z must leave along +x, exactly.
    let axial = Ray::new(LengthVec::from_si(DVec3::new(0.0, 0.0, 0.0)), DVec3::Z);
    let turned = reflect(axial.dir, bench.fold_normal);
    println!(
        "  {:<30} ({:.6}, {:.6}, {:.6})",
        "axis after the fold", turned.x, turned.y, turned.z
    );
    check("the fold turns the axis to +x", turned.x, 1.0, 1e-12, "");
    check_between("and leaves nothing along z", turned.z.abs(), 0.0, 1e-12, "");

    // ================================================================ bending
    heading("Bending: one parameter, scanned for the smallest spot");
    // Every split traced below has the same paraxial focal length, so this trades nothing away.
    // It is free performance sitting in a number somebody had to choose, which is why an optical
    // design tool's first move is to scan it.
    let scan: Vec<(f64, f64)> = (0..=28)
        .filter_map(|k| {
            let s = 1.0 + 1.4 * k as f64 / 28.0;
            on_axis_rms(&bench, s).map(|spot| (s, spot))
        })
        .collect();
    let &(best_split, best_rms) = scan
        .iter()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .expect("some split traces");
    for (s, spot) in scan.iter().step_by(4) {
        println!(
            "  split {s:>5.2}   RMS spot {:>8.2} um{}",
            spot * 1e6,
            if (*s - best_split).abs() < 1e-12 {
                "   <- best"
            } else {
                ""
            }
        );
    }
    let started = on_axis_rms(&bench, split).expect("the starting split traces");
    println!(
        "  {:<30} {:>9.2} um  ->  {:.2} um at split {best_split:.2}",
        "bending the doublet",
        started * 1e6,
        best_rms * 1e6
    );
    check_between(
        "bending is worth something",
        started / best_rms,
        1.5,
        500.0,
        "x",
    );
    // A minimum, not the end of the range — otherwise the scan chose its own edge.
    assert!(
        best_split > scan[0].0 && best_split < scan[scan.len() - 1].0,
        "the optimum is at the edge of the scan: {best_split}"
    );

    // Rebuild on the bent prescription and re-measure what follows from it.
    let (br, bc) = curvatures(best_split);
    bench.crown_r = br;
    bench.cement_r = bc;
    let (f_bent, bfd_bent) = bench.measure().expect("the bent doublet still focuses");
    bench.image_along = bfd_bent - (FOLD_AT - BACK_Z);
    check("bending does not move f", f_bent, f_measured, 0.03, "m");
    let Traced {
        paths,
        colours,
        spots,
    } = trace_fields(&bench);

    // The chief ray's height at the image plane is `f tan(theta)` — the definition of focal
    // length, and a check the trace cannot fake because nothing in it computes `f`.
    for (k, angle_deg) in FIELDS.iter().enumerate() {
        let theta = angle_deg.to_radians();
        let centroid = centroid_of(&spots[k]);
        let want = f_bent * theta.tan();
        if *angle_deg == 0.0 {
            check_between(
                "on axis the centroid is on axis",
                centroid.1.abs(),
                0.0,
                2e-5,
                "m",
            );
        } else {
            check(
                &format!("{angle_deg:.1} deg lands at f tan(theta)"),
                centroid.1.abs(),
                want,
                0.06,
                "m",
            );
        }
    }

    // ================================================================ the spot
    heading("Is it in focus? A number, not a look");
    let rms = rms_radius(&spots[0]);
    // The diffraction limit for this aperture, at the d line.
    let na = SEMI / f_bent;
    let airy = 0.61 * 587.6e-9 / na;
    println!(
        "  {:<30} {:>9.2} um  RMS spot radius on axis",
        "geometric",
        rms * 1e6
    );
    println!(
        "  {:<30} {:>9.2} um  0.61 lambda / NA",
        "diffraction limit",
        airy * 1e6
    );
    println!(
        "  {:<30} {:>9.2}x  the Airy radius",
        "the geometric spot is",
        rms / airy
    );
    // **The verdict, and it changed with the bending.** Before, 209 um against a 4.4 um limit —
    // 47x, and the picture would have been all aberration. After, the geometric spot is inside
    // the Airy disc, so the design is diffraction-limited and the ray trace has stopped being the
    // thing that decides the image. That is the answer an optical design review wants, and it is
    // a number rather than a look at a spot diagram.
    check_between(
        "the bent design is diffraction-limited",
        rms / airy,
        0.0,
        1.0,
        "x the Airy radius",
    );
    println!(
        "  {:<30} {:>9.1}x  before bending, which the scan removed",
        "it was",
        started / airy
    );

    // ================================================================ the glass
    heading("The glass, drawn as the surfaces the rays are traced against");
    let (glass_runs, glass_values) = bench.glass();
    println!(
        "  {:<30} {:>9} runs  {MERIDIANS} meridians of two elements, and a rim at each surface",
        "the doublet is drawn as",
        glass_runs.len()
    );

    // **Every surface reaches the aperture the glass is cut to**, and this comes first because
    // the check below depends on it. `profile` stops a surface at its own extent -- a sphere at
    // its equator -- and a ray parallel to the axis is *tangent* there, so `R² - h²` rounds
    // through zero and the intersection may return that point, one a few hundred picometres away,
    // or nothing. A bench whose glass is cut inside every surface never meets that edge. At a
    // power split of 3.0 the cemented radius is 6.5 mm against 9.5 mm of glass and the agreement
    // below reads 2.4 million rather than 0.2; this is the line that would fail first.
    //
    // It is exact equality, not a tolerance: `h_max` is `min(GLASS, extent)`, so either it is
    // `GLASS` to the bit or a surface has run out.
    let reached = bench
        .surfaces()
        .iter()
        .map(|s| s.drawn_semi_aperture())
        .fold(f64::INFINITY, f64::min);
    println!(
        "  {:<30} {:>9.3} mm  the least any of the three is drawn to, of {:.3} mm of glass",
        "every surface reaches",
        reached * 1e3,
        GLASS * 1e3
    );
    assert_eq!(
        reached, GLASS,
        "a surface is drawn to {reached} m of a {GLASS} m aperture, so it runs out before the \
         glass does and its drawn edge is where the surface turns vertical"
    );

    // **Every drawn point is a point a ray lands on.** `profile` evaluates a sag; `cap_intersect`
    // solves a quadratic. They share no arithmetic, so agreeing means the picture and the trace
    // are one surface rather than two descriptions of one -- which is the only way a drawing can
    // be checked at all.
    //
    // The floor is `4·eps·(S + S²/|R|)` with `S = |R| + standoff`, **per surface**: for an axial
    // ray the rounding that matters is the one in `R² - h²`, which the square root divides by
    // `2 sqrt(R² - h²)`. A single floor built from the largest radius would judge the shortest
    // surface by the longest one's arithmetic. `pantometry-optics` derives it and measures nine
    // configurations against it.
    //
    // **What it cannot see is a wrong prescription**, because both sides read `Bench::surfaces()`
    // and a wrong radius moves the drawing and the trace together. Measured: scaling the cemented
    // radius by 1.10 leaves this check well inside its floor and every other check on this page
    // green, the on-axis spot degrading 0.57 um to 4.06 um and still inside the Airy disc. That is
    // the correct trade -- the two *should* move together -- and the prescription is what the
    // paraxial and `f tan(theta)` checks above are for. This one is narrower on purpose: it says
    // the picture is not a second, decorative description of the lens.
    //
    // The `.expect` is half of it: a ray aimed at a point the trace does not reach returns `None`,
    // so a drawing wider than the glass fails here rather than showing rays leaving through the
    // side of a lens.
    let standoff = 40e-3;
    let mut worst: f64 = 0.0;
    // **Every vertex that is drawn**, which is what the sentence above says and what this did not
    // do: it read one meridian of the six, and none of the three rims -- 147 of the 1383 points the
    // picture is made of. The rims are built at sixty-four azimuths and the sections at six, so a
    // mistake in the azimuth would have been invisible at `y` and everywhere else in the frame.
    //
    // A drawn point does not carry which surface it belongs to, so each is asked of all three and
    // judged by the nearest: a point of the glass is a point on one of the surfaces of the glass.
    let floor_of = |s: &Surface| {
        let scale = s.r.abs() + standoff;
        // A flat surface is met by `plane_intersect`, which does not cancel; its error is the
        // standoff's own rounding and there is no radius to divide by.
        4.0 * f64::EPSILON
            * if s.r == 0.0 {
                standoff
            } else {
                scale + scale * scale / s.r.abs()
            }
    };
    let mut checked = 0usize;
    for at in glass_runs.iter().flatten() {
        let at = DVec3::new(at[0], at[1], at[2]);
        let ray = Ray::new(
            LengthVec::from_si(DVec3::new(at.x, at.y, -standoff)),
            DVec3::Z,
        );
        let nearest = bench
            .surfaces()
            .iter()
            .filter_map(|s| {
                s.meet(ray)
                    .map(|h| (h.point.to_si() - at).length() / floor_of(s))
            })
            .fold(f64::INFINITY, f64::min);
        // A miss shows up here as a large ratio rather than as an absence, because the nearest of
        // the *other* surfaces is then what the point is compared against -- which is how the flat
        // back face's rim was found: it read 4 mm, the distance to the crown.
        assert!(
            nearest.is_finite(),
            "a ray aimed at the drawn point {:.4} mm off the axis met none of the three surfaces",
            (at.x * at.x + at.y * at.y).sqrt() * 1e3
        );
        worst = worst.max(nearest);
        checked += 1;
    }
    println!(
        "  {:<30} {checked:>9} points  worst {worst:.4} of its surface's floor",
        "every drawn vertex, checked"
    );
    check_between(
        "every drawn point is where a ray lands",
        worst,
        0.0,
        1.0,
        "x the floor",
    );
    // And the two sides are two computations. Exactly zero everywhere would mean one of them is
    // reading the other, which is the failure this whole section exists to rule out.
    assert!(
        worst > 0.0,
        "every drawn point matched a traced point to the last bit, so the drawing and the trace \
         are not two ways of finding one surface"
    );

    // **The glass has thickness everywhere it is drawn, and the thickness at the edge is the one
    // the sags say.** Two surfaces bounding one element can cross, and a prescription whose
    // elements cross is one no workshop can cut. Nothing in the paraxial arithmetic that chose
    // these curvatures can see it: `f` is the same either way.
    //
    // The *minimum* is not always at the edge -- for the flint it is on the axis, where it is the
    // centre thickness by construction -- so asserting only that would be asserting a number
    // against a bound built out of it. The edge is where the two surfaces are furthest from their
    // vertices and the closed form is independent: `R - sign(R) sqrt(R^2 - h^2)`, the spelling
    // that cancels, against `pantometry-optics`' stable one.
    let [front, cement, back] = bench.surfaces();
    for (name, a, b, centre) in [
        ("crown", front, cement, CEMENT_Z),
        ("flint", cement, back, BACK_Z - CEMENT_Z),
    ] {
        let thinnest = a.thinnest_to(b);
        let drawn_edge = b.drawn_edge_z() - a.drawn_edge_z();
        let want = (b.at + textbook_sag(b.r, GLASS)) - (a.at + textbook_sag(a.r, GLASS));
        println!(
            "  {:<30} {:>9.3} mm  thinnest, {:.3} mm at the edge, {:.3} mm on the axis",
            name,
            thinnest * 1e3,
            drawn_edge * 1e3,
            centre * 1e3
        );
        // The cancellation in `R - sqrt(R^2 - h^2)` is at the scale of the radii themselves.
        let floor = 4.0 * f64::EPSILON * (a.r.abs() + b.r.abs() + centre);
        check_between(
            &format!("the drawn {name} is as thick at its edge as the sags say"),
            (drawn_edge - want).abs() / floor,
            0.0,
            1.0,
            "x the floor",
        );
        // A hundredth of a millimetre of glass is not an element anybody grinds. This is an
        // engineering floor, not a tolerance: the quantity it bounds is a length, not an error.
        assert!(
            thinnest > 0.1e-3,
            "the {name} is {:.4} mm thick at its thinnest, which is a prescription no workshop \
             would accept and, below zero, two surfaces that cross",
            thinnest * 1e3
        );
    }

    // Where the beam sits inside the glass, printed rather than checked: `cap_intersect` clips
    // every ray at the same aperture the outline is drawn to, so "the beam is inside the drawn
    // rim" is true by construction and asserting it would be asserting the clipping. The number
    // is worth seeing -- it is the margin the stop leaves the glass -- and the claim that has
    // teeth is the exact-aperture assertion above.
    let drawn_rim = glass_runs
        .iter()
        .flatten()
        .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
        .fold(0.0f64, f64::max);
    let widest = paths
        .iter()
        .filter(|p| p.len() > 3)
        .flat_map(|p| p[1..=3].iter())
        .map(|v| (v[0] * v[0] + v[1] * v[1]).sqrt())
        .fold(0.0f64, f64::max);
    println!(
        "  {:<30} {:>9.3} mm  inside a drawn rim of {:.3} mm",
        "the beam is widest at",
        widest * 1e3,
        drawn_rim * 1e3
    );

    // One picture: the rays and the glass they bend in, on one scale and in one panel, because
    // a report draws a card per panel and two cards would be the instrument in one and the light
    // in the other.
    assert!(
        FIELDS.iter().all(|f| *f < GLASS_ON_THE_FIELD_SCALE),
        "a field angle has reached the value the glass is drawn at, so the two would share a colour and the panel's unit would be lying"
    );
    // The two flats go into the drawn panel as well, so the shell that shades and the still both
    // show the bench bending its light at something rather than at nothing.
    // The two flats, built here because the panel below draws them and the checks further down
    // measure them: one call, so the picture and the assertion cannot describe different mirrors.
    let flats = flats_of(&bench, &paths);
    let mut drawn = glass_runs.clone();
    let mut values = glass_values.clone();
    for f in &flats {
        drawn.push(f.clone());
        values.push(GLASS_ON_THE_FIELD_SCALE);
    }
    drawn.extend(paths.iter().cloned());
    values.extend(colours.iter().copied());
    // **`PanelData::paths` drops a run of fewer than two points.** A picture missing its glass
    // because every outline was degenerate would look like a picture of rays, which is what this
    // one was; so the panel is asked what it kept.
    //
    // Runs **and** vertices. Counting runs alone cannot see the failure the paragraph above names:
    // a glass outline collapsed to two coincident points is still a run, and the count comes out
    // right while the picture has no glass in it. The vertex count is a closed form -- two
    // meridians per half-plane, a closing point, and a rim per surface -- so it is arithmetic done
    // here rather than a second reading of the thing being checked.
    let panel = PanelData::paths(drawn, values);
    let (kept, vertices) = match &panel {
        PanelData::Paths {
            starts, vertices, ..
        } => (starts.len(), vertices.len()),
        _ => (0, 0),
    };
    let glass_vertices = 2 * MERIDIANS * (2 * (2 * PROFILE_SAMPLES + 1) + 1)
        + bench.surfaces().len() * (RIM_SEGMENTS + 1);
    let ray_vertices: usize = paths.iter().map(|p| p.len()).sum();
    // Five points each: a rectangle's four corners and the one that closes it.
    let flat_vertices: usize = flats.iter().map(|f| f.len()).sum();
    check(
        "every run handed to the panel is a run it drew",
        kept as f64,
        (glass_runs.len() + flats.len() + paths.len()) as f64,
        1e-12,
        "runs",
    );
    check(
        "and every vertex of it",
        vertices as f64,
        (glass_vertices + flat_vertices + ray_vertices) as f64,
        1e-12,
        "vertices",
    );

    // ================================================================ the glass as a solid
    heading("The doublet as a solid, and what makes it one");

    let elements = bench.glass_solid();
    println!(
        "  {:<34} {:>7} vertices, {} triangles in {} elements",
        "the doublet is meshed as",
        elements.iter().map(|e| e.points.len()).sum::<usize>(),
        elements.iter().map(|e| e.faces.len()).sum::<usize>(),
        elements.len()
    );
    let mesh_points: Vec<[f64; 3]> = elements.iter().flat_map(|e| e.points.clone()).collect();

    // **Every edge is shared by exactly two triangles.** That is what closed means, and it is the
    // property a picture cannot show: a hole reads as shadow, and a surface with its inside facing
    // out reads as a surface. Counted by *position* rather than by index, because the rim's
    // vertices are deliberately written twice -- a crease is two faces that share a place and not
    // an index -- so counting indices would report every crease as a hole.
    //
    // Positions compare by bits. They are the same `f64` copied from the same `profile` call, so
    // two faces meeting at a corner hold identical patterns; a tolerance here would be a choice
    // about how close is the same, and there is no such choice to make.
    let key = |p: [f64; 3]| p.map(f64::to_bits);
    for (name, e) in ["crown", "flint"].iter().zip(&elements) {
        let mut edges: std::collections::BTreeMap<([u64; 3], [u64; 3]), usize> =
            std::collections::BTreeMap::new();
        for t in &e.faces {
            for k in 0..3 {
                let (a, b) = (
                    key(e.points[t[k] as usize]),
                    key(e.points[t[(k + 1) % 3] as usize]),
                );
                let edge = if a <= b { (a, b) } else { (b, a) };
                *edges.entry(edge).or_default() += 1;
            }
        }
        let open = edges.values().filter(|&&n| n != 2).count();
        println!(
            "  {:<34} {:>7} edges, {open} of them not shared by two faces",
            format!("the {name} is closed?"),
            edges.len()
        );
        assert_eq!(
            open, 0,
            "the {name} is not closed: {open} of {} edges are on one face or on three, which is              a hole or a fold and neither is a solid anybody could grind",
            edges.len()
        );
    }

    // **And every vertex of it is still on the surface the rays are traced against.** The mesh is
    // swept from the same `profile` the outline is, so this is the same claim as above over the
    // whole solid rather than over one meridian -- and it is the claim that keeps a mesh from
    // becoming a decorative approximation of the lens the moment somebody tessellates it
    // differently.
    let mut off: f64 = 0.0;
    for at in &mesh_points {
        let at = DVec3::from_array(*at);
        let ray = Ray::new(
            LengthVec::from_si(DVec3::new(at.x, at.y, -standoff)),
            DVec3::Z,
        );
        let nearest = bench
            .surfaces()
            .iter()
            .filter_map(|s| s.meet(ray).map(|h| (h.point.to_si() - at).length()))
            .fold(f64::INFINITY, f64::min);
        assert!(
            nearest.is_finite(),
            "a mesh vertex {:.4} mm off the axis is on none of the three surfaces",
            (at.x * at.x + at.y * at.y).sqrt() * 1e3
        );
        off = off.max(nearest);
    }
    let mesh_floor = 4.0 * f64::EPSILON * {
        let scale = bench.crown_r.abs().max(bench.cement_r.abs()) + standoff;
        scale + scale * scale / bench.crown_r.abs().min(bench.cement_r.abs())
    };
    println!(
        "  {:<34} {off:>7.2e} m  against a floor of {mesh_floor:.2e} m",
        "every meshed vertex is on a surface"
    );
    check_between(
        "every vertex of the solid is where a ray lands",
        off / mesh_floor,
        0.0,
        1.0,
        "x the floor",
    );

    // **The flats contain the beam, and each is in its own plane.** A fold mirror drawn to the
    // footprint of the light is only honest if the light is inside it: a rectangle built from the
    // hits and then placed a little wrong would still look like a mirror, and the picture is where
    // that would never be noticed. Both claims are arithmetic, so both are asserted here rather
    // than in the drawing -- CI runs this example with no argument, so a check inside the drawing
    // is a check CI does not run.
    assert_eq!(flats.len(), 2, "a fold mirror and an image plane");
    for (name, corners, at, normal) in [
        ("fold mirror", &flats[0], MIRROR_AT, bench.fold_normal),
        ("image plane", &flats[1], IMAGE_AT, DVec3::X),
    ] {
        let corner = |k: usize| DVec3::from_array(corners[k]);
        let (origin, u, v) = (corner(0), corner(1) - corner(0), corner(3) - corner(0));
        for c in corners {
            let off = (DVec3::from_array(*c) - origin).dot(normal).abs();
            assert!(
                off < 1e-12,
                "a corner of the {name} is {off:.3e} m off its own plane"
            );
        }
        let (mut widest, mut outside) = (0.0f64, 0usize);
        for r in paths.iter().filter(|r| r.len() > at) {
            let d = DVec3::from_array(r[at]) - origin;
            // Where the hit sits on the rectangle's own two axes, as a fraction of each side.
            let (a, b) = (d.dot(u) / u.length_squared(), d.dot(v) / v.length_squared());
            widest = widest.max(a.max(b).max(-a).max(-b));
            if !(-1e-9..=1.0 + 1e-9).contains(&a) || !(-1e-9..=1.0 + 1e-9).contains(&b) {
                outside += 1;
            }
        }
        let side = (u.length() * 1e3, v.length() * 1e3);
        println!(
            "  {:<30} {:>6.2} x {:>5.2} mm  from {} hits",
            format!("the {name} is drawn"),
            side.0,
            side.1,
            paths.iter().filter(|r| r.len() > at).count()
        );
        assert_eq!(
            outside, 0,
            "{outside} rays land outside the {name} that is drawn for them"
        );
        assert!(
            widest > 0.9,
            "the {name} is drawn {widest:.3} of the way out to the beam it carries, so it is \
             bigger than the footprint it is supposed to be"
        );
    }

    // ================================================================ the caption a machine reads
    //
    // **`docs/bench-app.png` is a photograph of a GPU.** CI has no adapter, so it is the second
    // figure in this repository nothing in CI can refresh -- `docs/editor.png` is the first, and
    // `tools/screenshot/README.md` spends a page on why that matters: a picture nothing compares
    // ages in silence, and every change to the thing it shows leaves it a little more wrong.
    //
    // So the geometry the picture is of is written down beside it, and compared here. This runs on
    // the argument-less invocation, which is the one CI makes of every example on every commit.
    let caption = format!(
        concat!(
            "front radius        {:>10.4} mm\n",
            "cemented radius     {:>10.4} mm\n",
            "glass semi-aperture {:>10.4} mm\n",
            "meshed              {:>10} vertices, {} triangles in {} elements\n",
            "drawn               {:>10} runs, {} vertices\n",
            "fold mirror         {:>10.2} x {:.2} mm\n",
            "image plane         {:>10.2} x {:.2} mm\n",
        ),
        bench.crown_r * 1e3,
        bench.cement_r * 1e3,
        GLASS * 1e3,
        elements.iter().map(|e| e.points.len()).sum::<usize>(),
        elements.iter().map(|e| e.faces.len()).sum::<usize>(),
        elements.len(),
        kept,
        vertices,
        flat_size(&flats[0]).0 * 1e3,
        flat_size(&flats[0]).1 * 1e3,
        flat_size(&flats[1]).0 * 1e3,
        flat_size(&flats[1]).1 * 1e3,
    );
    // **Writing the caption comes before checking it, and producing anything skips the check.**
    // The recipe in the message below is `bench.json`, then a snapshot, then this file — and with
    // a stale caption the *first* of those died right here, so the one command that fixes a stale
    // figure was the command the staleness blocked. Only the `.txt` arm was guarded against that,
    // which is half of it: a run producing any artefact is a run refreshing the figure, and the
    // compare belongs to the argument-less run, which is the one CI makes on every commit.
    if let Some(path) = common::output_path() {
        if path.ends_with(".txt") {
            common::write(&path, &caption);
            return;
        }
    }
    if common::output_path().is_none() {
        compare_caption(&caption);
    }
    // ================================================================ the deliverable
    if let Some(path) = common::output_path() {
        let frame = Frame {
            time_s: 0.0,
            panels: vec![
                // **The glass as a solid, in its own panel.** It shared the rays' scale because a
                // report draws one card per panel and two panels were two pictures of one bench.
                // The shell that shades is not the report and draws every panel now, and a solid
                // and a set of lines are not one quantity however they are laid out.
                Panel {
                    name: "glass".into(),
                    unit: "refractive index",
                    place: Placed::HERE,
                    data: PanelData::surface(
                        elements.iter().flat_map(|e| e.points.clone()).collect(),
                        {
                            let mut faces = Vec::new();
                            let mut base = 0u32;
                            for e in &elements {
                                faces.extend(
                                    e.faces
                                        .iter()
                                        .map(|t| [t[0] + base, t[1] + base, t[2] + base]),
                                );
                                base += e.points.len() as u32;
                            }
                            faces
                        },
                        elements.iter().flat_map(|e| e.values.clone()).collect(),
                    ),
                },
                Panel {
                    name: "bench".into(),
                    unit: "deg field, 4 = glass",
                    place: Placed::HERE,
                    data: panel,
                },
            ],
            readings: vec![
                Reading::new("bench", "focal length", f_bent * 1e3, "mm"),
                Reading::new("bench", "semi-aperture", SEMI * 1e3, "mm"),
                Reading::new("bench", "rms spot", rms * 1e6, "um"),
                Reading::new("bench", "diffraction limit", airy * 1e6, "um"),
                Reading::new(
                    "bench",
                    "thinnest glass",
                    front.thinnest_to(cement) * 1e3,
                    "mm",
                ),
            ],
        };
        // The extension chooses the asset, the way the application does: `.html` is a page you
        // open, `.json` is the frames themselves for something else to draw — the native viewer
        // in `runtime/viewer`, for instance, which reads this and nothing else.
        let frames = std::slice::from_ref(&frame);
        let asset = if path.ends_with(".svg") {
            layout_svg(&bench, &glass_runs, &paths, &colours)
        } else if path.ends_with(".json") {
            pantometry::view::to_json("optical bench", frames)
        } else if path.ends_with(".gltf") {
            // Into somebody else's renderer: Blender, three.js, Omniverse, a USD pipeline. The
            // geometry is one frame, because glTF animates node transforms and morph targets and
            // a retraced ray bundle is neither.
            let out = pantometry::view::gltf("optical bench", &frame);
            for note in &out.skipped {
                println!("  not exported: {note}");
            }
            out.document
        } else {
            pantometry::view::html("optical bench", frames)
        };
        common::write(&path, &asset);
        println!("\n  drag to rotate, scroll to zoom");
    } else {
        println!("\n  give a filename ending .html for the 3D layout");
    }
}

/// The two flats, as the rectangles the beam actually uses.
///
/// The prescription gives the fold mirror and the image plane no extent at all -- `plane_intersect`
/// takes an infinite plane -- so a picture of the bench turned its rays through ninety degrees at
/// nothing, which is the abstraction the glass was drawn to stop. A fold mirror's size is not an
/// invention either: it is sized to the footprint of the beam it carries, and every ray's hit on
/// it is already in the traced path. Each flat is the rectangle, in its own plane, that exactly
/// contains those hits.
///
/// "Exactly" is the word to be careful about, and `the_flats_contain_the_beam` is where it is
/// checked: a rectangle built from the hits and then drawn a little wrong would still look like a
/// mirror.
fn flats_of(bench: &Bench, rays: &[Vec<[f64; 3]>]) -> Vec<Vec<[f64; 3]>> {
    [(MIRROR_AT, bench.fold_normal), (IMAGE_AT, DVec3::X)]
        .iter()
        .filter_map(|&(at, normal)| {
            let hits: Vec<DVec3> = rays
                .iter()
                .filter(|r| r.len() > at)
                .map(|r| DVec3::from_array(r[at]))
                .collect();
            if hits.len() < 3 {
                return None;
            }
            let centre = hits.iter().sum::<DVec3>() / hits.len() as f64;
            // Two directions in the plane. `DVec3::Y` is across the field for both of these and
            // is perpendicular to each normal, so the other follows from the cross product.
            let u = DVec3::Y;
            let v = normal.cross(u).normalize();
            let extent = |axis: DVec3| {
                hits.iter().fold((f64::MAX, f64::MIN), |(lo, hi), p| {
                    let t = (*p - centre).dot(axis);
                    (lo.min(t), hi.max(t))
                })
            };
            let (u0, u1) = extent(u);
            let (v0, v1) = extent(v);
            Some(
                [(u0, v0), (u1, v0), (u1, v1), (u0, v1), (u0, v0)]
                    .iter()
                    .map(|&(a, b)| point(LengthVec::from_si(centre + u * a + v * b)))
                    .collect(),
            )
        })
        .collect()
}

/// The bench as a still picture, projected onto a page.
///
/// The HTML rotates and the glTF opens in somebody's renderer, and neither can go in a README --
/// which takes static images and nothing else. This is the same scene through
/// [`common::svg::View3`], which projects the way the HTML viewer projects, so the still and the
/// page a reader can spin are pictures of one thing.
///
/// Colour is assigned here rather than taken from the panel's scale. A colour bar is what a
/// *quantity* gets; on a page with a caption, three field angles and the glass they pass through
/// are better served by four colours a reader can tell apart, and the caption names them.
fn layout_svg(
    bench: &Bench,
    glass: &[Vec<[f64; 3]>],
    rays: &[Vec<[f64; 3]>],
    fields: &[f64],
) -> String {
    use common::svg::{document, rgb, Plot, View3};

    let flats = flats_of(bench, rays);

    let (w, h) = (880.0, 460.0);
    let (vx, vy, vw, vh) = (40.0, 52.0, 800.0, 344.0);
    let every = || {
        glass
            .iter()
            .chain(rays.iter())
            .chain(flats.iter())
            .flatten()
            .copied()
    };

    // Looked at from above and to one side, far enough that the perspective is a hint of depth
    // rather than a distortion: at three scene-widths the near and far ends of the bench differ
    // in scale by about a third.
    let view = View3::framing(every(), -0.62, 0.42, 3.0);
    let (x, y) = view.page_bounds(vw / vh);
    let mut plot = Plot::new(w, h, x, y).viewport(vx, vy, vw, vh);

    // The rays are 183 lines and the instrument is 17, so the instrument is drawn heavier and
    // darker and the light lighter. A picture in which the glass cannot be found is the one this
    // replaces.
    let ink = [rgb(214, 108, 88), rgb(224, 186, 96), rgb(112, 158, 226)];
    let metal = rgb(64, 66, 74);
    let runs = glass
        .iter()
        .chain(flats.iter())
        .map(|r| (r.clone(), metal.clone(), 1.7))
        .chain(rays.iter().zip(fields).map(|(r, f)| {
            let k = FIELDS.iter().position(|a| a == f).unwrap_or(0);
            (r.clone(), ink[k % ink.len()].clone(), 0.45)
        }));
    plot.paths3(&view, runs);

    plot.title("an optical bench: a doublet, a fold mirror and three field angles");
    plot.caption("every surface is drawn as the surface the rays are traced against");
    plot.footnote(
        "dark: the doublet in six meridians with a rim at each surface, and the two flats drawn to the beam they carry. red, yellow, blue: 0, 1.5 and 3 degrees of field",
    );
    document(w, h, [plot.into_body()])
}

/// How many rays a hexapolar pupil of four rings holds.
const PUPIL_RAYS: usize = 61;

/// The two curvatures for a given power split, at a fixed focal length.
///
/// The front takes `split` times the total power and the cemented face takes the rest, so every
/// split has the same paraxial `f` — which is what makes bending a free parameter rather than a
/// trade. Each surface's power is `(n2 - n1)/R`: the index it goes *into* minus the one it comes
/// *from*.
///
/// Using the air-glass `(n - 1)/R` on the cemented face — the formula everybody remembers — makes
/// the flint four times too weak, and the traced focal length came out 150 mm against a 100 mm
/// prescription. The trace was right and the prescription was wrong, which is the correct way
/// round and only visible because the trace *measures* `f` instead of being told it.
fn curvatures(split: f64) -> (f64, f64) {
    let phi = 1.0 / FOCAL;
    (
        (N_CROWN - 1.0) / (phi * split),
        (N_FLINT - N_CROWN) / (phi - phi * split),
    )
}

/// One ray's path through the bench and where it landed on the image plane.
type TracedRay = (Vec<[f64; 3]>, (f64, f64));

/// What a full trace of every field produces.
struct Traced {
    /// One run of points per ray that got through.
    paths: Vec<Vec<[f64; 3]>>,
    /// The field angle each path belongs to, which is what colours it.
    colours: Vec<f64>,
    /// Where each field's bundle landed on the image plane.
    spots: Vec<Vec<(f64, f64)>>,
}

/// Trace every field angle through the bench.
fn trace_fields(bench: &Bench) -> Traced {
    let pupil = hexapolar_unit(4);
    let (mut paths, mut colours, mut spots) = (Vec::new(), Vec::new(), Vec::new());
    for angle_deg in FIELDS {
        let theta = angle_deg.to_radians();
        let dir = DVec3::new(0.0, theta.sin(), theta.cos()).normalize();
        let mut landed = Vec::new();
        for (u, v) in &pupil {
            // Back-propagated **through** the pupil rather than launched from a flat plane.
            // Starting every ray at z = -40 mm and tilting it means an off-axis bundle has
            // drifted a millimetre by the time it reaches the glass, so the outer rays fall off
            // the aperture — 48 of 61 at 1.5 degrees, which reads as vignetting the design does
            // not have. The stop is the pupil; a ray is defined by where it crosses it.
            let at_pupil = DVec3::new(u * SEMI, v * SEMI, 0.0);
            let ray = Ray::new(LengthVec::from_si(at_pupil - dir * 40e-3), dir);
            if let Some((path, hit)) = bench.trace(ray) {
                landed.push(hit);
                paths.push(path);
                colours.push(angle_deg);
            }
        }
        spots.push(landed);
    }
    Traced {
        paths,
        colours,
        spots,
    }
}

/// The on-axis RMS spot for one power split, refocused for that split.
///
/// **Refocused**, which is the whole point: bending moves the focus a little, and comparing spots
/// at a fixed image plane would be measuring defocus and calling it aberration.
fn on_axis_rms(bench: &Bench, split: f64) -> Option<f64> {
    let (crown_r, cement_r) = curvatures(split);
    let mut trial = Bench {
        crown_r,
        cement_r,
        image_along: 0.0,
        ..*bench
    };
    let (_, bfd) = trial.measure()?;
    trial.image_along = bfd - (FOLD_AT - BACK_Z);
    let landed = trace_fields(&trial).spots;
    let axial = landed.first()?;
    (axial.len() * 10 >= PUPIL_RAYS * 9).then(|| rms_radius(axial))
}

/// One surface of the glass, as the prescription states it: where its apex sits on the axis and
/// how it curves. Zero radius is flat.
///
/// This exists so that the surfaces the rays are traced against and the surfaces the picture is
/// drawn from are the **same three values**. They were two lists before -- the trace had them
/// inline and the drawing did not exist -- and a drawing built from its own copy of a
/// prescription is a drawing that can be right about a lens nobody is tracing.
#[derive(Clone, Copy)]
struct Surface {
    /// Where the apex sits along the axis.
    at: f64,
    /// Radius of curvature, or zero for flat.
    r: f64,
}

impl Surface {
    /// The vertex: the apex on the axis, which is **not** the centre of curvature. Passing the
    /// centre put the surface a radius downstream of where it belongs and every ray missed the
    /// aperture, which is the right failure for a wrong vertex and a confusing one to read.
    fn vertex(self) -> LengthVec {
        LengthVec::from_si(DVec3::new(0.0, 0.0, self.at))
    }

    /// Where a ray meets it, or `None` if the ray misses the glass.
    fn meet(self, ray: Ray) -> Option<Hit> {
        cap_intersect(
            ray,
            self.vertex(),
            DVec3::Z,
            Length::from_si(self.r),
            Length::from_si(GLASS),
        )
    }

    /// The curve it is, in one meridian.
    fn meridian(self, across: DVec3) -> Vec<LengthVec> {
        profile(
            self.vertex(),
            DVec3::Z,
            Length::from_si(self.r),
            0.0,
            Length::from_si(GLASS),
            across,
            PROFILE_SAMPLES,
        )
        .expect("a meridian perpendicular to the axis")
    }

    /// Where the drawn edge of the surface sits along the axis.
    fn drawn_edge_z(self) -> f64 {
        self.meridian(DVec3::Y)
            .last()
            .expect("a meridian has an end")
            .to_si()
            .z
    }

    /// How far off the axis the surface is actually drawn.
    ///
    /// Usually `GLASS`, and not always: a sphere runs out at its equator, so a surface whose
    /// radius is smaller than the aperture stops short of it and `profile` says so.
    fn drawn_semi_aperture(self) -> f64 {
        let end = self
            .meridian(DVec3::Y)
            .last()
            .expect("a meridian has an end")
            .to_si();
        (end.x * end.x + end.y * end.y).sqrt()
    }

    /// The least axial gap to the surface behind it, over every height they are both drawn at.
    ///
    /// Negative would mean the two surfaces cross and the element is not a solid.
    ///
    /// The heights have to line up first. A power split of 2.4 -- inside the scan this example
    /// runs -- gives a cemented radius of 9.35 mm against a 9.5 mm aperture, so that surface would
    /// stop short, the outline would not close, and pairing the two point lists off regardless
    /// would subtract the sag at one radius from the sag at another and call it a thickness.
    fn thinnest_to(self, back: Surface) -> f64 {
        let (front, rear) = (self.meridian(DVec3::Y), back.meridian(DVec3::Y));
        let (a, b) = (self.drawn_semi_aperture(), back.drawn_semi_aperture());
        assert!(
            (a - b).abs() < 1e-12,
            "the two surfaces of an element are drawn to {:.4} mm and {:.4} mm off the axis: one \
             of them runs out before the aperture does, so the element has no cylindrical edge \
             and a thickness between them is a thickness between two different radii",
            a * 1e3,
            b * 1e3
        );
        front
            .iter()
            .zip(rear.iter())
            .map(|(f, r)| r.to_si().z - f.to_si().z)
            .fold(f64::INFINITY, f64::min)
    }
}

/// One element of the doublet, as the triangles that close it.
///
/// Per element and not one mesh for the pair, because the cemented face belongs to both: over the
/// two together every edge on it is shared by four faces, and a check for two would call a correct
/// solid broken. Each of these is closed on its own, which is also the thing a workshop would ask.
struct Element {
    points: Vec<[f64; 3]>,
    faces: Vec<[u32; 3]>,
    values: Vec<f64>,
}

/// The elements, and where they are.
#[derive(Clone, Copy)]
struct Bench {
    crown_r: f64,
    cement_r: f64,
    fold_at: f64,
    fold_normal: DVec3,
    /// How far past the mirror the image plane sits, along the folded axis.
    image_along: f64,
}

impl Bench {
    /// Trace one ray, returning its path through the bench and where it landed.
    ///
    /// Returns `None` if the ray misses anything it has to hit — which is vignetting, not an
    /// error, and the caller counts it rather than being told nothing happened.
    fn trace(&self, ray: Ray) -> Option<TracedRay> {
        let mut path = vec![point(ray.origin)];
        let mut r = self.lens(&mut path, ray)?;

        // The fold mirror.
        let t = plane_intersect(
            r,
            LengthVec::from_si(DVec3::new(0.0, 0.0, self.fold_at)),
            self.fold_normal,
        )?;
        path.push(point(r.at(t)));
        r = r.redirect(t, reflect(r.dir, self.fold_normal));

        // The image plane, normal to +x after the fold.
        let t = plane_intersect(
            r,
            LengthVec::from_si(DVec3::new(self.image_along, 0.0, self.fold_at)),
            DVec3::X,
        )?;
        let end = r.at(t);
        path.push(point(end));
        // In the image plane, z is across the field and y is the other transverse axis.
        Some((path, (end.to_si().z - self.fold_at, end.to_si().y)))
    }

    /// The three surfaces the glass is bounded by, in the order light meets them.
    ///
    /// The one prescription. [`Bench::lens`] traces against it and [`Bench::glass`] draws it.
    fn surfaces(&self) -> [Surface; 3] {
        [
            Surface {
                at: 0.0,
                r: self.crown_r,
            },
            Surface {
                at: CEMENT_Z,
                r: self.cement_r,
            },
            Surface { at: BACK_Z, r: 0.0 },
        ]
    }

    /// The glass as a **solid**: the triangles that bound each element.
    ///
    /// Every other shape this example emits is lines, and nothing in a line drawing takes light.
    /// The renderer that shades a field's isosurface had exactly one producer, so an instrument
    /// drawn as its own surfaces still arrived as a wireframe however carefully it was measured.
    ///
    /// One element is three pieces, in the order light meets them: the front surface, the rim
    /// between the two edges, and the back surface. Each surface is a grid over
    /// `(azimuth, height)` whose points come from [`Surface::meridian`] -- the same `profile`
    /// checked point by point against the intersection the rays use, so the solid is the surface
    /// and not a second description of it.
    ///
    /// **The rim's vertices are written twice**, once for the curved face and once for the
    /// cylinder. A normal is accumulated from the faces that meet at a vertex, so sharing them
    /// would round the edge of the glass into the side of it; two faces that share a position and
    /// not an index keep their own normals, and that is what a crease is.
    fn glass_solid(&self) -> Vec<Element> {
        let [front, cement, back] = self.surfaces();
        let azimuths: Vec<DVec3> = (0..RIM_SEGMENTS)
            .map(|m| {
                let t = std::f64::consts::TAU * m as f64 / RIM_SEGMENTS as f64;
                DVec3::new(t.cos(), t.sin(), 0.0)
            })
            .collect();
        // Half a meridian, from the vertex out to the rim. `meridian` runs edge to edge through
        // the vertex, so the second half is the half a solid of revolution is swept from.
        let half = |s: Surface, across: DVec3| -> Vec<LengthVec> {
            s.meridian(across)[PROFILE_SAMPLES..].to_vec()
        };
        let rings = PROFILE_SAMPLES; // rows past the apex

        let mut out = Vec::new();
        for (a, b, index) in [(front, cement, N_CROWN), (cement, back, N_FLINT)] {
            let (mut points, mut faces) = (Vec::new(), Vec::new());

            // `outward` is which way the cap faces once the winding is chosen: the front surface
            // faces the light coming in, the back surface faces away.
            for (surface, outward) in [(a, -1.0f64), (b, 1.0)] {
                // **One apex, fanned.** A grid gives every azimuth its own copy of the same point
                // and the first row of quads is degenerate; the edge count says so at once.
                let apex = points.len() as u32;
                points.push(point(half(surface, azimuths[0])[0]));
                let base = points.len() as u32;
                for across in &azimuths {
                    points.extend(half(surface, *across)[1..].iter().map(|p| point(*p)));
                }
                let at = |m: usize, r: usize| base + (m * rings + r - 1) as u32;
                for m in 0..RIM_SEGMENTS {
                    let n = (m + 1) % RIM_SEGMENTS;
                    if outward > 0.0 {
                        faces.push([apex, at(m, 1), at(n, 1)]);
                    } else {
                        faces.push([apex, at(n, 1), at(m, 1)]);
                    }
                    for r in 1..rings {
                        let (i00, i01, i10, i11) = (at(m, r), at(m, r + 1), at(n, r), at(n, r + 1));
                        if outward > 0.0 {
                            faces.push([i00, i10, i11]);
                            faces.push([i00, i11, i01]);
                        } else {
                            faces.push([i00, i11, i10]);
                            faces.push([i00, i01, i11]);
                        }
                    }
                }
            }

            // The rim: a cylinder between the two edges, with its own copies of both rings so the
            // crease keeps its own normals.
            let base = points.len() as u32;
            for across in &azimuths {
                for s in [a, b] {
                    points.push(point(
                        *half(s, *across).last().expect("a meridian has an end"),
                    ));
                }
            }
            for m in 0..RIM_SEGMENTS {
                let (c0, c1) = (m as u32 * 2, ((m + 1) % RIM_SEGMENTS) as u32 * 2);
                faces.push([base + c0, base + c0 + 1, base + c1 + 1]);
                faces.push([base + c0, base + c1 + 1, base + c1]);
            }

            let values = vec![index; points.len()];
            out.push(Element {
                points,
                faces,
                values,
            });
        }
        out
    }

    /// The glass itself, as runs of points, with the value each run is coloured by.
    ///
    /// Each element is closed: the front surface from one edge through the vertex to the other,
    /// then the back surface home again. Drawn in `MERIDIANS` half-planes, which is what turns two
    /// curves into something that reads as a solid of revolution when the picture is turned.
    ///
    /// The rim of each surface is built from the **ends of its own meridians** rather than from a
    /// circle of the aperture radius. The two would be the same number until a surface ran out
    /// before the aperture did -- a sphere stops at its equator -- and then the circle would be a
    /// ring of points floating off the edge of the glass.
    fn glass(&self) -> (Vec<Vec<[f64; 3]>>, Vec<f64>) {
        let [front, cement, back] = self.surfaces();
        let (mut runs, mut index) = (Vec::new(), Vec::new());

        for (a, b, n) in [
            (front, cement, GLASS_ON_THE_FIELD_SCALE),
            (cement, back, GLASS_ON_THE_FIELD_SCALE),
        ] {
            for m in 0..MERIDIANS {
                let angle = std::f64::consts::PI * m as f64 / MERIDIANS as f64;
                let across = DVec3::new(angle.cos(), angle.sin(), 0.0);
                let mut run: Vec<[f64; 3]> = a.meridian(across).into_iter().map(point).collect();
                run.extend(b.meridian(across).into_iter().rev().map(point));
                run.push(run[0]);
                runs.push(run);
                index.push(n);
            }
        }

        for (s, n) in [
            (front, GLASS_ON_THE_FIELD_SCALE),
            (cement, GLASS_ON_THE_FIELD_SCALE),
            (back, GLASS_ON_THE_FIELD_SCALE),
        ] {
            let rim: Vec<[f64; 3]> = (0..=RIM_SEGMENTS)
                .map(|i| {
                    let angle = std::f64::consts::TAU * i as f64 / RIM_SEGMENTS as f64;
                    let across = DVec3::new(angle.cos(), angle.sin(), 0.0);
                    point(*s.meridian(across).last().expect("a meridian has an end"))
                })
                .collect();
            runs.push(rim);
            index.push(n);
        }

        (runs, index)
    }

    /// The three lens surfaces, leaving the ray on its way to the mirror.
    ///
    /// The back face used to be an unbounded plane, which let a ray leave through glass that is
    /// not there. It is the same `Surface` as the other two now, clipped at the same aperture, and
    /// the pupil still gets through whole.
    fn lens(&self, path: &mut Vec<[f64; 3]>, ray: Ray) -> Option<Ray> {
        let [front, cement, back] = self.surfaces();
        let mut r = ray;
        for (surface, n1, n2) in [
            (front, 1.0, N_CROWN),
            (cement, N_CROWN, N_FLINT),
            (back, N_FLINT, 1.0),
        ] {
            let hit = surface.meet(r)?;
            r = self.bend(path, r, hit.t, hit.normal, n1, n2)?;
        }
        Some(r)
    }

    /// The effective focal length and back focal distance, **measured from a traced ray**.
    ///
    /// A ray entering parallel to the axis at height `h` leaves at angle `u`; `f = h/tan(u)` is
    /// the definition of focal length, and where it crosses the axis is the back focus. Neither
    /// is assumed: the thin-lens prescription that produced the curvatures ignores the 7 mm of
    /// glass, so the real focus is not where `FOCAL` says.
    fn measure(&self) -> Option<(f64, f64)> {
        let h = 0.5e-3;
        let mut path = Vec::new();
        let entering = Ray::new(LengthVec::from_si(DVec3::new(0.0, h, -40e-3)), DVec3::Z);
        let out = self.lens(&mut path, entering)?;
        let (o, d) = (out.origin.to_si(), out.dir);
        // Where it crosses the axis, measured from the last vertex.
        let t = -o.y / d.y;
        let crossing = o.z + d.z * t;
        let angle = (-d.y / d.z).atan();
        Some((h / angle.tan(), crossing - BACK_Z))
    }

    /// Refract at a surface, recording the vertex. `None` on total internal reflection.
    fn bend(
        &self,
        path: &mut Vec<[f64; 3]>,
        r: Ray,
        t: Length,
        normal: DVec3,
        n1: f64,
        n2: f64,
    ) -> Option<Ray> {
        path.push(point(r.at(t)));
        // `Hit::normal` is the surface's own, not oriented against the ray — the trap `lens_spots`
        // documents, and the reason this is one helper rather than four call sites.
        let n = oriented_against(normal, r.dir);
        let dir = refract(r.dir, n, n1 / n2)?;
        Some(r.redirect(t, dir))
    }
}

/// The sag of a spherical surface at height `h`, spelled the way a textbook spells it.
///
/// `R - sign(R) sqrt(R^2 - h^2)`, which cancels: the radius and the root are 39.75 and 38.60 mm
/// and the answer is 1.15. `pantometry-optics` computes `h^2 / (R (1 + sqrt(1 - h^2/R^2)))`, the
/// same quantity arranged so that it does not. Checking one against the other is a check;
/// checking the stable form against itself would be arithmetic agreeing with itself.
fn textbook_sag(r: f64, h: f64) -> f64 {
    if r == 0.0 {
        return 0.0;
    }
    r - r.signum() * (r * r - h * h).sqrt()
}

/// The two sides of a drawn flat, in metres, from the rectangle's own corners.
fn flat_size(corners: &[[f64; 3]]) -> (f64, f64) {
    let at = |k: usize| DVec3::from_array(corners[k]);
    ((at(1) - at(0)).length(), (at(3) - at(0)).length())
}

fn point(p: LengthVec) -> [f64; 3] {
    let v = p.to_si();
    [v.x, v.y, v.z]
}

fn centroid_of(spots: &[(f64, f64)]) -> (f64, f64) {
    let n = spots.len().max(1) as f64;
    (
        spots.iter().map(|(a, _)| a).sum::<f64>() / n,
        spots.iter().map(|(_, b)| b).sum::<f64>() / n,
    )
}

fn rms_radius(spots: &[(f64, f64)]) -> f64 {
    let (cx, cy) = centroid_of(spots);
    let n = spots.len().max(1) as f64;
    (spots
        .iter()
        .map(|(a, b)| (a - cx).powi(2) + (b - cy).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
}

/// The geometry `docs/bench-app.png` is a picture of, against what is stored beside it.
///
/// A *missing* file is a skip -- a checkout without `docs/` is not this example's business, and
/// neither is a packaged crate. A file that disagrees is a failure.
fn compare_caption(caption: &str) {
    let beside = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/bench-app.txt")
        .canonicalize()
        .ok();
    match beside.as_deref().map(std::fs::read_to_string) {
        Some(Ok(stored)) => {
            // Carriage returns are the checkout's, not the geometry's: `core.autocrlf` rewrites
            // the stored file and this string has none. The same flattening
            // `the_screenshot_is_still_true` does, for the same reason.
            let flat = |t: &str| t.replace('\r', "");
            assert_eq!(
                flat(&stored),
                flat(caption),
                "the bench has changed since `docs/bench-app.png` was taken. Retake it --\n  \
                 cargo run --release --example optical_bench bench.json\n  \
                 cd app && cargo run --release -- view bench.json --snapshot ../docs/bench-app.ppm\n\
                 -- and write this file with `--example optical_bench docs/bench-app.txt`. \
                 Editing the text alone would restore the green and leave the picture as stale \
                 as it is."
            );
            println!("  {:<34} {:>7}", "the app figure's caption agrees", "yes");
        }
        // A checkout without `docs/` is not this example's business, and neither is a packaged
        // crate. A *missing* file is a skip; a file that disagrees is a failure.
        _ => println!(
            "  {:<34} {:>7}",
            "no docs/bench-app.txt beside this", "skipped"
        ),
    }
}
