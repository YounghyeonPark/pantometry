//! A real lens, traced ray by ray in three dimensions, drawn as the surfaces it is, and checked
//! against the formulae an optical designer would use to sanity-check it.
//!
//! Every other example here evaluates a closed form. This one *traces*: a hexapolar bundle of
//! rays leaves a pupil, refracts at four glass surfaces, and lands on a plane, and what it
//! shows — a focal length, a spot size, a colour spread — is only ever the consequence of where
//! those rays went. Nothing below computes a spot diagram directly.
//!
//! The system is a **cemented achromatic doublet**: a positive crown element and a negative
//! flint one, chosen so their dispersions oppose and the red and blue foci nearly coincide. It
//! is the oldest interesting design in optics and it is still what a first lens looks like.
//!
//! ```sh
//! cargo run --release --example lens_spots
//! cargo run --release --example lens_spots out.svg
//! ```
//!
//! Three checks, each against something the trace did not compute:
//!
//! - the **paraxial focal length** from the lensmaker's equation, against where the marginal ray
//!   actually crosses the axis;
//! - **spherical aberration**, which must make the marginal focus *shorter* than the paraxial one
//!   for a positive singlet, and by a smaller margin once the doublet is bent to correct it;
//! - the **achromatic condition**, that the F and C foci land closer together than either does
//!   to the d focus of the equivalent singlet.

mod common;
use common::svg::{document, rgb, ticks, Plot};
use common::{check_between, heading};

use glam::DVec3;
use pantometry::prelude::*;
use pantometry_core::oriented_against;
use pantometry_optics::geometry::{cap_intersect, hexapolar_unit, profile, refract, Hit, Ray};
use pantometry_optics::material::Material;

/// One refracting surface: a sphere of a given curvature with glass behind it.
struct Surface {
    /// Vertex position on the axis, in metres.
    z: f64,
    /// Radius of curvature, metres. Positive means the centre is to the right of the vertex.
    /// `f64::INFINITY` is a flat surface.
    radius: f64,
    /// Clear semi-aperture.
    semi: f64,
    /// What is *after* this surface, going left to right.
    after: Material,
}

impl Surface {
    /// The vertex, on the axis.
    fn vertex(&self) -> LengthVec {
        LengthVec::m(0.0, 0.0, self.z)
    }

    /// Radius as the geometry wants it: zero is flat, where this type spells it `INFINITY`.
    fn curvature(&self) -> Length {
        Length::m(if self.radius.is_finite() {
            self.radius
        } else {
            0.0
        })
    }

    /// Where a ray meets it.
    fn meet(&self, ray: Ray) -> Option<Hit> {
        cap_intersect(
            ray,
            self.vertex(),
            DVec3::Z,
            self.curvature(),
            Length::m(self.semi),
        )
    }

    /// **The curve the surface is**, in the meridional plane, as `(z, height)` in millimetres.
    ///
    /// The picture drew each of these as a straight line at the vertex: a lens as three vertical
    /// strokes, which is a symbol for a lens and not a lens. `profile` samples the same sag
    /// [`cap_intersect`] solves against, so what is drawn is the surface the rays meet -- checked
    /// point by point in `main`, because a drawing that is merely plausible is the easiest wrong
    /// answer in optics.
    fn outline(&self) -> Vec<(f64, f64)> {
        profile(
            self.vertex(),
            DVec3::Z,
            self.curvature(),
            0.0,
            Length::m(self.semi),
            DVec3::X,
            PROFILE_SAMPLES,
        )
        .expect("a meridian perpendicular to the axis")
        .into_iter()
        .map(|p| {
            let v = p.to_si();
            (v.z * 1e3, v.x * 1e3)
        })
        .collect()
    }
}

/// Points per half-meridian of a drawn surface.
const PROFILE_SAMPLES: usize = 24;

/// A cemented achromatic doublet, **solved** rather than guessed, for a 100 mm focal length.
///
/// The first version of this had radii I asserted were "the classical crown-first shape". They
/// were not achromatic: the colour spread came out only 1.8x better than a singlet, and the
/// check that was supposed to confirm the design instead caught me making it up.
///
/// So the prescription is derived from the condition itself. For two thin elements in contact,
/// colour cancels when the powers are inversely proportional to the Abbe numbers:
///
/// ```text
///     phi_1 / V_1 + phi_2 / V_2 = 0        with  phi_1 + phi_2 = 1/f
///
///     phi_1 =  (1/f) V_1 / (V_1 - V_2)
///     phi_2 = -(1/f) V_2 / (V_1 - V_2)
/// ```
///
/// Both `V` come from the catalogue rather than from here, so the design follows the glass. The
/// crown's front radius is a free parameter — the *bending*, which trades spherical aberration
/// against nothing else — and the remaining two follow from the element powers.
fn doublet() -> Vec<Surface> {
    let crown = Material::from_catalog("N-BK7").expect("N-BK7 is in the catalogue");
    let flint = Material::from_catalog("F2").expect("F2 is in the catalogue");
    let (v1, v2) = (crown.abbe(), flint.abbe());
    let (n1, n2) = (crown.n_d(), flint.n_d());

    let power = 1.0 / 0.100;
    let phi1 = power * v1 / (v1 - v2);
    let phi2 = -power * v2 / (v1 - v2);

    // Bending: chosen so the crown is close to equiconvex, which is near the minimum-spherical
    // shape for a positive element in this configuration.
    let r1 = 0.0619;
    // phi_1 = (n1 - 1)(1/r1 - 1/r2)  ->  r2, the cemented surface.
    let r2 = 1.0 / (1.0 / r1 - phi1 / (n1 - 1.0));
    // phi_2 = (n2 - 1)(1/r2 - 1/r3)  ->  r3, the back.
    let r3 = 1.0 / (1.0 / r2 - phi2 / (n2 - 1.0));

    println!(
        "  solved from V = {v1:.1} and {v2:.1}:  R1 {:.2} mm, R2 {:.2} mm, R3 {:.2} mm",
        r1 * 1e3,
        r2 * 1e3,
        r3 * 1e3
    );
    vec![
        Surface {
            z: 0.0,
            radius: r1,
            semi: 0.010,
            after: crown,
        },
        Surface {
            z: 0.0040,
            radius: r2,
            semi: 0.010,
            after: flint,
        },
        Surface {
            z: 0.0065,
            radius: r3,
            semi: 0.010,
            after: Material::air(),
        },
    ]
}

/// A single crown element of roughly the same power, for comparison.
///
/// Equiconvex N-BK7, bent to land near the doublet's focal length so the two are compared at the
/// same job rather than at the same radii.
fn singlet_like(_doublet: &[Surface]) -> Vec<Surface> {
    let crown = Material::from_catalog("N-BK7").expect("N-BK7 is in the catalogue");
    // (n-1)(1/r - 1/-r) = 2(n-1)/r = 1/f
    let r = 2.0 * (crown.n_d() - 1.0) * 0.100;
    vec![
        Surface {
            z: 0.0,
            radius: r,
            semi: 0.010,
            after: crown,
        },
        Surface {
            z: 0.004,
            radius: -r,
            semi: 0.010,
            after: Material::air(),
        },
    ]
}

/// Trace one ray through the stack and return where it crosses the axis, in metres from the
/// last vertex.
///
/// `None` if it misses an aperture or is totally internally reflected — which is a real outcome
/// and not an error, so it is reported as absence rather than as a panic.
fn focus_of(height: f64, wavelength: Length, surfaces: &[Surface]) -> Option<f64> {
    traced(height, wavelength, surfaces).map(|(focus, _)| focus)
}

/// The same trace, keeping the points it passes through.
///
/// `focus_of` computed every one of these and threw them away, so the picture drew a ray as two
/// straight pieces with the glass jumped over -- "the bending inside is a few tenths of a
/// millimetre and drawing it would add ink and no information", which was true of a picture whose
/// glass was three vertical lines. Beside the real surfaces it is the information: it is where the
/// bending happens.
fn traced(height: f64, wavelength: Length, surfaces: &[Surface]) -> Option<(f64, Vec<(f64, f64)>)> {
    let mut ray = Ray::new(LengthVec::m(height, 0.0, -0.05), DVec3::new(0.0, 0.0, 1.0));
    let mut before = Material::air();
    let mut path = vec![(-50.0, height * 1e3)];
    for s in surfaces {
        let Some(hit) = s.meet(ray) else {
            if std::env::var("TRACE").is_ok() {
                eprintln!("  MISS at z={}", s.z);
            }
            return None;
        };
        let eta = before.index(wavelength) / s.after.index(wavelength);
        // `Hit::normal` is the surface's own geometric normal and is **not** oriented against
        // the ray — the field's documentation says so, and points at this helper. Its sign
        // follows the curvature, so at the concave second surface it pointed downstream,
        // `refract` saw a ray leaving the glass it was entering, and the trace turned round and
        // walked back out of the lens. A surface cannot know which side you approached from,
        // so orienting is the caller's job and skipping it fails silently as a miss.
        let n = oriented_against(hit.normal, ray.dir);
        // The normal `cap_intersect` reports points back along the ray; `refract` wants it on
        // the incoming side, which is the same convention.
        let Some(dir) = refract(ray.dir, n, eta) else {
            if std::env::var("TRACE").is_ok() {
                eprintln!("  TIR at z={} eta={eta}", s.z);
            }
            return None;
        };
        if std::env::var("TRACE").is_ok() {
            eprintln!(
                "  z={} t={:.6} normal={:?} dir={:?}",
                s.z,
                hit.t.to_si(),
                hit.normal,
                dir
            );
        }
        let at = hit.point.to_si();
        path.push((at.z * 1e3, at.x * 1e3));
        ray = ray.redirect(hit.t, dir);
        before = s.after.clone();
    }
    // Where it crosses the axis: solve x + t·dx = 0 in the meridional plane.
    let p = ray.origin.to_si();
    if ray.dir.x.abs() < 1e-15 {
        return None;
    }
    let t = -p.x / ray.dir.x;
    let z = p.z + ray.dir.z * t;
    path.push((z * 1e3, 0.0));
    Some((z - surfaces.last()?.z, path))
}

/// The lensmaker's equation for a thin doublet: the sum of the elements' powers.
///
/// Thin, so it ignores the 6.5 mm of glass the trace does not — which is exactly why it is a
/// *check* and not the answer. Agreement to a per cent or so is what a designer expects between
/// a thin-lens estimate and a real trace.
fn thin_lens_focal(wavelength: Length, surfaces: &[Surface]) -> f64 {
    let mut power = 0.0;
    let mut before = Material::air().index(wavelength);
    for s in surfaces {
        let after = s.after.index(wavelength);
        power += (after - before) / s.radius;
        before = after;
    }
    1.0 / power
}

fn main() {
    let d = Length::nm(587.6); // the sodium d line, where focal lengths are quoted
    let f_line = Length::nm(486.1); // hydrogen F, blue
    let c_line = Length::nm(656.3); // hydrogen C, red
    let lens = doublet();

    heading("Where the rays say the focus is, against the lensmaker's equation");
    // The *paraxial* focus, from a ray so close to the axis that aberration cannot reach it.
    let paraxial = focus_of(1e-4, d, &lens).expect("a paraxial ray gets through");
    let thin = thin_lens_focal(d, &lens);
    check_between(
        "traced paraxial focal length",
        paraxial,
        thin * 0.95,
        thin * 1.05,
        "m",
    );
    println!(
        "    thin-lens estimate {:.4} m, traced {:.4} m — the thin form ignores 6.5 mm of glass",
        thin, paraxial
    );

    heading("Spherical aberration: which way each design bends the edge of the pupil");
    // A positive **singlet** always focuses its marginal ray short. That is a law, not a design
    // choice: the surfaces bend outer rays too strongly and nothing in a single element opposes
    // it. So the singlet carries the sign assertion.
    //
    // The **doublet** does not, and the first version of this example asserted it did. A cemented
    // achromat's spherical depends on its bending and can land either side of zero — this one
    // comes out *over*corrected. What is universal is that it has less of it, which is the
    // comparison below.
    let singlet = singlet_like(&lens);
    let s_para = focus_of(1e-4, d, &singlet).expect("paraxial through the singlet");
    let s_marg = focus_of(0.0095, d, &singlet).expect("marginal through the singlet");
    check_between(
        "singlet: marginal focus minus paraxial",
        s_marg - s_para,
        -0.01,
        0.0,
        "m",
    );

    let marginal = focus_of(0.0095, d, &lens).expect("the marginal ray gets through");
    let longitudinal = marginal - paraxial;
    println!(
        "    singlet {:+.3} mm, doublet {:+.3} mm — the singlet undercorrects, this doublet overcorrects",
        (s_marg - s_para) * 1e3,
        longitudinal * 1e3
    );
    check_between(
        "doublet spherical, as a fraction of the singlet's",
        longitudinal.abs() / (s_marg - s_para).abs(),
        0.0,
        0.6,
        "x",
    );

    heading("The achromat earns its name");
    // Crown and flint disperse oppositely, so F and C should land near each other. Compared
    // against the spread a single N-BK7 element of the same power would give.
    let (ff, fc) = (
        focus_of(1e-4, f_line, &lens).expect("blue gets through"),
        focus_of(1e-4, c_line, &lens).expect("red gets through"),
    );
    let doublet_spread = (ff - fc).abs();
    let singlet_spread = (focus_of(1e-4, f_line, &singlet).unwrap()
        - focus_of(1e-4, c_line, &singlet).unwrap())
    .abs();
    check_between(
        "doublet F-to-C spread, as a fraction of the singlet's",
        doublet_spread / singlet_spread,
        0.0,
        0.25,
        "x",
    );
    println!(
        "    doublet {:.3} mm against singlet {:.3} mm — {:.0}x tighter",
        doublet_spread * 1e3,
        singlet_spread * 1e3,
        singlet_spread / doublet_spread
    );

    heading("The spot the whole pupil makes, at the paraxial focus");
    let plane = lens.last().unwrap().z + paraxial;
    let mut radii = Vec::new();
    for (u, v) in hexapolar_unit(6) {
        let h = (u * u + v * v).sqrt() * 0.0095;
        if h < 1e-9 {
            continue;
        }
        if let Some(f) = focus_of(h, d, &lens) {
            // How far off axis this ray is when it reaches the paraxial plane.
            let z_focus = lens.last().unwrap().z + f;
            let miss = h * (plane - z_focus) / (z_focus - lens.last().unwrap().z).max(1e-12);
            radii.push(miss.abs());
        }
    }
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let rms = (radii.iter().map(|r| r * r).sum::<f64>() / radii.len() as f64).sqrt();
    println!(
        "  {} rays through, RMS spot radius {:.1} um",
        radii.len(),
        rms * 1e6
    );

    // A perfect lens of this aperture cannot beat the Airy radius, so a spot smaller than that
    // would mean the trace is wrong rather than the lens is good.
    let airy = 1.22 * d.to_si() * paraxial / (2.0 * 0.0095);
    println!(
        "    the diffraction limit for this aperture is {:.1} um",
        airy * 1e6
    );
    assert!(
        rms > airy * 0.1,
        "a geometric spot far under the Airy radius means the trace is not tracing"
    );

    // ================================================================ the drawn glass
    heading("The glass the picture draws, against the surfaces the rays are traced against");

    // **Every drawn point is a point a ray lands on.** `profile` evaluates a sag; `cap_intersect`
    // solves a quadratic. Agreeing means the outline on the front page is the surface the trace
    // uses, rather than a second and merely plausible description of it.
    //
    // The floor is the quadratic's own error: for an axial ray the rounding that matters is the
    // one in `R^2 - h^2`, which the square root divides by `2 sqrt(R^2 - h^2)`, so it is
    // `4 eps (S + S^2/|R|)` with `S = |R| + standoff`, per surface. `pantometry-optics` derives it
    // and measures nine configurations against it.
    let standoff = 0.05;
    let mut worst = 0.0f64;
    for s in &lens {
        let scale = s.radius.abs() + standoff;
        let floor = 4.0 * f64::EPSILON * (scale + scale * scale / s.radius.abs());
        for (z_mm, h_mm) in s.outline() {
            let at = DVec3::new(h_mm * 1e-3, 0.0, z_mm * 1e-3);
            let ray = Ray::new(
                LengthVec::m(at.x, 0.0, -standoff),
                DVec3::new(0.0, 0.0, 1.0),
            );
            let hit = s
                .meet(ray)
                .expect("a ray aimed at a drawn point meets the surface it was drawn from");
            worst = worst.max((hit.point.to_si() - at).length() / floor);
        }
    }
    println!(
        "  {:<34} {:>8} points on {} surfaces",
        "the doublet is drawn as",
        lens.len() * (2 * PROFILE_SAMPLES + 1),
        lens.len()
    );
    check_between(
        "every drawn point is where a ray lands",
        worst,
        0.0,
        1.0,
        "x the floor",
    );
    // Exactly zero everywhere would mean one of the two is reading the other.
    assert!(
        worst > 0.0,
        "every drawn point matched a traced point to the last bit, so the drawing and the trace \
         are not two ways of finding one surface"
    );

    // **The glass has thickness where it is drawn.** Two surfaces bounding one element can cross,
    // and a prescription whose elements cross is one no workshop can cut -- which the lensmaker's
    // equation above has no opinion on, because `f` is the same either way.
    for (name, pair) in ["crown", "flint"].iter().zip(lens.windows(2)) {
        let (front, back) = (pair[0].outline(), pair[1].outline());
        let thinnest = front
            .iter()
            .zip(back.iter())
            .map(|((zf, _), (zb, _))| zb - zf)
            .fold(f64::INFINITY, f64::min);
        println!("  {:<34} {thinnest:>8.3} mm at its thinnest", name);
        assert!(
            thinnest > 0.1,
            "the {name} is {thinnest:.4} mm thick at its thinnest, which is a prescription no \
             workshop would accept and, below zero, two surfaces that cross"
        );
    }

    if let Some(path) = std::env::args().nth(1) {
        let svg = draw(&lens, d, f_line, c_line);
        common::write(&path, &svg);
        println!("\n  wrote {path}");
    } else {
        println!("\n  give a filename to draw the ray fan");
    }
}

/// The meridional ray fan, which is the picture an optical designer actually reads.
fn draw(lens: &[Surface], d: Length, f: Length, c: Length) -> String {
    let z_end = (lens.last().unwrap().z + focus_of(1e-4, d, lens).unwrap()) * 1e3;
    let mut plot = Plot::new(880.0, 420.0, (-12.0, z_end * 1.06), (-12.0, 12.0))
        .viewport(64.0, 54.0, 800.0, 320.0);

    // Red, yellow, blue: the C, d and F lines, which is the order they focus in.
    for (wavelength, colour) in [
        (c, rgb(226, 78, 52)),
        (d, rgb(238, 196, 62)),
        (f, rgb(74, 132, 238)),
    ] {
        for k in -6..=6 {
            let h = k as f64 * 0.0095 / 6.0;
            if h.abs() < 1e-9 {
                continue;
            }
            // **The path the trace took**, not two straight pieces with the glass jumped over.
            // Every vertex here is a point `cap_intersect` returned, so the ray is drawn bending
            // where it bends -- which is only worth ink because the glass it bends in is drawn.
            let Some((_, path)) = traced(h, wavelength, lens) else {
                continue;
            };
            plot.polyline(
                path.into_iter().map(|(z, x)| (z.max(-12.0), x)),
                &colour,
                1.0,
            );
        }
    }

    // **The elements, as the surfaces they are.** This drew a straight segment at each vertex --
    // three vertical strokes for a doublet whose three surfaces are curved, which is a symbol for
    // a lens rather than a lens, and it sat on the front page. Each element is closed here: its
    // front surface from one edge through the vertex to the other, then the surface behind it
    // home again.
    for pair in lens.windows(2) {
        let mut run = pair[0].outline();
        run.extend(pair[1].outline().into_iter().rev());
        run.push(run[0]);
        plot.polyline(run, &rgb(120, 122, 132), 1.8);
    }
    plot.polyline(
        [(-12.0, 0.0), (z_end * 1.06, 0.0)],
        &rgb(150, 152, 162),
        0.8,
    );

    plot.axes(
        &ticks(-12.0, z_end * 1.06, 8),
        &ticks(-12.0, 12.0, 5),
        |v| format!("{v:.0}"),
        |v| format!("{v:.0}"),
    );
    plot.title("a cemented achromat, drawn as its surfaces and traced through them");
    plot.caption("z (mm) against ray height (mm) — red C, yellow d, blue F");
    plot.footnote(
        "the three colours cross the axis within a tenth of a millimetre of each other, which is what makes it an achromat",
    );
    document(880.0, 420.0, [plot.finish()])
}
