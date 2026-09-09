//! A patch of a face exposed, against the face entire.
//!
//! [`Solid3D::losing_from_within`] exposes a box on a face rather than the whole of it, because a
//! part is bolted at **pads**. A bracket carries heat from whatever is mounted on it, along its
//! own shape, into two or three bolt bosses; that path is the thermal question a shape answers,
//! and a part cooled over its entire footprint does not have one — every cell sheds where it
//! stands, the part falls together and the field is a constant.
//!
//! # Stating a smaller area is not the same thing, and that is what is measured here
//!
//! The stated area is divided among the cells on the face, so a tenth of the area is a tenth of
//! the conductance **spread over the whole face**. The right total in the wrong place. Two blocks
//! with identical total conductance — one over the face, one over a patch — must lose heat at the
//! same rate at a uniform temperature and must reach *different* fields once there is a gradient,
//! and both halves are asserted below.

use pantometry_core::{
    units::{Area, Energy, Length, Temperature, Time},
    Domain, Exchange, Substance,
};
use pantometry_thermal::{Environment, Face, Solid3D};

/// Aluminium with radiation switched off, so a balance is linear and a closed form is the
/// problem's rather than a linearisation's.
fn grey() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

/// A bar of `n` cells along x, one cell across y and z, of side `dx`.
fn bar(n: usize, dx: f64, initial_c: f64) -> Solid3D {
    Solid3D::new(
        "bar",
        grey(),
        (n, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(initial_c),
    )
}

fn film(h: f64, area: f64) -> Environment {
    Environment {
        ambient: Temperature::celsius(20.0),
        convection_w_per_m2_k: h,
        area: Area::from_si(area),
    }
}

/// March to steady state, checking that it arrived. Returns nothing; the caller reads the block.
fn settle(block: &mut Solid3D) {
    let mut bus = Exchange::new();
    let mut previous = block.mean_temperature().to_si();
    for round in 0..8000 {
        let mut window = 0.0;
        while window < 1.0 {
            let dt = block.max_stable_dt(Time::ZERO);
            block.step(Time::ZERO, dt, &mut bus).expect("a stable step");
            window += dt.to_si();
        }
        let now = block.mean_temperature().to_si();
        if round > 0 && (now - previous).abs() / window < 1e-7 {
            return;
        }
        previous = now;
    }
    panic!("it never settled");
}

/// **A patch sheds through the cells it names and no others.**
///
/// The load-bearing statement, and it is checkable without a solver: at a **uniform** temperature
/// every exposed cell sheds `h·a·(T − T∞)` with `a` its share, so the block's total loss is
/// `h·A·(T − T∞)` whatever the patch — the stated area is what is divided. What the patch decides
/// is *which* cells carry it.
///
/// So a plate cooled on a ninth of its face and one cooled on all of it, with the same stated
/// area, must lose heat at exactly the same rate for the first step from uniform — and the
/// **profile** must differ immediately after. A patch that changed the total would be a patch that
/// silently rescaled the boundary condition; one that changed nothing would not be a patch.
#[test]
fn a_patch_carries_the_stated_area_and_puts_it_where_it_says() {
    let build = |patched: bool| {
        let plate = Solid3D::new(
            "plate",
            grey(),
            (9, 9, 3),
            Length::mm(2.0),
            Temperature::celsius(200.0),
        );
        let area = 1e-4; // the same stated area either way
        if patched {
            // The middle three by three of the nine by nine face: a bolt pad.
            plate.losing_from_within(Face::ZMin, [3, 3, 6, 6], film(500.0, area))
        } else {
            plate.losing_from(Face::ZMin, film(500.0, area))
        }
    };
    let (mut whole, mut pad) = (build(false), build(true));
    assert_eq!(whole.patch_on(Face::ZMin), None);
    assert_eq!(pad.patch_on(Face::ZMin), Some([3, 3, 6, 6]));

    // One step from uniform: the same joules leave, because the stated area is the same.
    let dt = Time::from_si(pad.max_stable_dt(Time::ZERO).to_si() * 0.5);
    whole
        .step(Time::ZERO, dt, &mut Exchange::new())
        .expect("stable");
    pad.step(Time::ZERO, dt, &mut Exchange::new())
        .expect("stable");
    let (a, b) = (whole.lost_energy().to_si(), pad.lost_energy().to_si());
    println!("  first step: whole face lost {a:.9} J, the pad lost {b:.9} J");
    assert!(
        (a - b).abs() < 1e-12 * a.abs().max(1.0),
        "a patch changed the total conductance: {a:.12} against {b:.12}"
    );

    // And it left a different block behind. The pad's own cells are colder than the corner's;
    // over the whole face every cell fell together.
    let spread =
        |b: &Solid3D| b.temperature_at(0, 0, 0).to_si() - b.temperature_at(4, 4, 0).to_si();
    println!(
        "  after one step: whole face corner-to-centre {:.6} K, the pad {:.6} K",
        spread(&whole),
        spread(&pad)
    );
    assert!(
        spread(&whole).abs() < 1e-12,
        "a face cooled entire falls together: {:.3e} K",
        spread(&whole)
    );
    assert!(
        spread(&pad) > 1e-3,
        "a pad must cool its own cells first: {:.3e} K",
        spread(&pad)
    );
}

/// **A bar bolted at one end has the gradient its length says**, and cooled all over it has none.
///
/// The closed form is the one `a_block_that_generates` already uses: at steady state every watt
/// crosses every face and leaves through the film, so the rise from the source cell's centre is
/// `P·((n − ½)·dx/(k·A) + 1/(h·A))`. What is new is that the film is on a **patch** at the far
/// end of a bar whose other faces are insulated, which is what makes the conduction term the whole
/// point rather than a correction.
///
/// The comparison is the argument: the same bar with the same total conductance spread along its
/// long face instead is a **fin**, and its closed form is the other one every heat-transfer text
/// carries — `P/(kAm·tanh(mL))`, with `m = √(h′/kA)` and `h′` the loss per unit length per kelvin.
///
/// # Not a lump, which is what this was first written to assert
///
/// "Spread the conductance and the bar sits at the lumped balance `P/(hA)`" was measured and is
/// wrong: 263.37 K against 222.22. Heat generated at one end still has to reach the cells that
/// shed, and with only a sixteenth of the conductance under the source that path is most of the
/// answer. The fin equation is what describes it, and it agrees to **0.044%**.
#[test]
fn a_bar_bolted_at_one_end_has_a_gradient_and_one_cooled_all_over_does_not() {
    let (watts, n, dx, h) = (4.0, 16usize, 3e-3, 2000.0);
    let area = dx * dx;
    let k = grey()
        .thermal
        .expect("aluminium conducts")
        .conductivity
        .to_si();

    // Bolted at the far end: one cell of the x-max face, which is the whole face here, but the
    // long faces stay insulated — so the heat has sixteen cells to cross.
    let mut bolted = bar(n, dx, 20.0)
        .dissipating(watts, |i, _, _| i == 0)
        .losing_from(Face::XMax, film(h, area));
    settle(&mut bolted);
    let rise = bolted.temperature_at(0, 0, 0).to_si() - Temperature::celsius(20.0).to_si();
    let closed = watts * ((n as f64 - 0.5) * dx / (k * area) + 1.0 / (h * area));
    println!("  bolted at the end: {rise:.4} K, closed form {closed:.4} K");
    assert!(
        (rise / closed - 1.0).abs() < 2e-3,
        "a bar bolted at one end rose {rise:.4} K where the series says {closed:.4} K"
    );

    // The same stated conductance, spread along the bar's own long face instead of concentrated
    // at its end. Each of the `n` cells carries `h·A/n`, so the loss per unit length per kelvin
    // is `h·A/(n·dx)` and the bar is a fin with an insulated tip.
    let mut spread = bar(n, dx, 20.0)
        .dissipating(watts, |i, _, _| i == 0)
        .losing_from(Face::YMin, film(h, area));
    settle(&mut spread);
    let along = spread.temperature_at(0, 0, 0).to_si() - Temperature::celsius(20.0).to_si();
    let per_length = h * area / (n as f64 * dx);
    let m = (per_length / (k * area)).sqrt();
    let fin = watts / (k * area * m * (m * (n as f64 * dx)).tanh());
    println!("  cooled along its length: {along:.4} K, the fin says {fin:.4} K");
    // The tolerance is the discretisation: the fin equation is continuous and this is sixteen
    // cells, with the source and the tip each half a cell inside their own ends.
    assert!(
        (along / fin - 1.0).abs() < 5e-3,
        "a fin of {n} cells rose {along:.4} K where `P/(kAm tanh(mL))` says {fin:.4} K"
    );
    // And it is cheaper than bolting the same conductance to the far end, which is why a heatsink
    // has fins: the heat does not have to cross the whole bar to reach any of it.
    assert!(
        along < rise,
        "spreading the same conductance must cost less: {along:.4} against {rise:.4}"
    );
    // Cheaper, and *not* the lumped balance either — which is what this was first asserted to be.
    let lumped = watts / (h * area);
    assert!(
        along > 1.1 * lumped,
        "a fin is not a lump: {along:.4} K against the lumped {lumped:.4} K"
    );
}

/// **A patch on one face leaves the other five insulated**, and an empty one loses nothing.
///
/// The two ends. A patch that reached past its bounds would shed from cells the caller excluded —
/// the quiet direction, because the block still cools and still conserves — and a patch that
/// selects nothing must lose exactly zero rather than a very small amount, since a block that
/// cannot shed has no steady state and one that sheds a little has one far away.
#[test]
fn a_patch_reaches_no_further_than_its_bounds() {
    let build = |within: [usize; 4]| {
        Solid3D::new(
            "plate",
            grey(),
            (6, 6, 2),
            Length::mm(2.0),
            Temperature::celsius(200.0),
        )
        .losing_from_within(Face::ZMin, within, film(500.0, 1e-4))
    };

    // A corner cell only. After a while the far corner has cooled by conduction and the near one
    // by conduction *and* the film, so they cannot be equal — and the excluded corner must not
    // have shed anything of its own.
    let mut corner = build([0, 0, 1, 1]);
    let dt = Time::from_si(corner.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..400 {
        corner
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    let (near, far) = (
        corner.temperature_at(0, 0, 0).to_si(),
        corner.temperature_at(5, 5, 0).to_si(),
    );
    println!("  one-cell patch: the bolted corner {near:.4} K, the far one {far:.4} K");
    assert!(
        far - near > 1.0,
        "the patch cooled the whole face: {near:.4} against {far:.4}"
    );

    // A patch that selects nothing sheds nothing at all, to the bit.
    let mut shut = build([2, 2, 2, 2]);
    let start = shut.mean_temperature().to_si();
    shut.deposit(0, 0, 0, Energy::from_si(1.0));
    for _ in 0..400 {
        shut.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    println!(
        "  empty patch: {:.9} J lost, mean {start:.6} -> {:.6} K",
        shut.lost_energy().to_si(),
        shut.mean_temperature().to_si()
    );
    assert!(
        shut.lost_energy().to_si() == 0.0,
        "a patch that names no cells lost {:.6e} J",
        shut.lost_energy().to_si()
    );
}

/// **A patch tightens the stability limit, and the limit it reports is one the sweep survives.**
///
/// The dangerous direction, and it is not the one a sabotage finds first. Ignoring a patch when
/// the limit is computed makes the limit *shorter* than it needs to be — wasteful and safe, and
/// three tests above pass with that mistake in place, correctly.
///
/// What is not safe is a limit that comes out **longer**. A patch concentrates the stated area on
/// fewer cells, so each carries a larger share: nine cells of a bolt pad hold what six hundred
/// held between them. A limit computed from the face's count and a flux applied to the patch's
/// would be loose by that ratio — and an explicit boundary past its limit does not diverge loudly.
/// It oscillates about the air it is losing to while the conservation audit stays perfectly happy,
/// because the heat really did leave. `loss_conductance_at` and `film_flux` read one `cells_on`,
/// which `count_faces` fills with the patched count; this is what says so.
#[test]
fn a_patched_face_reports_a_limit_the_march_survives() {
    // A hard film on one cell of a big face: the worst ratio between a face's count and a
    // patch's that this grid can make.
    let mut plate = Solid3D::new(
        "plate",
        grey(),
        (12, 12, 2),
        Length::mm(2.0),
        Temperature::celsius(400.0),
    )
    .losing_from_within(Face::ZMin, [0, 0, 1, 1], film(8000.0, 5e-4));

    let unpatched = Solid3D::new(
        "plate",
        grey(),
        (12, 12, 2),
        Length::mm(2.0),
        Temperature::celsius(400.0),
    )
    .losing_from(Face::ZMin, film(8000.0, 5e-4));
    let (tight, loose) = (
        plate.max_stable_dt(Time::ZERO).to_si(),
        unpatched.max_stable_dt(Time::ZERO).to_si(),
    );
    println!("  limit with the pad {tight:.6e} s, over the whole face {loose:.6e} s");
    assert!(
        tight < loose,
        "concentrating the area on one cell must tighten the limit: {tight:.3e} against {loose:.3e}"
    );

    // And the reported limit is one the march survives. An over-long limit oscillates about
    // ambient rather than diverging, so what is asserted is that the cooled cell falls
    // monotonically — an oscillating boundary cell goes down, past, and back up.
    let mut previous = plate.temperature_at(0, 0, 0).to_si();
    for step in 0..600 {
        let dt = plate.max_stable_dt(Time::ZERO);
        plate
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("its own limit is a step it can take");
        let now = plate.temperature_at(0, 0, 0).to_si();
        assert!(
            now <= previous + 1e-9,
            "the bolted cell rose at step {step}: {previous:.6} then {now:.6} K"
        );
        previous = now;
    }
    println!("  600 steps at its own limit: the bolted cell fell to {previous:.4} K");
    assert!(previous < 400.0 + 273.15, "nothing cooled: {previous:.4} K");
}
