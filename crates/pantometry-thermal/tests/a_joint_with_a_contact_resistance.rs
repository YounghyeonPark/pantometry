//! A face that carries a stated resistance, against the series it is supposed to be.
//!
//! [`Solid3D::joined`] puts a contact conductance `h` on an interior face, in series with the two
//! half cells already there. The whole claim is one line of algebra:
//!
//! ```text
//! 1/k_face = 1/k_series + 1/(dx·h)      so      1/G_face = 1/G_series + 1/(A·h)
//! ```
//!
//! and that second form is what is checked here — the resistance a joint **adds** is exactly
//! `1/(A·h)`, whatever the grid, whatever the materials either side.
//!
//! # Why this is not a convenience over a thin layer of material
//!
//! A bolted joint, a thermal interface material, a solder layer, an oxide: all are tens of microns
//! in a part discretised at a millimetre. A layer of material is at least **one cell** thick, so
//! the thinnest resistance a grid can state at all is `dx/k` — and a 1.5 mm cell cannot say 100 µm
//! of solder at any bound it is given. Resolving one means cells a hundred times finer in every
//! direction, which is a million times the work to model a joint that has no interesting field
//! inside it.
//!
//! A contact has a resistance and **no thickness**, so there is no floor: the same joint is the
//! same number at every cell size, and the cell size is chosen for the part rather than for its
//! thinnest feature.
//!
//! `a_layer_has_a_floor_at_one_cell_and_a_contact_has_none` measures both sides of that rather
//! than describing it.

use pantometry_core::{
    units::{Area, Energy, Length, Temperature, ThermalConductivity, Time},
    Domain, Exchange, Substance,
};
use pantometry_thermal::{Axis, Environment, Face, Solid3D};

/// The conductivity of the bar's material, taken from the substance rather than written down, so
/// the closed forms below cannot drift from what the block is made of.
fn k_bulk() -> f64 {
    grey()
        .thermal
        .expect("aluminium conducts")
        .conductivity
        .to_si()
}

/// Aluminium with radiation switched off.
///
/// The film's conductance carries `4εσT³`, so with an emissivity the balance below is the root of
/// a quartic rather than a division and the closed form would be a linearisation's. `steady_state`
/// makes the same substance for the same reason.
fn grey() -> Substance {
    let mut s = Substance::aluminium_6061();
    if let Some(t) = s.thermal.as_mut() {
        t.emissivity = 0.0;
    }
    s
}

/// A bar of `n` cubic cells along x, of side `dx`.
fn bar(n: usize, dx: f64) -> Solid3D {
    Solid3D::new(
        "bar",
        grey(),
        (n, 1, 1),
        Length::from_si(dx),
        Temperature::celsius(20.0),
    )
}

/// The chain of face resistances from the first cell centre to the last, in K/W.
///
/// Read out of `face_conductance`, which the crate documents as "the number the sweep actually
/// uses" — so this is the resistance the march sees and not a restatement of the formula.
fn chain(b: &Solid3D) -> f64 {
    let (nx, _, _) = b.counts();
    (1..nx)
        .map(|i| {
            1.0 / b
                .face_conductance((i - 1, 0, 0), (i, 0, 0))
                .expect("neighbours along x")
                .to_si()
        })
        .sum()
}

/// **An unstated contact is the block that has none, bit for bit.**
///
/// `joined` adds a branch to the innermost loop of `resolve`, which rebuilds every face
/// conductivity in the block. If `h = +∞` did not return the harmonic mean *exactly*, every closed
/// form the domain is already checked against would have been quietly reinterpreted — and by an
/// amount too small to fail any of them, which is the failure that matters here.
///
/// # The conductivity is 100.018 W/m/K and that is not arbitrary
///
/// Without the guard, an infinite contact computes `1/(1/k + 1/(dx·∞))` = `1/(1/k)`, and a double
/// reciprocal is **usually** the identity: about one double in six is not, and aluminium's 167,
/// borosilicate's 1.114 and the harmonic mean between them all happen to be in the five that are.
/// This test was written with those three and **passed with the guard deleted** — it was measuring
/// the arithmetic's luck rather than the code.
///
/// 100.018 W/m/K makes a face whose double reciprocal is one unit in the last place away, and the
/// test asserts that before it asserts anything else: a check that cannot see the defect is worth
/// stating as a prerequisite rather than assuming.
#[test]
fn a_contact_of_infinity_leaves_the_block_exactly_as_it_was() {
    // Chosen for its reciprocal, and held to it. A later edit that rounds this number puts the
    // test back to measuring luck, and this is what says so.
    const AWKWARD: f64 = 100.018;
    let face = 2.0 * AWKWARD * AWKWARD / (AWKWARD + AWKWARD);
    assert!(
        1.0 / (1.0 / face) != face,
        "{face:.17e} survives a double reciprocal, so this test cannot see a missing guard"
    );

    let build = || {
        let mut m = grey();
        if let Some(t) = m.thermal.as_mut() {
            t.conductivity = ThermalConductivity::w_per_m_k(AWKWARD);
        }
        let mut b = Solid3D::new(
            "bar",
            m,
            (9, 1, 1),
            Length::from_si(2e-3),
            Temperature::celsius(20.0),
        );
        // A second material as well, so the two-material face goes through the same branch.
        b.fill(Substance::borosilicate_crown(), |i, _, _| i >= 7);
        b
    };
    let plain = build();
    // Every interior face named, so this is not passing by never being consulted.
    let joined = build().joined(|_, _, _, _| Some(f64::INFINITY));

    for i in 1..9 {
        let a = plain
            .face_conductance((i - 1, 0, 0), (i, 0, 0))
            .expect("neighbours")
            .to_si();
        let b = joined
            .face_conductance((i - 1, 0, 0), (i, 0, 0))
            .expect("neighbours")
            .to_si();
        assert!(a - b == 0.0, "face {i}: {a:.17e} against {b:.17e}");
    }

    let (mut plain, mut joined) = (plain, joined);
    plain.deposit(0, 0, 0, Energy::from_si(5.0));
    joined.deposit(0, 0, 0, Energy::from_si(5.0));
    let dt = Time::from_si(plain.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..400 {
        plain
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
        joined
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    let (a, b) = (
        plain.peak_temperature().to_si(),
        joined.peak_temperature().to_si(),
    );
    println!("  peak {a:.17e} against {b:.17e}");
    assert!(a - b == 0.0, "an infinite contact moved the answer");
    // And a face nobody stated anything about reports no contact.
    assert_eq!(joined.contact_at(Axis::X, 4, 0, 0), None);
}

/// **A contact adds exactly `1/(A·h)` to the face it is on, and nothing to any other.**
///
/// The load-bearing equality. `A·(R_with − R_without)` is the areal resistance the joint
/// contributes, and it has a closed form with no discretisation in it at all: the two half cells
/// are unchanged, so the whole difference is the contact.
///
/// Checked across five decades of `h` — from a dry bolted joint at 500 W·m⁻²·K⁻¹ to a soldered one
/// at 5e6 — because a formula with `dx` on the wrong side of a division agrees at exactly one
/// value and this is how that is caught.
#[test]
fn a_contact_adds_the_areal_resistance_its_conductance_says() {
    let (n, dx) = (9usize, 2e-3);
    let area = dx * dx;
    for h in [5e2, 5e3, 5e4, 5e5, 5e6] {
        let plain = bar(n, dx);
        let joined = bar(n, dx).joined(|axis, i, _, _| (axis == Axis::X && i == 5).then_some(h));
        let added = area * (chain(&joined) - chain(&plain));
        let closed = 1.0 / h;
        println!("  h = {h:>9.0}: added {added:.6e}, closed form {closed:.6e}");
        assert!(
            (added / closed - 1.0).abs() < 1e-9,
            "h = {h}: a contact added {added:.6e} m²K/W where 1/h is {closed:.6e}"
        );
        // On that face and no other: the rest of the chain is untouched.
        for i in 1..n {
            if i == 5 {
                continue;
            }
            let a = plain
                .face_conductance((i - 1, 0, 0), (i, 0, 0))
                .expect("neighbours")
                .to_si();
            let b = joined
                .face_conductance((i - 1, 0, 0), (i, 0, 0))
                .expect("neighbours")
                .to_si();
            assert!(a - b == 0.0, "face {i} moved and no contact was put on it");
        }
        assert_eq!(joined.contact_at(Axis::X, 5, 0, 0), Some(h));
    }
}

/// **A layer's resistance has a floor at one cell, and a contact has none.**
///
/// The reason this capability exists, measured rather than asserted. Both blocks model the same
/// physical joint at the same physical place; only the way it is written differs.
///
/// - As a **contact**: the areal resistance is `1/h` at every cell size, exactly.
/// - As a **layer one cell thick**: `dx/k`, which is a different resistance at every cell size —
///   half as much when the cell halves, measured below.
///
/// # This is a floor and not a convergence failure, and the difference took a measurement
///
/// A caller who states a layer by its **bounds** rather than by a cell keeps its thickness when
/// the grid is refined, and `pantometry-world`'s resolution sweep does exactly that: it doubles a
/// region's `from` and `to` along with the counts. So a 1.5 mm layer stays 1.5 mm and the answer
/// does not move — which is right, and is not the point.
///
/// The point is that **1.5 mm was already the thinnest thing that grid could say.** A layer is at
/// least one cell, so a 1.5 mm cell cannot state 100 µm of solder at any bound; the file says
/// 1.5 mm because that is the floor, and carries fifteen times the resistance the joint has. This
/// test measures the floor by holding the layer at one cell and halving it — the shape a caller
/// gets from `fill(|_, _, k| k == 5)`, which is how a layer thinner than the grid is usually
/// written — and the contact beside it, which has no floor to measure.
#[test]
fn a_layer_has_a_floor_at_one_cell_and_a_contact_has_none() {
    let (n, dx, h) = (8usize, 2e-3, 5e3);
    // The conductivity that makes a single cell of width `dx` worth `1/h`: `dx/k = 1/h`.
    let k_layer = dx * h;

    let areal_contact = |cells: usize, side: f64, at: usize| {
        let plain = bar(cells, side);
        let joined =
            bar(cells, side).joined(|axis, i, _, _| (axis == Axis::X && i == at).then_some(h));
        side * side * (chain(&joined) - chain(&plain))
    };
    // The same joint as a cell of poor material, and its areal resistance measured the same way.
    let areal_layer = |cells: usize, side: f64, at: usize| {
        let mut thin = grey();
        if let Some(t) = thin.thermal.as_mut() {
            t.conductivity = ThermalConductivity::w_per_m_k(k_layer);
        }
        let plain = bar(cells, side);
        let mut layered = bar(cells, side);
        layered.fill(thin, |i, _, _| i == at);
        side * side * (chain(&layered) - chain(&plain))
    };

    // The coarse grid puts the joint on the face before cell 4; the fine grid halves every cell,
    // so the same physical plane is the face before cell 8 and the same physical layer is cell 8.
    let (coarse_c, fine_c) = (areal_contact(n, dx, 4), areal_contact(2 * n, dx / 2.0, 8));
    let (coarse_l, fine_l) = (areal_layer(n, dx, 4), areal_layer(2 * n, dx / 2.0, 8));
    println!(
        "  contact: {coarse_c:.6e} at {n} cells, {fine_c:.6e} at {} cells",
        2 * n
    );
    println!(
        "  layer:   {coarse_l:.6e} at {n} cells, {fine_l:.6e} at {} cells",
        2 * n
    );
    println!(
        "  the layer gave back {:.1}% of its resistance",
        100.0 * (1.0 - fine_l / coarse_l)
    );

    // The contact does not move. An equality, because refining the grid does not touch it at all.
    assert!(
        (fine_c / coarse_c - 1.0).abs() < 1e-9,
        "a contact moved on refinement: {coarse_c:.6e} then {fine_c:.6e}"
    );
    assert!(
        (coarse_c * h - 1.0).abs() < 1e-9,
        "and it is 1/h at both: {coarse_c:.6e} against {:.6e}",
        1.0 / h
    );
    // The layer halves, which is what `dx/k` does when `dx` halves. Not a tolerance on a
    // discretisation — a statement that the answer is proportional to the cell.
    assert!(
        (fine_l / coarse_l - 0.5).abs() < 0.02,
        "a one-cell layer went from {coarse_l:.6e} to {fine_l:.6e}, a factor of {:.3} rather than 2",
        coarse_l / fine_l
    );
}

/// The same bar, generating at one end and shedding at the other, with or without a joint.
fn driven(watts: f64, n: usize, dx: f64, h_film: f64, ambient: f64, joint: Option<f64>) -> Solid3D {
    let block = match joint {
        Some(h) => bar(n, dx).joined(|axis, i, _, _| (axis == Axis::X && i == n / 2).then_some(h)),
        None => bar(n, dx),
    };
    block.dissipating(watts, |i, _, _| i == 0).losing_from(
        Face::XMax,
        Environment {
            ambient: Temperature::celsius(ambient),
            // No radiation, so the balance is linear and the closed form is the problem's
            // rather than a linearisation's.
            convection_w_per_m2_k: h_film,
            area: Area::from_si(dx * dx),
        },
    )
}

/// **A bar with a joint drops what the series says**, marched to steady state.
///
/// The check above is on the operator; this is on the answer the sweep produces from it. At steady
/// state every watt generated at one end crosses every face and leaves through the film at the
/// other, so the rise from the first cell's centre to ambient is
///
/// ```text
/// P · ( (n − ½)·dx/(k·A)  +  1/(h·A)  +  1/(h_film·A) )
/// ```
///
/// with every term computed here from the constants.
///
/// **The half cell is the part that was measured wrong first.** `n − 1` cell widths is the path
/// from the source cell's centre to the last cell's centre, and the first version of this stopped
/// there — 222.96 K against a settled 226.69, low by 3.73 K. The film does not act at that centre:
/// `loss_conductance_at` puts it in series with the half cell of solid between the centre and the
/// surface, which is documented there and is what the missing 3.73 K is, to four figures.
#[test]
fn a_bar_with_a_joint_drops_what_the_series_says() {
    let (watts, n, dx, h) = (5.0, 12usize, 4e-3, 4e3);
    let (area, h_film, ambient) = (dx * dx, 5000.0, 20.0);

    let mut joined = driven(watts, n, dx, h_film, ambient, Some(h));
    settle(&mut joined);
    let rise = joined.temperature_at(0, 0, 0).to_si() - Temperature::celsius(ambient).to_si();
    let closed = watts
        * ((n as f64 - 0.5) * dx / (k_bulk() * area) + 1.0 / (h * area) + 1.0 / (h_film * area));
    println!("  rise {rise:.4} K, closed form {closed:.4} K");
    assert!(
        (rise / closed - 1.0).abs() < 2e-3,
        "a bar with a joint rose {rise:.4} K where the series says {closed:.4} K"
    );

    // And the joint is a real part of it: without the contact the same bar is cooler by exactly
    // the resistance the joint adds, `P/(hA)`.
    let mut plain = driven(watts, n, dx, h_film, ambient, None);
    settle(&mut plain);
    let without = plain.temperature_at(0, 0, 0).to_si() - Temperature::celsius(ambient).to_si();
    let jump = rise - without;
    let expected = watts / (h * area);
    println!("  the joint itself: {jump:.4} K, closed form {expected:.4} K");
    assert!(
        (jump / expected - 1.0).abs() < 5e-3,
        "the joint accounted for {jump:.4} K where P/(hA) is {expected:.4} K"
    );
}

/// **A contact of zero carries nothing**, which is a clearance.
///
/// The other end of the range, and it has to be exact rather than very small: `1/(dx·0)` is an
/// infinity and the guard is on the conductance rather than on the reciprocal. A block cut in two
/// this way has a half that cannot reach the other, and the far half stays where it started to the
/// last bit.
#[test]
fn a_contact_of_zero_is_a_clearance() {
    let (n, dx) = (8usize, 2e-3);
    let mut split = bar(n, dx)
        .joined(|axis, i, _, _| (axis == Axis::X && i == 4).then_some(0.0))
        .dissipating(10.0, |i, _, _| i == 0);
    assert!(
        split
            .face_conductance((3, 0, 0), (4, 0, 0))
            .expect("neighbours")
            .to_si()
            == 0.0,
        "a zero contact still conducts"
    );
    let cold_before = split.temperature_at(n - 1, 0, 0).to_si();
    let dt = Time::from_si(split.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..2000 {
        split
            .step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    let hot = split.temperature_at(0, 0, 0).to_si();
    let cold = split.temperature_at(n - 1, 0, 0).to_si();
    println!("  the fed half reached {hot:.2} K; the far half is at {cold:.6} K");
    assert!(
        hot > cold_before + 50.0,
        "the fed half did not heat: {hot:.2}"
    );
    assert!(
        cold - cold_before == 0.0,
        "heat crossed a clearance: {cold:.17e} from {cold_before:.17e}"
    );
}

/// March to steady state, checking that it arrived rather than marching a while and calling it
/// settled — the same driver `a_block_that_generates` uses, and for the reasons stated there.
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

/// **A negative contact conductance is refused, not quietly turned into an insulator.**
///
/// A face that carries a negative conductance would deliver heat up its own gradient, which is not
/// a poor joint but a mistake — a sign error in a caller, or a scene file with a minus in it. The
/// first version of this guard clamped such a value to zero, and an insulator is plausible enough
/// that nobody would have looked: a joint stated as `-5000` would have modelled a **clearance**
/// and the run would have completed and reported four figures.
///
/// It becomes `NaN`, which reaches `worst_rate` and refuses the sweep — the road a substance that
/// does not say what it conducts already takes. **The message is not that road's**, though, and
/// that took a second pass: a bad joint and a substance with no conductivity leave the same
/// unsayable limit, and reporting "substance has no diffusivity" sends a caller to their
/// materials, which is the one place the fault is not. `bad_contact` is what makes the refusal
/// name the joint, and the wording is pinned below.
#[test]
fn a_negative_contact_is_refused_rather_than_insulating() {
    let (n, dx) = (8usize, 2e-3);
    let mut wrong = bar(n, dx).joined(|axis, i, _, _| (axis == Axis::X && i == 4).then_some(-5e3));
    // `max_stable_dt` answers `+∞` for a block it cannot say anything about, not `NaN` — see its
    // own implementation. Infinity is the loud value here: a driver that halves it and steps is
    // refused on the next line rather than handed a plausible number.
    assert!(
        wrong.max_stable_dt(Time::ZERO).to_si().is_infinite(),
        "a negative conductance left a usable step: {}",
        wrong.max_stable_dt(Time::ZERO).to_si()
    );
    let refused = wrong.step(Time::ZERO, Time::from_si(1e-6), &mut Exchange::new());
    println!("  the sweep said: {refused:?}");
    let said = refused.expect_err("a negative conductance was stepped with");
    assert!(
        said.quantity.contains("contact conductance") && said.quantity.contains("negative"),
        "the refusal blamed something else: {}",
        said.quantity
    );

    // And a `NaN` takes the same road, which is what a scene file's missing number becomes.
    let mut absent =
        bar(n, dx).joined(|axis, i, _, _| (axis == Axis::X && i == 4).then_some(f64::NAN));
    let said = absent
        .step(Time::ZERO, Time::from_si(1e-6), &mut Exchange::new())
        .expect_err("a NaN conductance was stepped with");
    assert!(
        said.quantity.contains("not a number"),
        "the refusal blamed something else: {}",
        said.quantity
    );
}

/// **A mounting adds exactly `1/(h·A)` to the boundary path**, marched to steady state.
///
/// `joined` cannot reach the one face that matters most in a real design: the outer one, where a
/// part is bolted to its heatsink. `mounted_on` is that face's half, and the check is the same
/// series with one more term:
///
/// ```text
/// P · ( (n − ½)·dx/(k·A)  +  1/(h_mount·A)  +  1/(h_film·A) )
/// ```
///
/// Measured across three decades of mounting, because a term dropped or inverted agrees at one
/// value. The bar is grey, so the film has no radiative half and the closed form is the problem's.
#[test]
fn a_mounting_adds_the_resistance_its_conductance_says() {
    let (watts, n, dx) = (5.0, 10usize, 4e-3);
    let (area, h_film, ambient) = (dx * dx, 5000.0, 20.0);
    for mount in [1e3, 1e4, 1e5] {
        let mut block = bar(n, dx)
            .dissipating(watts, |i, _, _| i == 0)
            .losing_from(
                Face::XMax,
                Environment {
                    ambient: Temperature::celsius(ambient),
                    convection_w_per_m2_k: h_film,
                    area: Area::from_si(area),
                },
            )
            .mounted_on(Face::XMax, mount);
        assert_eq!(block.mounting_of(Face::XMax), Some(mount));
        settle(&mut block);
        let rise = block.temperature_at(0, 0, 0).to_si() - Temperature::celsius(ambient).to_si();
        let closed = watts
            * ((n as f64 - 0.5) * dx / (k_bulk() * area)
                + 1.0 / (mount * area)
                + 1.0 / (h_film * area));
        println!("  mount {mount:>8.0}: rise {rise:.4} K, closed form {closed:.4} K");
        assert!(
            (rise / closed - 1.0).abs() < 2e-3,
            "mounted at {mount}: rose {rise:.4} K where the series says {closed:.4} K"
        );
    }
}

/// **An unmounted face is the face that was there**, and a mounting of zero cannot lose heat.
///
/// The two ends. Without the guard the first is a `1/∞` that has to come out as the film exactly,
/// and the second is a face bolted on through a perfect insulator — which must shed nothing at
/// all rather than a very small amount, because a block that cannot lose heat has no steady state
/// and a block that loses a little has one a long way off.
#[test]
fn a_face_with_no_mounting_is_unchanged_and_a_mounting_of_zero_sheds_nothing() {
    let (watts, n, dx) = (5.0, 10usize, 4e-3);
    let (area, h_film, ambient) = (dx * dx, 5000.0, 20.0);
    let build = |mount: Option<f64>| {
        let b = bar(n, dx).dissipating(watts, |i, _, _| i == 0).losing_from(
            Face::XMax,
            Environment {
                ambient: Temperature::celsius(ambient),
                convection_w_per_m2_k: h_film,
                area: Area::from_si(area),
            },
        );
        match mount {
            Some(m) => b.mounted_on(Face::XMax, m),
            None => b,
        }
    };
    let (mut plain, mut open) = (build(None), build(Some(f64::INFINITY)));
    settle(&mut plain);
    settle(&mut open);
    let (a, b) = (
        plain.temperature_at(0, 0, 0).to_si(),
        open.temperature_at(0, 0, 0).to_si(),
    );
    println!("  no mounting {a:.17e}, an infinite one {b:.17e}");
    assert!(
        (a - b).abs() < 1e-9,
        "an infinite mounting moved the answer: {a:.17e} against {b:.17e}"
    );
    assert_eq!(plain.mounting_of(Face::XMax), None);

    // A perfect insulator: the block keeps every joule it makes, so its books say nothing left.
    let mut shut = build(Some(0.0));
    let dt = Time::from_si(shut.max_stable_dt(Time::ZERO).to_si() * 0.5);
    for _ in 0..4000 {
        shut.step(Time::ZERO, dt, &mut Exchange::new())
            .expect("stable");
    }
    println!(
        "  shut: {:.6} J lost after {:.1} s",
        shut.lost_energy().to_si(),
        4000.0 * dt.to_si()
    );
    assert!(
        shut.lost_energy().to_si() == 0.0,
        "a face mounted through a perfect insulator lost {:.6e} J",
        shut.lost_energy().to_si()
    );
}
