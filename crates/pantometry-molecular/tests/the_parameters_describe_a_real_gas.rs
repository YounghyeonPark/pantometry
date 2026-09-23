//! **`σ` and `ε` are fits, and this is what makes them more than that.**
//!
//! They come from a gas's second virial coefficient and its viscosity. The Lennard-Jones model's
//! own triple point is a simulation result for the potential — `T* = 0.694`, `ρ* = 0.84` — and the
//! two together predict where each gas freezes and how dense its liquid is there. A triple point
//! is measured by putting the gas in a cell and watching it, which has nothing to do with either
//! fit.
//!
//! So a substance that was invented, or transcribed wrongly, or quoted in the wrong units, fails
//! here. That is the whole reason the catalogue carries the measured triple point beside the
//! fitted parameters rather than only the parameters.
//!
//! # And where the model stops, printed rather than asserted tightly
//!
//! The same arithmetic gives the critical temperature, and there the potential is **4 to 6 per
//! cent** out on every gas. That is a sphere with a well standing in for an atom in the one place
//! where fluctuations on every scale decide the answer. The band below is where it has always
//! been; tightening it would be asserting the model is better than it is, and loosening it would
//! stop it saying anything.

use pantometry_molecular::substance::{
    CATALOGUE, LJ_CRITICAL_T, LJ_TRIPLE_DENSITY, LJ_TRIPLE_T, NAMES, UNSHIPPED,
};
use pantometry_molecular::Substance;

/// The band the model's own error is allowed to occupy, as a fraction.
///
/// Shared by the test that holds the catalogue inside it and by the one that measures what it lets
/// through, so tightening it cannot quietly leave the second one asserting nothing.
const MODEL_BAND: f64 = 0.05;

/// The two residuals of a substance's predicted triple point, **signed**, as fractions.
///
/// Signed because the sign is half of what they say. Argon's density sits *below* the measurement
/// and krypton's *above*, and that is the measurement which says no single correction to
/// `LJ_TRIPLE_DENSITY` exists: one constant cannot move both towards zero. The best single value
/// across the two, 0.8325, takes argon from 0.39% to 1.27% to take krypton from 2.21% to 1.31% --
/// a redistribution rather than a correction, and the reason this file has never had one.
fn residuals(s: &Substance) -> (f64, f64) {
    (
        (s.temperature(LJ_TRIPLE_T) - s.triple_point_k) / s.triple_point_k,
        (s.density(LJ_TRIPLE_DENSITY) - s.triple_point_density) / s.triple_point_density,
    )
}

/// **Every gas's own triple point, from its fitted parameters and the model's.**
#[test]
fn the_parameters_predict_each_gases_triple_point() {
    for s in CATALOGUE {
        let temperature = s.temperature(LJ_TRIPLE_T);
        let density = s.density(LJ_TRIPLE_DENSITY);
        let (signed_t, signed_d) = residuals(&s);
        let (dt, dd) = (signed_t.abs(), signed_d.abs());
        println!(
            "  {:<8} T {:>7.2} K against {:>7.2} ({:>5.2}%)   rho {:>7.0} against {:>7.0} ({:>5.2}%)",
            s.name,
            temperature,
            s.triple_point_k,
            100.0 * dt,
            density,
            s.triple_point_density,
            100.0 * dd
        );
        // **Five per cent**, which is the size of the disagreement between the model's triple
        // point and a real one across this whole column. It catches the gross transcription
        // errors and **not** a slipped digit: this comment used to say that "a digit, a unit, a
        // gas" was "tens of per cent or more, and every one of those fails here", and two of the
        // three are. A digit is 4.24%, measured in
        // `a_slipped_digit_passes_the_band_that_holds_the_model`, and the pinned residuals in
        // `each_gases_residual_is_where_it_was_left` are what catches that instead.
        assert!(
            dt < MODEL_BAND,
            "{}'s parameters put its triple point at {temperature:.2} K, {:.1}% from the measured \
             {:.2} K",
            s.name,
            100.0 * dt,
            s.triple_point_k
        );
        // The density is the sharper of the two, because it depends on `σ³`: a one per cent error
        // in the length shows here as three.
        assert!(
            dd < MODEL_BAND,
            "{}'s parameters put its triple-point liquid at {density:.0} kg/m³, {:.1}% from the \
             measured {:.0}",
            s.name,
            100.0 * dd,
            s.triple_point_density
        );
    }
}

/// **A one per cent error in `σ` is a three per cent error in the density**, and that is why the
/// density check is the one that catches a mistyped length.
///
/// Not a claim about the gases: a claim about the arithmetic, so the test above is known to be
/// sensitive to the number it is most likely to be wrong about.
#[test]
fn the_density_check_is_three_times_as_sharp_as_the_length_it_depends_on() {
    let argon = Substance::named("argon").expect("argon is in the catalogue");
    let mut stretched = argon;
    stretched.sigma *= 1.01;
    let (was, now) = (
        argon.density(LJ_TRIPLE_DENSITY),
        stretched.density(LJ_TRIPLE_DENSITY),
    );
    let moved = (was - now) / was;
    println!("  1% on sigma moves the density {:.2}%", 100.0 * moved);
    assert!(
        (0.029..0.030).contains(&moved),
        "a 1% length should move a density by 1 - 1/1.01³ = 2.94%, and moved it {:.3}%",
        100.0 * moved
    );
}

/// **The critical point is where this potential is worst, and the number is written down.**
#[test]
fn the_model_is_a_few_per_cent_out_at_the_critical_point() {
    let mut worst: f64 = 0.0;
    for s in CATALOGUE {
        let predicted = s.temperature(LJ_CRITICAL_T);
        let out = (predicted - s.critical_k).abs() / s.critical_k;
        println!(
            "  {:<8} Tc {:>7.2} K against {:>7.2} ({:>5.2}%)",
            s.name,
            predicted,
            s.critical_k,
            100.0 * out
        );
        worst = worst.max(out);
    }
    // **A band, not a bound.** Below 3% would be claiming a sphere with a well does better near a
    // critical point than it does; above 8% would stop the number saying anything at all. The
    // measured spread is 4.1% to 5.5%.
    assert!(
        (0.03..0.08).contains(&worst),
        "the worst critical-temperature error is {:.1}%, which is outside the band this model has \
         always been in — either the parameters moved or the model did",
        100.0 * worst
    );
}

/// **The two gases that did not pass, kept so the boundary is a measurement.**
///
/// They are out for opposite reasons and the de Boer parameter tells them apart. **Neon** is the
/// most quantum of the four at `Λ* = 0.594`, three times argon's, and no choice of `σ` and `ε`
/// fixes an atom whose position is spread over a quarter of its own well. **Xenon** is the *least*
/// quantum at 0.063 — so its 10.8% cannot be the potential, and it is not: it is a pair written
/// down here that could not be sourced.
///
/// Substituting a pair that passes is the one thing this whole file exists to stop, so xenon stays
/// out until a sourced one arrives.
#[test]
fn a_gas_this_model_cannot_describe() {
    let miss = |s: &Substance| {
        (s.density(LJ_TRIPLE_DENSITY) - s.triple_point_density).abs() / s.triple_point_density
    };
    let worst = CATALOGUE.iter().map(miss).fold(0.0f64, f64::max);
    for s in UNSHIPPED {
        println!(
            "  {:<8} out by {:>6.2}% on density, de Boer {:.3}",
            s.name,
            100.0 * miss(&s),
            s.de_boer()
        );
        assert!(
            miss(&s) > 2.0 * worst,
            "{} is out by {:.2}% and the catalogue's worst by {:.2}%; if those have converged then \
             either the model improved or a shipped entry has drifted",
            s.name,
            100.0 * miss(&s),
            100.0 * worst
        );
    }

    // **Which of the two reasons it is, measured.** Quantum spread against the well's own width:
    // everything shipped is below a quarter, neon is more than half, and xenon is the lowest of
    // all four — which is what says its miss is arithmetic and not physics.
    for s in CATALOGUE {
        println!("  {:<8} de Boer {:.3}", s.name, s.de_boer());
        assert!(
            s.de_boer() < 0.25,
            "{} has a de Boer parameter of {:.3}, which is where classical Lennard-Jones stops",
            s.name,
            s.de_boer()
        );
    }
    let neon = UNSHIPPED[0];
    let xenon = UNSHIPPED[1];
    assert!(
        neon.de_boer() > 0.5,
        "neon's de Boer parameter is {:.3}; the argument for leaving it out was that it is large",
        neon.de_boer()
    );
    assert!(
        xenon.de_boer() < CATALOGUE.iter().map(|s| s.de_boer()).fold(1.0f64, f64::min),
        "xenon is supposed to be the most classical of the four, which is why its miss is the \
         parameters' and not the potential's"
    );
}

/// **The names are the catalogue's, in order, and a name that is not there is refused.**
#[test]
fn a_substance_is_found_by_name_or_not_at_all() {
    let mut sorted: Vec<&str> = CATALOGUE.iter().map(|s| s.name).collect();
    sorted.sort_unstable();
    assert_eq!(
        sorted, NAMES,
        "`NAMES` is what a scene's error message lists"
    );
    for name in NAMES {
        assert!(
            Substance::named(name).is_some(),
            "{name} is in the catalogue"
        );
    }
    // Not a default, and not the nearest match: a misspelling that silently became argon would be
    // a run of the wrong gas reported in the right units.
    assert!(Substance::named("Argon").is_none(), "the names are exact");
    assert!(
        Substance::named("nitrogen").is_none(),
        "and this is not a sphere"
    );
    for s in UNSHIPPED {
        assert!(
            Substance::named(s.name).is_none(),
            "{} is measured in this file and is not offered to a scene",
            s.name
        );
    }
}

/// **A slipped digit in `σ` passes the five per cent band**, which is why the residuals are pinned.
///
/// `3.405` written as `3.450`, the last two digits transposed, is the likeliest transcription
/// error a table of four-figure lengths can suffer, and it moves argon's triple-point density by
/// **4.24%** against a bound of 5%. The two errors the bound does catch are measured beside it, so
/// the claim that replaced "a digit, a unit, a gas is tens of per cent or more" is itself a
/// measurement rather than a second sentence.
#[test]
fn a_slipped_digit_passes_the_band_that_holds_the_model() {
    let argon = Substance::named("argon").expect("argon is in the catalogue");
    let mut slipped = argon;
    slipped.sigma = 3.450e-10;
    let (_, d) = residuals(&slipped);
    println!(
        "  sigma 3.405 -> 3.450 moves the density residual to {:+.2}%",
        100.0 * d
    );
    assert!(
        d.abs() < MODEL_BAND,
        "a transposed digit is {:.2}% and the band is {:.0}%. If this now fails the band has been \
         tightened past a transposition, and the pinned residuals may no longer be the only thing \
         catching one -- check that before deleting anything",
        100.0 * d.abs(),
        100.0 * MODEL_BAND
    );

    // The two the old comment was right about, so what replaced it is held to the same standard.
    let mut wrong_gas = argon;
    wrong_gas.sigma = 4.10e-10;
    let mut wrong_unit = argon;
    wrong_unit.sigma = 3.405e-9;
    for (what, s) in [
        ("xenon's length on argon", wrong_gas),
        ("nanometres for angstroms", wrong_unit),
    ] {
        let (_, d) = residuals(&s);
        println!("  {what} moves it to {:+.1}%", 100.0 * d);
        assert!(
            d.abs() > MODEL_BAND,
            "{what} is {:.1}% and is the size of error this bound exists to catch",
            100.0 * d.abs()
        );
    }
}

/// **Each gas's residual is where it was left**, which is what catches the digit the band cannot.
///
/// Not a closed form, and not a claim that the parameters are right --
/// `the_parameters_predict_each_gases_triple_point` is that. This is a pin, of the kind the
/// determinism digest is: the residual is a deterministic function of six `f64` constants in
/// `substance.rs`, nothing samples and nothing consults anything, so **any** edit to `σ`, `ε`, the
/// mass or either measured value moves it and has to be re-recorded by somebody who looked at the
/// new number.
///
/// # Where the tolerance comes from
///
/// `1e-6` on a fraction, and it traces to the precision the parameters are quoted at. The smallest
/// change any of them can suffer is one digit in its last place: argon's `σ` is `3.405` to
/// `0.001 Å`, which is `8.8e-4` in the density residual, and its `ε/k` is `119.8` to `0.1`, which
/// is `8.3e-4` in the temperature one. The bound is **eight hundred times tighter than the
/// smallest transcription error this file can hold**, and far looser than the arithmetic's own
/// reproducibility, which is exact.
#[test]
fn each_gases_residual_is_where_it_was_left() {
    // Printed by the code and copied back, not computed alongside it -- a second implementation of
    // the same arithmetic would agree with a wrong one.
    const RECORDED: [(&str, f64, f64); 2] = [
        ("argon", -0.007_979_955, -0.003_901_637),
        ("krypton", -0.016_963_206, 0.022_141_977),
    ];
    assert_eq!(
        RECORDED.len(),
        CATALOGUE.len(),
        "a gas was added to the catalogue without a pinned residual, so nothing here would have \
         noticed its parameters being wrong by less than the model's own error"
    );
    for (name, want_t, want_d) in RECORDED {
        let s = Substance::named(name).unwrap_or_else(|| panic!("{name} is in the catalogue"));
        let (t, d) = residuals(&s);
        println!(
            "  {name:<8} T {:+.5}% (pinned {:+.5}%)   rho {:+.5}% (pinned {:+.5}%)",
            100.0 * t,
            100.0 * want_t,
            100.0 * d,
            100.0 * want_d
        );
        assert!(
            (t - want_t).abs() < 1e-6,
            "{name}'s triple-point temperature residual is {:+.6}% and was pinned at {:+.6}%. A \
             parameter moved: the new number is to be looked at, not re-pinned",
            100.0 * t,
            100.0 * want_t
        );
        assert!(
            (d - want_d).abs() < 1e-6,
            "{name}'s triple-point density residual is {:+.6}% and was pinned at {:+.6}%. Either \
             sigma moved -- it enters cubed, so a slipped length shows up here three times over \
             -- or the measurement it is compared against did",
            100.0 * d,
            100.0 * want_d
        );
    }
}
