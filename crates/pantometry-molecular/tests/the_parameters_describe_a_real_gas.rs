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

/// **Every gas's own triple point, from its fitted parameters and the model's.**
#[test]
fn the_parameters_predict_each_gases_triple_point() {
    for s in CATALOGUE {
        let temperature = s.temperature(LJ_TRIPLE_T);
        let density = s.density(LJ_TRIPLE_DENSITY);
        let dt = (temperature - s.triple_point_k).abs() / s.triple_point_k;
        let dd = (density - s.triple_point_density).abs() / s.triple_point_density;
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
        // point and a real one across this whole column. A transcription error in `σ` or `ε` — a
        // digit, a unit, a gas — is tens of per cent or more, and every one of those fails here.
        assert!(
            dt < 0.05,
            "{}'s parameters put its triple point at {temperature:.2} K, {:.1}% from the measured \
             {:.2} K",
            s.name,
            100.0 * dt,
            s.triple_point_k
        );
        // The density is the sharper of the two, because it depends on `σ³`: a one per cent error
        // in the length shows here as three.
        assert!(
            dd < 0.05,
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
