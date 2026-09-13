//! The elastic network against things that are true of it whatever the structure is.
//!
//! Convention 1 asks for a closed form rather than a second implementation, and a network model
//! has an unusual number of them — because most of what it claims is a theorem about the Hessian
//! rather than a number that has to be computed to be known.
//!
//! | | why it is exact |
//! | --- | --- |
//! | Six null modes | `V` is a function of distances and a rigid motion changes none of them |
//! | A straight chain is `γ · 4sin²(jπ/2N)` | Its Hessian *is* the one-dimensional chain, times `γ`, with the two transverse directions carrying no stiffness at all |
//! | A ligand raises every eigenvalue | Weyl: adding springs adds a positive semi-definite matrix |
//! | ... and lowers every fluctuation | `H₂ ⪰ H₁` gives `H₂⁺ ⪯ H₁⁺` on the shared range |
//! | `(3n − r) k_BT` of energy | Equipartition, one `k_BT` per mode, and the amplitude was chosen to put it there |
//! | Doubling `γ` halves every fluctuation | `γ` multiplies `H`, so it divides `H⁺`, and cancels from every shape |
//!
//! # The structure these use
//!
//! An alpha helix: 2.3 Å radius, 1.5 Å rise and 100° of turn per residue. It is the real geometry
//! rather than a lattice, it makes a connected network at any sensible cutoff, and it is not
//! symmetric enough for the spectrum to be accidentally degenerate — a cube or a ring would be,
//! and a degenerate spectrum hides a wrong Hessian behind a subspace nothing can distinguish
//! inside.

use pantometry_core::{Domain, Exchange};
use pantometry_protein::modes::correlation;
use pantometry_protein::spectrum::Symmetric;
use pantometry_protein::{Modes, Network, Protein, Structure};
use pantometry_units::{Length, Mass, Qty, Temperature, Time};

/// Ångström, in metres.
const A: f64 = 1e-10;

/// One kcal/mol/Å², the stiffness the anisotropic network model is usually quoted with, in N/m.
fn gamma() -> Qty<0, 1, -2, 0, 0, 0, 0> {
    Qty::from_si(0.695)
}

fn cutoff(angstrom: f64) -> Length {
    Qty::from_si(angstrom * A)
}

/// `n` residues of an alpha helix.
fn helix(n: usize) -> Structure {
    Structure::from_positions(
        (0..n)
            .map(|i| {
                let turn = i as f64 * 100.0f64.to_radians();
                [
                    2.3 * A * turn.cos(),
                    2.3 * A * turn.sin(),
                    1.5 * A * i as f64,
                ]
            })
            .collect(),
    )
}

/// `‖Hv‖ / (‖H‖_F ‖v‖)`, so a "this is zero" claim is relative to the matrix it is zero in.
fn relative_annihilation(h: &Symmetric, v: &[f64]) -> f64 {
    let hv = h.multiply(v);
    let norm = |x: &[f64]| x.iter().map(|a| a * a).sum::<f64>().sqrt();
    norm(&hv) / (h.frobenius_squared().sqrt() * norm(v))
}

/// **The six rigid-body motions are in the null space, exactly.**
///
/// Three translations and three rotations, written down from the geometry and multiplied through
/// the Hessian. This is the check the assembly has to pass, and it is sensitive to every way of
/// getting the assembly wrong: a diagonal block added twice, a sign flipped on the off-diagonal
/// block, the outer product transposed, a pair counted once instead of contributing to both
/// nodes. None of those leaves a matrix that annihilates a translation.
#[test]
fn a_hessian_annihilates_every_rigid_motion() {
    let s = helix(30);
    let network = Network::new(&s, cutoff(10.0), gamma());
    let h = network.hessian();
    let n = network.nodes();

    let mut worst: f64 = 0.0;
    for axis in 0..3 {
        // A translation: the same direction at every node.
        let mut v = vec![0.0; 3 * n];
        for i in 0..n {
            v[3 * i + axis] = 1.0;
        }
        worst = worst.max(relative_annihilation(&h, &v));

        // A rotation about the same axis through the origin: `ω × r`.
        let mut w = vec![0.0; 3 * n];
        for i in 0..n {
            let r = network.at(i);
            let (a, b, c) = (axis, (axis + 1) % 3, (axis + 2) % 3);
            w[3 * i + a] = 0.0;
            w[3 * i + b] = -r[c];
            w[3 * i + c] = r[b];
        }
        worst = worst.max(relative_annihilation(&h, &w));
    }
    assert!(
        worst < 1e-14,
        "a rigid motion left {worst:e} of the Hessian behind, so the assembly is wrong"
    );

    // And the spectrum agrees about how many there are, which is the same fact counted a second
    // way — through the eigensolver rather than through the geometry.
    let modes = Modes::of(&network);
    assert_eq!(modes.rigid_body_modes(), 6, "a helix is not collinear");
    assert!(
        modes.separation() > 1e6,
        "the six zeros are only {:e} below the first real mode, which is not a gap",
        modes.separation()
    );
    println!(
        "  helix of 30: rigid motions leave {worst:.3e}, separation {:.3e}",
        modes.separation()
    );
}

/// **A straight chain's network *is* the one-dimensional chain, times `γ`.**
///
/// Springs along a line have `ê` along the line, so `ê êᵀ` has one nonzero entry and the Hessian
/// is block diagonal: the longitudinal block is `γ` times the path graph's Laplacian, whose
/// spectrum `4 sin²(jπ/2N)` `spectrum.rs` already checks against, and the two transverse blocks
/// are **zero**.
///
/// That last part is the physics, not an artefact: a line of beads on springs has no resistance
/// to being bent, so a collinear structure has `2N + 1` zero modes and not six. The five rigid
/// ones are among them. A model that reported six here would be reporting a stiffness it does
/// not have.
#[test]
fn a_straight_chain_is_the_one_dimensional_chain_and_bends_for_free() {
    let n = 12;
    let s = Structure::from_positions((0..n).map(|i| [0.0, 0.0, 3.8 * A * i as f64]).collect());
    // Just past one spacing, so only neighbours are joined.
    let network = Network::new(&s, cutoff(4.0), gamma());
    assert_eq!(
        network.springs(),
        n - 1,
        "more than the neighbours are joined"
    );

    let spectrum = network.hessian().eigen();
    let g = gamma().to_si();
    let scale = network.hessian().frobenius_squared().sqrt();

    // `2N` transverse zeros and the one longitudinal zero, then the chain's own spectrum.
    let modes = Modes::of(&network);
    assert_eq!(
        modes.rigid_body_modes(),
        2 * n + 1,
        "a line of beads resists bending, which it must not"
    );
    let mut worst: f64 = 0.0;
    for j in 1..n {
        let exact = g
            * 4.0
            * (j as f64 * std::f64::consts::PI / (2.0 * n as f64))
                .sin()
                .powi(2);
        worst = worst.max((spectrum.values()[2 * n + j] - exact).abs());
    }
    let floor = (spectrum.sweeps() * 3 * n) as f64 * f64::EPSILON * scale;
    assert!(
        worst <= floor,
        "the longitudinal spectrum is off by {worst:e} from 4 gamma sin^2(j pi / 2N), over {floor:e}"
    );
    println!(
        "  straight chain: {} zeros, longitudinal worst {worst:.3e}",
        modes.rigid_body_modes()
    );
}

/// **A network in two pieces has twelve zeros, and says so before anything downstream divides.**
///
/// Two clusters further apart than the cutoff are two structures, and each brings its own six.
/// The failure this guards is silent: the twelfth zero is not an error, the thirteenth eigenvalue
/// is tiny rather than absent, and the "softest collective mode" is one cluster leaving.
///
/// The counts run everywhere. The refusal is checked only where a panic unwinds, for the reason
/// written beside it.
#[test]
fn a_network_in_pieces_is_caught_before_it_is_used() {
    let mut points: Vec<[f64; 3]> = helix(10).residues().iter().map(|r| r.at).collect();
    points.extend(
        helix(10)
            .residues()
            .iter()
            .map(|r| [r.at[0] + 100.0 * A, r.at[1], r.at[2]]),
    );
    let apart = Network::new(&Structure::from_positions(points), cutoff(10.0), gamma());
    assert_eq!(apart.components(), 2);
    assert_eq!(Modes::of(&apart).rigid_body_modes(), 12);

    // **Only where a panic unwinds.** `catch_unwind` returns `Err` by catching the unwind, and
    // `wasm32` panics by aborting: this compiled on the local gate, which runs the wasm tests
    // with `--no-run`, and trapped in CI with `wasm trap: unreachable` and exit 134 -- the
    // assertion firing exactly as designed and taking the process with it.
    //
    // Everything above runs on every target. What is cfg'd out is the refusal, and the reason is
    // a property of the target's panic strategy rather than of the check.
    #[cfg(panic = "unwind")]
    {
        let refused = std::panic::catch_unwind(|| {
            Protein::new(
                "split",
                &apart,
                Qty::from_si(300.0),
                Qty::from_si(1.8e-25),
                1,
            )
        });
        let message = refused.expect_err("a protein was built out of two structures");
        let said = message
            .downcast_ref::<String>()
            .map(String::as_str)
            .unwrap_or("");
        assert!(
            said.contains("falls into 2 pieces"),
            "it refused for some other reason: {said:?}"
        );
    }
    println!("  two clusters: 2 components, 12 zeros, and Protein::new refuses");
}

/// **A ligand lowers every fluctuation — and does not raise every eigenvalue.**
///
/// The first half is a theorem, and it is a theorem about a matrix that has to be asked for.
/// The effective Hessian on the protein's coordinates, with the ligand integrated out, is the
/// bound potential minimised over where the ligand can go; the coupling terms are squares, so
/// that minimum is at least the bare potential, and the covariance `k_BT H⁺` therefore falls
/// everywhere. A residue that got **more** mobile on binding would be violating it.
///
/// Taking it from the complex's own Hessian instead — which is the obvious thing to do and what
/// this crate did first — breaks it: a pseudo-inverse removes the null space it is handed, the
/// complex's rigid motions are not the protein's, and residue 4 of this helix comes out
/// **1.003** times as mobile after binding. Through [`Modes::of_the_protein`] the same residue
/// comes out at **0.828**. The counter-example is kept below.
///
/// The second half is what this test was first written to assert, and it is false. A ligand adds
/// *nodes*, so the bound Hessian is bigger than the bare one and "every eigenvalue rises" is not
/// Weyl's inequality — Weyl compares two matrices on one space. The measured counter-example is
/// kept here rather than deleted, because the claim is a plausible one that a reviewer will reach
/// for again.
#[test]
fn binding_a_ligand_stiffens_everything_and_loosens_nothing() {
    let s = helix(24);
    let bare = Network::new(&s, cutoff(10.0), gamma());
    // Three atoms tucked against the middle of the helix.
    let middle = s.residues()[12].at;
    let bound = Network::new(&s, cutoff(10.0), gamma()).with_ligand(vec![
        [middle[0] + 4.0 * A, middle[1], middle[2]],
        [middle[0] + 4.0 * A, middle[1] + 3.0 * A, middle[2]],
        [middle[0] + 4.0 * A, middle[1], middle[2] + 3.0 * A],
    ]);
    assert!(
        bound.springs() > bare.springs(),
        "the ligand touched nothing"
    );

    let (a, b) = (Modes::of(&bare), Modes::of(&bound));
    assert_eq!(a.rigid_body_modes(), 6);
    assert_eq!(b.rigid_body_modes(), 6);
    // The protein's own coordinates, which is what a residue's fluctuation is a property of.
    let (pa, pb) = (Modes::of_the_protein(&bare), Modes::of_the_protein(&bound));
    assert_eq!(pa, a, "with no ligand the two Hessians are the same matrix");
    assert_eq!(pb.nodes(), s.len());
    assert_eq!(pb.rigid_body_modes(), 6);

    // The spectrum, which does *not* simply rise: find the counter-example rather than assert
    // its absence, so that a change making the claim true again fails here and gets read.
    let fell = (0..a.len().min(b.len()))
        .find(|&k| b.eigenvalue(k).to_si() < a.eigenvalue(k).to_si() * (1.0 - 1e-9));
    let fell = fell.expect(
        "no mode fell, so either the ligand is not attached or the bound Hessian is the bare one \
         with extra rows -- either way the comparison below is measuring nothing",
    );
    println!(
        "  spectrum: mode {fell} falls from {:.4e} to {:.4e} N/m, because the space got bigger",
        a.eigenvalue(fell).to_si(),
        b.eigenvalue(fell).to_si()
    );

    let t = Qty::from_si(300.0);
    let (before, after) = (pa.fluctuations(t), pb.fluctuations(t));
    let mut worst_rise: f64 = 0.0;
    let mut biggest_drop = (0usize, 0.0f64);
    for i in 0..s.len() {
        let ratio = after[i] / before[i];
        worst_rise = worst_rise.max(ratio);
        if 1.0 - ratio > biggest_drop.1 {
            biggest_drop = (i, 1.0 - ratio);
        }
    }
    assert!(
        worst_rise <= 1.0 + 1e-12,
        "a residue got {worst_rise} times as mobile after binding, which the effective Hessian \
         forbids"
    );

    // And the counter-example, so that the wrong matrix stays visibly wrong: through the
    // complex's own Hessian a residue does come out more mobile, which is how the two were
    // told apart in the first place.
    let naive = b.fluctuations(t);
    let bare_full = a.fluctuations(t);
    let risen = (0..s.len())
        .map(|i| naive[i] / bare_full[i])
        .fold(0.0f64, f64::max);
    assert!(
        risen > 1.0,
        "the complex's Hessian no longer disagrees with the protein's, so the distinction \
         `Modes::of_the_protein` exists for is untested"
    );
    println!("  the complex's own Hessian would have said {risen:.6} for the worst residue");
    println!(
        "  ligand: {} springs added, every mode stiffer, residue {} lost {:.1}% of its motion",
        bound.springs() - bare.springs(),
        biggest_drop.0,
        100.0 * biggest_drop.1
    );
    assert!(
        biggest_drop.1 > 0.01,
        "the largest effect anywhere was {:.3}%, so this ligand is not bound to anything",
        100.0 * biggest_drop.1
    );
}

/// **The stiffness sets the scale and changes no shape.**
///
/// `γ` multiplies the Hessian, so it divides every fluctuation by exactly the same factor. That
/// is why a correlation against measured B-factors is a prediction and not a fit: there is no
/// value of `γ` that improves it.
#[test]
fn the_spring_constant_divides_every_fluctuation_and_moves_nothing_else() {
    let s = helix(20);
    let t = Qty::from_si(300.0);
    let soft = Modes::of(&Network::new(&s, cutoff(10.0), gamma()));
    let stiff = Modes::of(&Network::new(&s, cutoff(10.0), Qty::from_si(2.0 * 0.695)));

    let (a, b) = (soft.fluctuations(t), stiff.fluctuations(t));
    let mut worst: f64 = 0.0;
    for i in 0..s.len() {
        worst = worst.max((a[i] / b[i] - 2.0).abs());
    }
    assert!(worst < 1e-10, "the ratio strayed from two by {worst:e}");
    assert_eq!(
        correlation(&a, &b).map(|c| (c * 1e12).round()),
        Some(1e12),
        "doubling the stiffness reordered the residues"
    );
    println!("  stiffness: ratio held to {worst:.3e}, correlation exactly one");
}

/// **The structure moving in time carries `(3n − 6) k_BT`, and the motion averages to the
/// fluctuations.**
///
/// Two claims, and the second is the one that ties the time evolution to the theory. The
/// amplitudes were chosen from equipartition, so the *static* prediction
/// [`Modes::fluctuations`] and the *time average* of the moving structure have to agree — if
/// they do not, either the amplitude or the fluctuation formula is wrong and neither would
/// otherwise say so.
///
/// # The sampling is derived, not picked
///
/// Two separate requirements, and the first version of this test met only one of them:
///
/// - **The window** has to hold many periods of the *softest* mode, because the finite-window
///   error in `⟨cos²⟩` falls as `1/M` with the number of periods `M`, and the softest mode
///   carries the most weight (`1/λ`).
/// - **The step** has to resolve the *fastest* one. A discrete average of `cos²` sampled near its
///   own period aliases, and the error does not care that there are many samples.
///
/// Four thousand samples over thirty-seven soft periods satisfied the first and failed the
/// second — the fastest mode was getting about one sample per period — and the time average came
/// out **9.0%** off.
///
/// # And it is an ensemble average, not one run's
///
/// Fixing the sampling left **4.7%**, and that part is real. The mean square displacement of one
/// realisation carries cross terms `⟨cos(ω_k t + φ_k) cos(ω_l t + φ_l)⟩` between different modes,
/// and those decay as `1/((ω_k − ω_l)T)` — which for a near-degenerate pair is no decay at all.
/// They vanish on averaging over the phases, not over time, because that is what the prediction
/// is an average over. Measured on a ten-residue helix, worst residue: **2.5%** at one seed,
/// **2.3%** at four, **0.58%** at sixteen, **0.82%** at sixty-four — a `1/√S` with a `1/√S`
/// worth of noise on it.
#[test]
fn the_motion_averages_to_the_fluctuations_it_was_built_from() {
    let s = helix(10);
    let network = Network::new(&s, cutoff(10.0), gamma());
    let t = Qty::from_si(300.0);
    let mass: Mass = Qty::from_si(1.8e-25);
    let modes = Modes::of(&network);
    let kt = 1.380_649e-23 * 300.0;
    let expected = modes.len() as f64 * kt;
    let held = Protein::new("helix", &network, t, mass, 7)
        .ledger()
        .get("energy")
        .expect("an energy ledger");
    assert!(
        (held / expected - 1.0).abs() < 1e-12,
        "the ledger holds {held:e} and equipartition says {expected:e}"
    );
    assert_eq!(modes.len(), 3 * s.len() - 6);

    // Sample the motion over many periods of the softest mode and compare the mean square
    // displacement per node against the static prediction.
    let predicted = modes.fluctuations(t);
    let rate = |k: usize| (modes.eigenvalue(k).to_si() / mass.to_si()).sqrt();
    let (slowest, fastest) = (rate(0), rate(modes.len() - 1));
    // Fifty periods of the softest mode, eight samples per period of the fastest.
    let windows = 50.0;
    let step_s = std::f64::consts::TAU / (8.0 * fastest);
    let samples = (windows * std::f64::consts::TAU / slowest / step_s).ceil() as usize;
    let step: Time = Qty::from_si(step_s);
    let seeds = 64u64;
    let mut mean = vec![0.0; s.len()];
    let mut bus = Exchange::new();
    let mut protein = Protein::new("helix", &network, t, mass, 0);
    for seed in 0..seeds {
        protein = Protein::new("helix", &network, t, mass, seed);
        let mut clock = 0.0;
        for _ in 0..samples {
            protein
                .step(Qty::from_si(clock), step, &mut bus)
                .expect("the protein steps");
            clock += step.to_si();
            for (i, slot) in mean.iter_mut().enumerate() {
                *slot += protein.displacement(i).to_si().powi(2) / (samples as f64 * seeds as f64);
            }
        }
    }
    let mut worst: f64 = 0.0;
    for i in 0..s.len() {
        worst = worst.max((mean[i] / predicted[i] - 1.0).abs());
    }
    // `1/M` with `M` the periods of the softest mode: the window's own truncation of a cosine.
    let floor = 1.0 / windows;
    assert!(
        worst < floor,
        "the time average is off the static prediction by {:.2}% at its worst, over {:.2}%",
        100.0 * worst,
        100.0 * floor
    );
    // The energy did not move while all that happened.
    assert_eq!(protein.ledger().get("energy"), Some(expected));
    println!(
        "  ensemble average against equipartition: worst {:.2}% of {:.2}%, {seeds} seeds x \
         {samples} samples, frequency ratio {:.1}",
        100.0 * worst,
        100.0 * floor,
        fastest / slowest
    );
}

/// **The same seed gives the same run, and a different one gives a different member of the same
/// ensemble.**
#[test]
fn the_phases_are_reproducible_and_the_statistics_are_not_the_phases() {
    let s = helix(14);
    let network = Network::new(&s, cutoff(10.0), gamma());
    let (t, m): (Temperature, Mass) = (Qty::from_si(300.0), Qty::from_si(1.8e-25));
    let make = |seed| Protein::new("helix", &network, t, m, seed);
    let (a, b, c) = (make(1), make(1), make(2));
    let at = |p: &Protein| {
        (0..s.len())
            .map(|i| p.displacement(i).to_si())
            .collect::<Vec<_>>()
    };

    assert_eq!(at(&a), at(&b), "the same seed moved differently");
    assert_ne!(at(&a), at(&c), "two seeds gave the same phases");
    // Different phases, same ensemble: the equilibrium spread does not depend on the draw.
    assert_eq!(a.spread(), c.spread());
    assert_eq!(a.ledger().get("energy"), c.ledger().get("energy"));
    println!("  seeds 1 and 2: different motion, identical spread and energy");
}
