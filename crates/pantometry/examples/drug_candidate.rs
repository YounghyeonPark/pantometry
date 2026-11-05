//! A molecule judged as a drug candidate, on what this repository has measured about it.
//!
//! ```text
//! cargo run --release --example drug_candidate
//! ```
//!
//! [`pantometry::candidate::assess`] reads the molecule and takes everything else as evidence:
//! log P, the binding free energy and the fraction unbound in plasma, each with its standard
//! error. Here the evidence is this repository's own free-energy calculations, each of which ran
//! for a day or more and is recorded in the CHANGELOG:
//!
//! - **benzene**: log P 1.92 ± 0.23 from the direct legs over two seeds, against experiment's 2.13.
//!   It binds T4 lysozyme L99A at −8.975 ± 0.279 kcal/mol with UFF throughout, and at −4.61 ± 0.29
//!   with GAFF's van der Waals, ff14SB's protein and AM1-BCC charges, against experiment's
//!   −5.19 ± 0.16;
//! - **p-xylene**: log P 3.06 ± 0.22, against 3.15. Its binding was not computed;
//! - **aspirin**: nothing measured. What can be read off the molecule is all it gets.
//!
//! No fraction unbound was measured for any of them. The tissue partitions are therefore printed
//! at `f_u = 1`, the largest there is, and labelled as assumed. It is not quite "no plasma
//! binding": Rodgers and Rowland's equation still has the plasma's own lipid taking some up, and at
//! `f_u = 1` that makes its protein term slightly negative, 2% of p-xylene's liver `K_p`.
//!
//! # What it checks
//!
//! - **The descriptors** against counts made by hand from each dictionary entry. Aspirin is
//!   C9H8O4, with one O–H and four oxygens; its weight is 180.16.
//! - **Each verdict** against what the module's rules give for these inputs, worked by hand below.
//!   Benzene passes both screens with either force field. p-xylene's verdict is "not measured",
//!   because no binding was given, and so is aspirin's.
//! - **The two binding free energies give two very different `K_d`**, and the verdict cannot tell
//!   them apart. Benzene is small, so even the weaker binding is efficient per heavy atom. The
//!   factor between the two `K_d` is `exp(Δ/RT)`, checked against that closed form, and it is
//!   1 583, which is the number the module documentation quotes.
//!
//! # The thresholds are this example's choice
//!
//! Two violations to raise the absorption alert, and 0.3 kcal/mol per heavy atom for ligand
//! efficiency. Both are the usual values. Neither was checked against its paper here, and that is
//! why the library has no default for either. See [`pantometry::candidate::Thresholds`].

use pantometry::candidate::{
    assess, Assessment, BindingFreeEnergy, Call, Distribution, Evidence, Measured, Rule, Thresholds,
};
use pantometry::forcefield::{components, Component};
use pantometry::units::Temperature;

/// W3's and W4's temperature, and the experiments'.
const KELVIN: f64 = 298.15;

const THRESHOLDS: Thresholds = Thresholds {
    violations_for_alert: 2,
    ligand_efficiency: 0.3,
};

fn binding(kcal: f64, sigma: f64) -> Option<BindingFreeEnergy> {
    Some(BindingFreeEnergy::from_kcal_per_mol(
        Measured::new(kcal, sigma),
        Temperature::kelvin(KELVIN),
    ))
}

fn report(label: &str, a: &Assessment) {
    println!("{label}");
    print!("{a}");
    if matches!(a.distribution, Distribution::Reported(_)) {
        println!("  (K_p at f_u = 1, assumed: no fraction unbound was measured)");
    }
    println!();
}

fn main() {
    // The entries `pantometry-forcefield` ships, by their dictionary identifiers.
    let entry = |id: &str| Component::from_ccd(components::text(id).unwrap()).unwrap();
    let benzene = entry("BNZ");
    let p_xylene = entry("PXY");
    let aspirin = entry("AIN");

    let benzene_log_p = Some(Measured::new(1.92, 0.23));
    let uff = assess(
        &benzene,
        &Evidence {
            log_p: benzene_log_p,
            binding: binding(-8.975, 0.279),
            fraction_unbound: Some(1.0),
        },
        &THRESHOLDS,
    )
    .unwrap();
    let gaff = assess(
        &benzene,
        &Evidence {
            log_p: benzene_log_p,
            binding: binding(-4.61, 0.29),
            fraction_unbound: Some(1.0),
        },
        &THRESHOLDS,
    )
    .unwrap();
    let xylene = assess(
        &p_xylene,
        &Evidence {
            log_p: Some(Measured::new(3.06, 0.22)),
            binding: None,
            fraction_unbound: Some(1.0),
        },
        &THRESHOLDS,
    )
    .unwrap();
    let asa = assess(&aspirin, &Evidence::default(), &THRESHOLDS).unwrap();

    report("benzene, binding with UFF throughout", &uff);
    report("benzene, binding with GAFF, ff14SB and AM1-BCC", &gaff);
    report("p-xylene", &xylene);
    report("aspirin", &asa);

    // The descriptors, counted by hand from the dictionary entries.
    let counts = |a: &Assessment| {
        let d = a.descriptors;
        (d.donors, d.acceptors, d.heavy_atoms)
    };
    assert_eq!(counts(&uff), (0, 0, 6));
    assert_eq!(counts(&xylene), (0, 0, 8));
    assert_eq!(counts(&asa), (1, 4, 13));
    assert!((asa.descriptors.molecular_weight - 180.16).abs() <= 0.005);

    // Worked by hand. Benzene's log P interval is [1.46, 2.38], under 5, and every count is under
    // its limit, so absorption passes. Its ligand efficiency is 8.975/6 = 1.50 or 4.61/6 = 0.77,
    // each with an interval far above 0.3, so binding passes with either force field.
    for a in [&uff, &gaff] {
        assert_eq!(a.absorption.call, Call::Pass);
        assert_eq!(a.affinity.unwrap().call, Call::Pass);
        assert_eq!(a.verdict, Call::Pass);
    }
    // p-xylene's log P interval is [2.62, 3.50]: absorption passes, binding was not given.
    assert_eq!(xylene.absorption.call, Call::Pass);
    assert_eq!(xylene.verdict, Call::NotMeasured);
    // Aspirin's log P was not given, but one open limit cannot reach an alert at two.
    assert_eq!(asa.absorption.call, Call::Pass);
    assert_eq!(asa.absorption.open, 1);
    assert_eq!(asa.absorption.limit(Rule::LogP).call, Call::NotMeasured);
    assert_eq!(asa.verdict, Call::NotMeasured);
    assert_eq!(asa.distribution, Distribution::NotMeasured);
    // The hydrocarbons were given log P and f_u, and are neutral, so all three tissues are there.
    for a in [&uff, &gaff, &xylene] {
        assert!(matches!(a.distribution, Distribution::Reported(_)));
    }

    // The two force fields' K_d differ by exp((−4.61 + 8.975)/RT), with R written out here.
    let rt = 6.022_140_76e23 * 1.380_649e-23 / 4184.0 * KELVIN;
    let ratio =
        gaff.affinity.unwrap().dissociation_constant / uff.affinity.unwrap().dissociation_constant;
    let want = ((-4.61 + 8.975) / rt).exp();
    assert!(
        (ratio / want - 1.0).abs() <= 1e-12,
        "{ratio} against {want}"
    );
    // The factor the module documentation and EXAMPLES.md quote.
    assert_eq!(ratio.round(), 1583.0, "{ratio}");
    println!(
        "The two force fields' K_d differ by a factor of {ratio:.0}, and the verdict is the same."
    );
}
